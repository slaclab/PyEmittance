import numpy as np
import datetime
from epics import caget, PV
from pyemittance.image import Image
from pyemittance.saving_io import save_image, save_config
import os

import logging
logger = logging.getLogger(__name__)


def get_beamsizes_otrs(config_dict,
                       reject_bad_beam=True,
                       save_summary=False,
                         ):
    """Data acquisition from OTRS
    Additional option to reject bad beams
    Returns xrms, yrms, xrms_err, yrms_err
    """
    # Saving configs
    meas_pv_info = config_dict["meas_pv_info"]
    meas_read_pv = meas_pv_info["meas_device"]["pv"]["read"]
    savepaths = config_dict["savepaths"]
    # Load image processing setting info
    im_proc = config_dict["img_proc"]
    amp_threshold_x = im_proc["amp_threshold"]
    amp_threshold_y = im_proc["amp_threshold"]
    max_samples = im_proc["max_samples"]

    # Measurement PVs
    meas_pv_info = config_dict["meas_pv_info"]
    # in meters for emittance calc
    resolution = caget(meas_pv_info["diagnostic"]["pv"]["resolution"]) * 1e-6

    # Thresholds in meters
    min_sigma_meters = im_proc["min_sigma"] * resolution
    max_sigma_meters = im_proc["max_sigma"] * resolution

    xrms = np.nan
    yrms = np.nan
    xrms_err = np.nan
    yrms_err = np.nan
    xamp = np.nan
    yamp = np.nan
    beamsizes = [np.nan, np.nan, np.nan, np.nan, np.nan, np.nan]

    if reject_bad_beam:

        count = 0

        while (
            xrms <= min_sigma_meters
            or yrms <= min_sigma_meters
            or xrms > max_sigma_meters
            or yrms > max_sigma_meters
            or xamp < amp_threshold_x
            or yamp < amp_threshold_y
            or xamp * amp_threshold_x < 1500
            or yamp * amp_threshold_y < 1500
            or np.isnan(np.array(beamsizes[0:6])).any()
        ):


            if count == max_samples:
                # resample beamsize only max_samples times
                logger.info(f"Resampled {count-1} times (max_samples = {max_samples}, beam still out of bounds \n")

                logger.info(
                    f"xrms {xrms/1e-6:.2f} um, yrms {yrms/1e-6:.2f} um (threshold: min_rms {min_sigma_meters/1e-6:.2f} um, max_rms {max_sigma_meters/1e-6:.2f} um)"
                )
                logger.info(
                    f"xamp {xamp:.2f}, yamp {yamp:.2f} (amp_thresh: {amp_threshold_x}, in json)"
                )
                logger.info(
                    f"area_x {xamp*amp_threshold_x:.1f}, area_y {yamp*amp_threshold_y:.1f} (threshold: 1500, hardcoded)\n"
                )
                logger.info("Returning NaNs")

                return  {'xrms':    np.nan, 
                        'yrms':     np.nan,
                        'xrms_err': np.nan, 
                        'yrms_err': np.nan,
                           }

            if count > 0:
                logger.info("Low beam intensity/noisy or beam too small/large.")
                # logger.info("Waiting 1 sec and repeating measurement...")
                # time.sleep(1)         
                

            bdat = getbeamsizes_from_img(config_dict)
            
            # Extract beam sizes, convert to meters
            xrms = bdat['xrms'] * resolution
            yrms = bdat['yrms'] * resolution
            xrms_err = bdat['xrms_err'] * resolution
            yrms_err = bdat['yrms_err'] * resolution
            xamp = bdat['xamp']
            yamp = bdat['yamp']            

            # For legacy logic above
            beamsizes = [xrms, yrms, xrms_err, yrms_err, xamp, yamp]
                            
            count = count + 1

    else:

        bdat = getbeamsizes_from_img(config_dict)

    # Extract beam sizes, convert to meters
    xrms = bdat['xrms'] * resolution
    yrms = bdat['yrms'] * resolution
    xrms_err = bdat['xrms_err'] * resolution
    yrms_err = bdat['yrms_err'] * resolution
    xamp = bdat['xamp']
    yamp = bdat['yamp']

    if save_summary:
        timestamp = (datetime.datetime.now()).strftime("%Y-%m-%d_%H-%M-%S-%f")
        save_config(
            xrms,
            yrms,
            xrms_err,
            yrms_err,
            timestamp,
            meas_read_pv,
            configpath=savepaths["summaries"],
            impath=savepaths["images"],
        )

    return {'xrms': xrms, 
            'yrms': yrms,
            'xrms_err': xrms_err, 
            'yrms_err':yrms_err,
            'extra': bdat['image']
           }


def getbeamsizes_from_img(config_dict):
    """Returns xrms, yrms, xrms_err, yrms_err for multiple sampled images;
    can optionally average multiple images
    RETURNS IN RAW PIXEL UNITS-- NEED TO MULTIPLY BY RESOLUTION FOR METERS
    """
    # Load image processing setting info
    save_im_path = config_dict["savepaths"]["images"]
    im_proc = config_dict["img_proc"]
    subtract_bg = im_proc["subtract_bg"]
    bg_image = im_proc["background_im"]  # specify path to bg im in json file
    amp_threshold_x = im_proc["amp_threshold"]
    amp_threshold_y = im_proc["amp_threshold"]
    min_sigma = im_proc["min_sigma"]  # in raw units (pixels)
    max_sigma = im_proc["max_sigma"]
    max_samples = im_proc["max_samples"]
    avg_ims = im_proc["avg_ims"]
    num_images = im_proc["n_to_acquire"]

    # Measurement PVs
    meas_pv_info = config_dict["meas_pv_info"]
    # in meters for emittance calc
    n_col_pv = PV(meas_pv_info["diagnostic"]["pv"]["ncol"])
    n_row_pv = PV(meas_pv_info["diagnostic"]["pv"]["nrow"])

    xrms, yrms = [0] * num_images, [0] * num_images
    xrms_err, yrms_err = [0] * num_images, [0] * num_images
    xamp, yamp, im = [0] * num_images, [0] * num_images, [0] * num_images

    ncol, nrow = n_col_pv.get(), n_row_pv.get()
    
    for i in range(0, num_images):

        repeat = True
        count = 0

        # retake bad images
        while repeat:
            meas = get_beam_image(config_dict)
            xrms[i] = meas['xrms']
            yrms[i] = meas['yrms']
            xrms_err[i] = meas['xrms_err']
            yrms_err[i] = meas['yrms_err']
            xamp[i] = meas['xamp']
            yamp[i] = meas['yamp']
            im[i] = meas['proc_image']
            count = count + 1

            if (
                xamp[i] >= amp_threshold_x
                and yamp[i] >= amp_threshold_y
                and xrms[i] > min_sigma
                and yrms[i] > min_sigma
                and xrms[i] < max_sigma
                and yrms[i] < max_sigma
            ):
                # if conditions are met, stop resampling this image
                repeat = False
            elif count == max_samples:
                # if still bad after retries, return as is (reject in outer function)
                logger.info(
                    f"Beam params out of bounds in image {i} out of {num_images} samples"
                )
                repeat = False

    # average images before taking fits
    if avg_ims:

        idx = ~np.isnan(xrms)
        if True not in idx:
            logger.info(
                f"All {num_images} image NaNs are NaNs, trying to average images now."
            )
            all_nan = True
        else:
            all_nan = False

        im = np.mean(im, axis=0)
        im = Image(im, nrow, ncol, bg_image=bg_image)
        
        im.reshape_im()
        if subtract_bg:
            im.subtract_bg()
        im.get_im_projection()
        
        # Debug:
        #import matplotlib.pyplot as plt
        #plt.imshow(im.proc_image)

        meas = im.get_sizes(show_plots=False)
        mean_xrms = meas['xrms']
        mean_yrms = meas['yrms']
        mean_xrms_err = meas['xrms_err']
        mean_yrms_err = meas['yrms_err']
        mean_xamp = meas['xamp']
        mean_yamp = meas['yamp']

        if (
            mean_xamp < amp_threshold_x
            or mean_yamp < amp_threshold_y
            or mean_xrms < min_sigma
            and mean_yrms < min_sigma
            or mean_xrms > max_sigma
            and mean_yrms > max_sigma
        ):
            if not all_nan:
                logger.info(f"Beam params out of bounds in averaged image")
            else:
                logger.info(
                    f"Beam params out of bounds in averaged image, initial {num_images} all NaNs"
                )
                # return [np.nan, np.nan, np.nan, np.nan, np.nan, np.nan]

        # save_beam = list(np.array(beamsizes[0:4])*resolution/1e-6)

        timestamp = (datetime.datetime.now()).strftime("%Y-%m-%d_%H-%M-%S-%f")
        # pass beamsizes in um
        if save_im_path is None:
            pass
        elif os.path.exists(save_im_path):
            save_image(im.proc_image,  nrow, ncol, timestamp, impath=save_im_path, avg_img=True)
        else:
            logger.warning(f"Save image path does not exist: {save_im_path}, not saving")


    # average individual rms fits
    else:
        idx = ~np.isnan(xrms)

        if True not in idx:
            logger.info("All points are NaNs")
            return np.nan, np.nan, np.nan, np.nan, np.nan, np.nan

        mean_xrms = np.mean(np.array(xrms)[idx])
        mean_yrms = np.mean(np.array(yrms)[idx])
        mean_xrms_err = np.std(np.array(xrms)[idx]) / np.sqrt(len(idx))
        mean_yrms_err = np.std(np.array(yrms)[idx]) / np.sqrt(len(idx))
        mean_xamp = np.mean(np.array(xamp)[idx])
        mean_yamp = np.mean(np.array(yamp)[idx])

        im = Image(im[0], nrow, ncol, bg_image=bg_image)

        im.reshape_im()
        im.get_im_projection()

        
        
    return {
            'xrms':     mean_xrms,
            'yrms':     mean_yrms,
            'xrms_err': mean_xrms_err,
            'yrms_err': mean_yrms_err,
            'xamp':     mean_xamp,
            'yamp':     mean_yamp,
            'image':         im,
        }


def get_beam_image(config_dict):
    """Get beam image from screen
    Return beamsize info and processed image
    """
    # Load image processing setting info
    save_im_path = config_dict["savepaths"]["images"]
    im_proc = config_dict["img_proc"]
    subtract_bg = im_proc["subtract_bg"]
    bg_image = im_proc["background_im"]  # specify path to bg im in json file
    use_roi = im_proc["use_roi"]
    roi_xmin = im_proc["roi"]["xmin"]
    roi_ymin = im_proc["roi"]["ymin"]
    roi_xmax = im_proc["roi"]["xmax"]
    roi_ymax = im_proc["roi"]["ymax"]

    # Measurement PVs
    meas_pv_info = config_dict["meas_pv_info"]
    # in meters for emittance calc
    im_pv = PV(meas_pv_info["diagnostic"]["pv"]["image"])
    n_col_pv = PV(meas_pv_info["diagnostic"]["pv"]["ncol"])
    n_row_pv = PV(meas_pv_info["diagnostic"]["pv"]["nrow"])

    im = im_pv.get()
    ncol, nrow = n_col_pv.get(), n_row_pv.get()

    beam_image = Image(im, nrow, ncol, bg_image=bg_image)
    beam_image.reshape_im()

    if subtract_bg:
        beam_image.subtract_bg()

    if use_roi:
        beam_image.proc_image = beam_image.proc_image[
            roi_ymin:roi_ymax, roi_xmin:roi_xmax
        ]

    beam_image.get_im_projection()

    # fit the profile and return the beamsizes
    beamsize_data = beam_image.get_sizes(show_plots=False)

    # save_beam = list(np.array(beamsizes[0:4])*resolution/1e-6)

    timestamp = (datetime.datetime.now()).strftime("%Y-%m-%d_%H-%M-%S-%f")
    # savesummary(beamsizes,timestamp)# pass beamsizes in um
    if save_im_path is None:
        pass
    elif os.path.exists(save_im_path):
        save_image(im.proc_image,  nrow, ncol, timestamp, impath=save_im_path, avg_img=True)
    else:
        logger.warning(f"Save image path does not exist: {save_im_path}, not saving")
    logger.info(timestamp)

    #return list(beamsizes) + [beam_image.proc_image]

    # Attach image
    beamsize_data['proc_image'] = beam_image.proc_image
    return beamsize_data