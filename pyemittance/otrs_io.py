import numpy as np
from os import path
import json
import time
import datetime

from epics import caget, PV
from pyemittance.image import Image
from pyemittance.saving_io import save_image, numpy_save, save_config


this_dir, this_filename = path.split(__file__)
CONFIG_PATH = path.join(this_dir, "configs/LCLS2_OTR3")

# Load image processing setting info
im_proc = json.load(open(CONFIG_PATH+'/img_proc.json'))

subtract_bg = im_proc['subtract_bg']
bg_image = im_proc['background_im']  # specify path to bg im in json file
use_roi = im_proc['use_roi']
roi_xmin = im_proc['roi']['xmin']
roi_ymin = im_proc['roi']['ymin']
roi_xmax = im_proc['roi']['xmax']
roi_ymax = im_proc['roi']['ymax']
avg_ims = im_proc['avg_ims']
n_acquire = im_proc['n_to_acquire']
amp_threshold_x = im_proc['amp_threshold']
amp_threshold_y = im_proc['amp_threshold']
min_sigma = im_proc['min_sigma'] # in raw units (pixels)
max_sigma = im_proc['max_sigma']
max_samples = im_proc['max_samples']

# Measurement PVs
meas_pv_info = json.load(open(CONFIG_PATH + '/meas_pv_info.json'))

# in meters for emittance calc
resolution = caget(meas_pv_info['diagnostic']['pv']['resolution']) *1e-6
im_pv = PV(meas_pv_info['diagnostic']['pv']['image'])
n_col_pv = PV(meas_pv_info['diagnostic']['pv']['ncol'])
n_row_pv = PV(meas_pv_info['diagnostic']['pv']['nrow'])
x_size_pv = PV(meas_pv_info['diagnostic']['pv']['profmonxsize'])
y_size_pv = PV(meas_pv_info['diagnostic']['pv']['profmonysize'])

# thresholds in meters
min_sigma_meters = im_proc['min_sigma']*resolution
max_sigma_meters = im_proc['max_sigma']*resolution

def get_beamsizes_otrs(use_profmon):
    """Main function imported by beam_io
    Option to use ProfMon PVs to get beamsizes
    Returns xrms, yrms, xrms_err, yrms_err
    """
    beamsize = get_beamsizes(use_profMon=use_profmon)
    xrms = beamsize[0]
    yrms = beamsize[1]
    xrms_err = beamsize[2]
    yrms_err = beamsize[3]
    return xrms, yrms, xrms_err, yrms_err

def get_beamsizes(use_profMon=False, reject_bad_beam=True,
                  save_summary=False, post=None):
    """ Data acquisition from OTRS
    Option to use either ProfMon or image processing
    Additional option to reject bad beams
    Returns xrms, yrms, xrms_err, yrms_err
    """
    #use_profMon = True

    xrms = np.nan
    yrms = np.nan
    xrms_err = np.nan
    yrms_err = np.nan
    xamp = np.nan
    yamp = np.nan
    beamsizes = [np.nan, np.nan, np.nan, np.nan, np.nan, np.nan]
    im = None

    if reject_bad_beam:

        count = 0

        while xrms <= min_sigma_meters or yrms <= min_sigma_meters \
                or xrms > max_sigma_meters or yrms > max_sigma_meters \
                or xamp < amp_threshold_x or yamp < amp_threshold_y \
                or xamp*amp_threshold_x<1500 or yamp*amp_threshold_y<1500 \
                or np.isnan(np.array(beamsizes[0:6])).any():
                
            if count == max_samples:
                # resample beamsize only max_samples times
                print(f"Resampled {count-1} times, beam still out of bounds \n")
                print(f"xrms {xrms/1e-6:.2f} um, yrms {yrms/1e-6:.2f} um (threshold: min_rms {min_sigma_meters/1e-6:.2f} um, max_rms {max_sigma_meters/1e-6:.2f} um)")
                print(f"xamp {xamp:.2f}, yamp {yamp:.2f} (amp_thresh: {amp_threshold_x}, in json)")
                print(f"area_x {xamp*amp_threshold_x:.1f}, area_y {yamp*amp_threshold_y:.1f} (threshold: 1500, hardcoded)\n")
                print("Returning NaNs")
                return np.nan, np.nan, np.nan, np.nan
            
            if count > 0:
                print("Low beam intensity/noisy or beam too small/large.")
                print("Waiting 1 sec and repeating measurement...")
                time.sleep(1)

            # if this fails, make sure stats is checked on profmon gui
            if use_profMon:
                beamsizes = []
                xrms = x_size_pv.get() * 1e-6  # in meters
                yrms = y_size_pv.get() * 1e-6  # in meters
                # add 2% error on ProfMon measurement
                xrms_err = xrms*0.02
                yrms_err = yrms*0.02
                
                xamp, yamp = amp_threshold_x, amp_threshold_y # DOES NOT UPDATE W/PROFMON

                count = count + 1

            if not use_profMon:

                # if post:
                #    beamsizes = getbeamsizes_from_img(post = post)
                #:
                beamsizes = getbeamsizes_from_img()

                xrms = beamsizes[0]
                yrms = beamsizes[1]
                xrms_err = beamsizes[2]
                yrms_err = beamsizes[3]
                xamp = beamsizes[4]
                yamp = beamsizes[5]
                # im = beamsizes[6]

                # convert to meters
                xrms = xrms * resolution
                yrms = yrms * resolution
                xrms_err = xrms_err * resolution
                yrms_err = yrms_err * resolution

            count = count + 1


    else:

        # make sure stats is checked on profmon gui
        if use_profMon:
            xrms, xrms_err = x_size_pv.get() * 1e-6, 0  # in meters
            yrms, yrms_err = y_size_pv.get() * 1e-6, 0  # in meters
            xamp, yamp = amp_threshold_x, amp_threshold_y # DOES NOT UPDATE W/PROFMON
            #print("Using ProfMon beamsizes.")

        else:
            if post:
                beamsizes = getbeamsizes_from_img(post=post)
            else:
                beamsizes = getbeamsizes_from_img()

            xrms = beamsizes[0]
            yrms = beamsizes[1]
            xrms_err = beamsizes[2]
            yrms_err = beamsizes[3]
            xamp = beamsizes[4]
            yamp = beamsizes[5]
            # im = beamsizes[6]

            # convert to meters
            xrms = xrms * resolution
            yrms = yrms * resolution
            xrms_err = xrms_err * resolution
            yrms_err = yrms_err * resolution

    if save_summary:
        timestamp = (datetime.datetime.now()).strftime("%Y-%m-%d_%H-%M-%S-%f")
        save_config(xrms, yrms, xrms_err, yrms_err, timestamp)
        numpy_save(xrms, yrms, xrms_err, yrms_err, timestamp)

    return xrms, yrms, xrms_err, yrms_err

def getbeamsizes_from_img(num_images=n_acquire, avg=avg_ims,
                          subtract_bg=subtract_bg, post=None):
    """Returns xrms, yrms, xrms_err, yrms_err for multiple sampled images;
    can optionally average multiple images
    RETURNS IN RAW PIXEL UNITS-- NEED TO MULTIPLY BY RESOLUTION FOR METERS
    """

    xrms, yrms = [0] * num_images, [0] * num_images
    xrms_err, yrms_err = [0] * num_images, [0] * num_images
    xamp, yamp, im = [0] * num_images, [0] * num_images, [0] * num_images

    if post:
        ncol = post[0][1]
        nrow = post[0][2]
    else:
        ncol, nrow = n_col_pv.get(), n_row_pv.get()

    for i in range(0, num_images):
        
        repeat = True
        count = 0

        # retake bad images
        while repeat:
            meas = get_beam_image(subtract_bg, post)
            xrms[i], yrms[i] = meas[0:2]
            xrms_err[i], yrms_err[i] = meas[2:4]
            xamp[i], yamp[i], im[i] = meas[4:]

            # plt.imshow(im[i])
            # plt.show()

            count = count + 1

            if xamp[i] >= amp_threshold_x and yamp[i] >= amp_threshold_y \
                    and xrms[i] > min_sigma and yrms[i] > min_sigma \
                    and xrms[i] < max_sigma and yrms[i] < max_sigma:
                # if conditions are met, stop resampling this image
                repeat = False
            elif count == max_samples:
                # if still bad after retries, return as is (reject in outer function)
                print(f"Beam params out of bounds in image {i} out of {num_images} samples")
                # xrms[i], yrms[i], xrms_err[i], yrms_err[i], xamp[i], yamp[
                #     i] = np.nan, np.nan, np.nan, np.nan, np.nan, np.nan
                repeat = False

    # average images before taking fits
    if avg == True:
                      
        idx = ~np.isnan(xrms)
        if True not in idx:
            print(f"All {num_images} image NaNs are NaNs, trying to average images now.")
            all_nan = True
        else:
            all_nan = False

        im = np.mean(im, axis=0)

        im = Image(im, ncol, nrow, bg_image=bg_image)

        im.reshape_im()
        if subtract_bg:
            im.subtract_bg()
        im.get_im_projection()

        # plt.imshow(im.proc_image)

        meas =  im.get_sizes(show_plots=False)
        mean_xrms, mean_yrms = meas[0:2]
        mean_xrms_err, mean_yrms_err = meas[2:4]
        mean_xamp, mean_yamp = meas[4:]

        if mean_xamp < amp_threshold_x or mean_yamp < amp_threshold_y \
        or mean_xrms < min_sigma and mean_yrms < min_sigma \
        or mean_xrms > max_sigma and mean_yrms > max_sigma:
            if not all_nan:
                print(f"Beam params out of bounds in averaged image")
            else:
                print(f"Beam params out of bounds in averaged image, initial {num_images} all NaNs")
                #return [np.nan, np.nan, np.nan, np.nan, np.nan, np.nan]

        # save_beam = list(np.array(beamsizes[0:4])*resolution/1e-6)

        timestamp = (datetime.datetime.now()).strftime("%Y-%m-%d_%H-%M-%S-%f")
        # pass beamsizes in um
        save_image(im.proc_image, ncol, nrow, timestamp, avg_img=True)

        return [mean_xrms, mean_yrms,
                mean_xrms_err, mean_yrms_err,
                mean_xamp, mean_yamp, im.proc_image]


    # average individual rms fits
    else:
        idx = ~np.isnan(xrms)

        if True not in idx:
            print("All points are NaNs")
            return np.nan, np.nan, np.nan, np.nan, np.nan, np.nan

        mean_xrms = np.mean(np.array(xrms)[idx])
        mean_yrms = np.mean(np.array(yrms)[idx])
        mean_xrms_err = np.std(np.array(xrms)[idx]) / np.sqrt(len(idx))
        mean_yrms_err = np.std(np.array(yrms)[idx]) / np.sqrt(len(idx))
        mean_xamp = np.mean(np.array(xamp)[idx])
        mean_yamp = np.mean(np.array(yamp)[idx])

        im = Image(im[0], ncol, nrow, bg_image=bg_image)

        im.reshape_im()
        im.get_im_projection()

        return [mean_xrms, mean_yrms,
                mean_xrms_err, mean_yrms_err,
                mean_xamp, mean_yamp, im.proc_image]

def get_beam_image(subtract_bg=subtract_bg, post=None):
    """Get beam image from screen
    Return beamsize info and processed image
    """

    if post:
        im = post[0]
        ncol = post[1]
        nrow = post[2]
    else:
        im = im_pv.get()
        ncol, nrow = n_col_pv.get(), n_row_pv.get()

    beam_image = Image(im, ncol, nrow, bg_image=bg_image)
    beam_image.reshape_im()

    if subtract_bg:
        beam_image.subtract_bg()

    if use_roi:
        beam_image.proc_image = beam_image.proc_image[roi_ymin:roi_ymax, roi_xmin:roi_xmax]

    beam_image.get_im_projection()

    # fit the profile and return the beamsizes
    beamsizes = beam_image.get_sizes(show_plots=False)

    # save_beam = list(np.array(beamsizes[0:4])*resolution/1e-6)

    timestamp = (datetime.datetime.now()).strftime("%Y-%m-%d_%H-%M-%S-%f")
    # savesummary(beamsizes,timestamp)# pass beamsizes in um
    save_image(im, ncol, nrow, timestamp, avg_img=False)

    return list(beamsizes) + [beam_image.proc_image]