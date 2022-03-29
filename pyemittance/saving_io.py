import numpy as np
from os import path, makedirs
from os.path import exists
import errno, os
import json
import datetime

from epics import caget, caget_many

this_dir, this_filename = path.split(__file__)
CONFIG_PATH = path.join(this_dir, "configs")
pv_savelist = json.load(open(CONFIG_PATH+'/save_scalar_pvs.json'))
savepaths = json.load(open(CONFIG_PATH+'/savepaths.json'))
opt_pv_info = json.load(open(CONFIG_PATH+'/opt_pv_info.json'))
meas_pv_info = json.load(open(CONFIG_PATH+'/meas_pv_info.json'))

meas_read_pv = meas_pv_info['meas_device']['pv']['read']
opt_pvs = opt_pv_info['opt_vars']


def mkdir_p(path):
    """Set up dirs for results in working dir"""
    try:
        makedirs(path)
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise

# Make directories if needed
mkdir_p(savepaths['images'])
mkdir_p(savepaths['summaries'])
mkdir_p(savepaths['fits'])
mkdir_p(savepaths['raw_saves'])

# Start headings
file_exists = path.exists(savepaths['summaries'] + "image_acq_quad_info.csv")

if not file_exists:

    # todo add others as inputs
    f = open(savepaths['summaries'] + "image_acq_quad_info.csv", "a+")
    f.write(
        f"{'timestamp'},{'ncol'},{'nrow'},{'roi_xmin'},{'roi_xmax'}"
        f",{'roi_ymin'},{'roi_ymax'},{'resolution'},{'bact'},"
        f"{'x_size'},{'y_size'},{'beamsizes[0]'},{'beamsizes[1]'},"
        f"{'beamsizes[2]'},{'beamsizes[3]'}\n")
    f.close()

file_exists = path.exists(savepaths['summaries'] + "beamsize_config_info.csv")

if not file_exists:
    # todo add others as inputs
    f = open(savepaths['summaries'] + "beamsize_config_info.csv", "a+")
    f.write(
        f"{'timestamp'},{'varx_cur'},{'vary_cur'},{'varz_cur'},"
        f"{'bact_cur'},{'xrms'},{'yrms'},{'xrms_err'},{'yrms_err'}\n")
    f.close()


def isotime():
    return datetime.datetime.utcnow().replace(tzinfo=datetime.timezone.utc).\
        astimezone().replace(microsecond=0).isoformat()

def save_image(im, ncol, nrow, timestamp,
              impath=savepaths['images'], avg_img=True):
    """Saves images with col,row info and corresp. settings"""

    if avg_img:

        np.save(str(impath) + f'img_avg_{timestamp}.npy', im)
        np.save(str(impath) + f'ncol_avg_{timestamp}.npy', ncol)
        np.save(str(impath) + f'nrow_avg_{timestamp}.npy', nrow)

    else:

        np.save(str(impath) + f'img_{timestamp}.npy', im)
        np.save(str(impath) + f'ncol_{timestamp}.npy', ncol)
        np.save(str(impath) + f'nrow_{timestamp}.npy', nrow)


def numpy_save(xrms, yrms, xrms_err, yrms_err, timestamp=False,
               savelist=pv_savelist['scalars'],
               path=savepaths['raw_saves']):
    ts = isotime()
    x = caget_many(savelist)
    x.append(ts)
    if timestamp:
        x.append(timestamp)
    else:
        x.append(ts)
    x.append(xrms)
    x.append(yrms)
    x.append(xrms_err)
    x.append(yrms_err)

    np.save(path + ts + '_x_.npy', np.array(x))


def save_config(xrms, yrms, xrms_err, yrms_err, timestamp,
                config_path=savepaths['summaries'],
                 im=None,
                impath=savepaths['images']):
    f = open(savepaths['summaries'] + "beamsize_config_info.csv", "a+")

    # todo make more general, pandas etc
    varx_cur = caget(opt_pvs[0])
    vary_cur = caget(opt_pvs[1])
    varz_cur = caget(opt_pvs[2])
    bact_cur = caget(meas_read_pv)
    f.write(f"{timestamp},{varx_cur},{vary_cur},{varz_cur},"
            f"{bact_cur},{xrms},{yrms},{xrms_err},{yrms_err}\n")

    f.close()

    if im:
        np.save((str(impath) + f'img_config_{timestamp}.npy', im.proc_image))


# def savesummary(beamsizes, timestamp=(datetime.datetime.now()).strftime("%Y-%m-%d_%H-%M-%S-%f")):
#     """Saves summary info for beamsize fits"""
#
#     # todo add others as inputs
#     f = open(savepaths['summaries'] + "image_acq_quad_info.csv", "a+")
#     bact = quad_read_pv.get()
#     x_size = x_size_pv.get()
#     y_size = y_size_pv.get()
#     f.write(
#         f"{timestamp},{ncol},{nrow},{roi_xmin},{roi_xmax},"
#         f"{roi_ymin},{roi_ymax},{resolution},{bact},{x_size},{y_size},"
#         f"{beamsizes[0]},{beamsizes[1]},{beamsizes[2]},{beamsizes[3]}\n")
#     f.close()