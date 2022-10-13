import time
import datetime
from epics import PV
from pyemittance.saving_io import numpy_save, save_config


def get_beamsizes_wire(online=False, config_dict=None, save_summary=False):
    """Main function imported by beam_io
    Returns xrms, yrms, xrms_err, yrms_err
    """
    # Measurement PVs
    meas_pv_info = config_dict['meas_pv_info']
    # in meters for emittance calc
    scan_pv = PV(meas_pv_info['diagnostic']['pv']['scan'])
    x_size_pv = PV(meas_pv_info['diagnostic']['pv']['xsize'])
    y_size_pv = PV(meas_pv_info['diagnostic']['pv']['ysize'])

    # run wire scans
    get_beamsize(online=online, scan_pv=scan_pv)

    # read in PVs 
    xrms = x_size_pv.get()*1e-6
    yrms = y_size_pv.get()*1e-6
    # add some error estimate
    xrms_err = xrms*0.02
    yrms_err = yrms*0.02
    
    if save_summary:
        timestamp = (datetime.datetime.now()).strftime("%Y-%m-%d_%H-%M-%S-%f")
        save_config(xrms, yrms, xrms_err, yrms_err, timestamp)
        numpy_save(xrms, yrms, xrms_err, yrms_err, timestamp)
    
    return xrms, yrms, xrms_err, yrms_err


def get_beamsize(online, scan_pv):
    if online:
        scan_pv.put(1)
        time.sleep(1)  
        
        status = scan_pv.get()
        if status == 2:
            while scan_pv.get()!= 0:
                time.sleep(5) 
            time.sleep(3) # to not break the wire scanner
             
        else: 
            print(f"WS did not run. Status {status}.")