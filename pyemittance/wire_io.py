import time
import datetime
from epics import PV
from pyemittance.saving_io import numpy_save, save_config


def get_beamsizes_wire(online=False, 
                       config_dict=None, 
                       save_summary=False,
                       multiwire=False
                      ):
    """Main function imported by beam_io
    Returns xrms, yrms, xrms_err, yrms_err
    """
    # TODO: check which of configs below are actually needed for wirescans
    # Saving configs
    opt_pv_info = config_dict['opt_pv_info']
    meas_pv_info = config_dict['meas_pv_info']
    meas_read_pv = meas_pv_info['meas_device']['pv']['read']
    opt_pvs = opt_pv_info['opt_vars']
    savepaths = config_dict['savepaths']
    pv_savelist = config_dict['save_scalar_pvs']
    
    # Measurement PVs
    meas_pv_info = config_dict['meas_pv_info']
    
    # BS in meters for emittance calc
    if not multiwire:
        scan_pv = [PV(meas_pv_info['diagnostic']['pv']['scan'])]
        x_size_pv = [PV(meas_pv_info['diagnostic']['pv']['xsize'])]
        y_size_pv = [PV(meas_pv_info['diagnostic']['pv']['ysize'])]
    else:
        # lists of PVs
        scan_pv = [PV(wire_pv) for wire_pv in meas_pv_info['diagnostic']['pv']['scan']]
        x_size_pv = [PV(x_pv) for x_pv in meas_pv_info['diagnostic']['pv']['xsize']]
        y_size_pv = [PV(y_pv) for y_pv in meas_pv_info['diagnostic']['pv']['ysize']]

    # run wire scans
    # if multiwire, runs them back to back
    error_dict = get_beamsize(online=online, scan_pv=scan_pv)

    # Get list of results
    xrms, yrms = [], []
    xrms_err, yrms_err = [], []
    for x_pv, y_pv in zip(x_size_pv, y_size_pv):
        # read in PVs, returned in um
        xrms.append(x_pv.get()*1e-6)
        yrms.append(y_pv.get()*1e-6)
        # add some error estimate, no PV that returns this as of 11/2022
        xrms_err.append(xrms*0.02)
        yrms_err.append(yrms*0.02)
    
    # TODO: this saves scalars, update to saving lists better
    if save_summary:
        timestamp = (datetime.datetime.now()).strftime("%Y-%m-%d_%H-%M-%S-%f")
        save_config(xrms,
                    yrms,
                    xrms_err,
                    yrms_err,
                    timestamp,
                    meas_read_pv,
                    opt_pvs,
                    configpath=savepaths['summaries'],
                    impath=savepaths['images']
                    )
                    
        numpy_save(xrms,
                   yrms,
                   xrms_err,
                   yrms_err,
                   timestamp,
                   savelist=pv_savelist['scalars'],
                   path=savepaths['raw_saves']
                   )
    
    if multiwire:
        # Check that we have enough data 
        len_data = sum(error_dict.values())
        if len_data < 3:
            print("Less than 3 data points successful.")
            print("Do not trust results.")

    # These are now lists
    # TODO: pass error dict to output_dict and flag there
    return xrms, yrms, xrms_err, yrms_err


def get_beamsize(online, scan_pv):
    error = {}
    for wire_pv in scan_pv:
        # Track faults
        error[f"{wire_pv}"] = False
        
        # When measurement is done, associated
        # WIRE:*:XRMS and *:YRMS PVs are updated.
        if online:
            if scan_pv.get() != 0:
                raise NotImplementedError(f"WS {scan_pv} not ready for running.")
            scan_pv.put(1)
            time.sleep(1)  

            status = scan_pv.get()
            if status == 2:
                while scan_pv.get()!= 0:
                    time.sleep(5) 
                time.sleep(3)  # to not break the wire scanner
            elif status == 0:
                print(f"WS {wire_pv} acquired successfully.")
            else: 
                print(f"WS {wire_pv} did not run. Status {status}.")
                error[f"{wire_pv}"] = True
                
    return error
                