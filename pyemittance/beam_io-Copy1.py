""" This file constains all i/o to/from machine/sim and machine info """
import numpy as np
import sys, os, errno
import json
import time
import datetime
import matplotlib.pyplot as plt

from pyemittance.image import Image
from os.path import exists

from epics import caget, caput, PV
import epics

this_dir, this_filename = os.path.split(__file__)
DATA_PATH = os.path.join(this_dir, "configs")

        
def get_beamsizes_machine(quad_val):
    """Takes quad value as input,
     returns xrms, yrms, xrms_err, yrms_err"""
    setquad(quad_val)
    time.sleep(3)
    beamsize = get_beamsizes(use_profMon=False) 
    xrms = beamsize[0]
    yrms = beamsize[1]
    xrms_err = beamsize[2]
    yrms_err = beamsize[3]
    return xrms, yrms, xrms_err, yrms_err

def get_beam_image(subtract_bg = subtract_bg, post=None):
        """Get beam image from screen and return beamsize info and processed image"""

        if post:
            im = post[0]
            ncol = post[1]
            nrow = post[2]
        else:
            im = im_pv.get()
            ncol, nrow = n_col_pv.get(), n_row_pv.get()

        beam_image = Image(im, ncol, nrow, bg_image = bg_image)
        beam_image.reshape_im()

        if subtract_bg:
            beam_image.subtract_bg()

        if use_roi:

            beam_image.proc_image = beam_image.proc_image[roi_ymin:roi_ymax, roi_xmin:roi_xmax]

        beam_image.get_im_projection()

        # fit the profile and return the beamsizes
        beamsizes = beam_image.get_sizes(show_plots=False)

        #save_beam = list(np.array(beamsizes[0:4])*resolution/1e-6)
        
        timestamp = (datetime.datetime.now()).strftime("%Y-%m-%d_%H-%M-%S-%f") 
        #savesummary(beamsizes,timestamp)# pass beamsizes in um
        saveimage(im, ncol, nrow, timestamp, avg_img=False)
        
        return list(beamsizes) + [beam_image.proc_image]


def getbeamsizes_from_img(num_images = n_acquire, avg = avg_ims, subtract_bg = subtract_bg, post = None):
    """Returns xrms, yrms, xrms_err, yrms_err for multiple sampled images;
    can optionally average multiple images
    RETURNS IN RAW UNITS-- NEED TO MULTIPLY BY RESOLUTION FOR M"""
    
    xrms, yrms, xrms_err, yrms_err, xamp, yamp, im = [0]*num_images, [0]*num_images, [0]*num_images, [0]*num_images, [0]*num_images, [0]*num_images, [0]*num_images
    
    if post:
        ncol = post[0][1]
        nrow = post[0][2]
    else:
        ncol, nrow = n_col_pv.get(), n_row_pv.get()

    for i in range(0,num_images):
        
        
        repeat = True
        count = 0
        # retake bad images 3 times
       # print(i)
        
        while repeat:
            xrms[i], yrms[i], xrms_err[i], yrms_err[i], xamp[i], yamp[i], im[i] = get_beam_image(subtract_bg,post)
            
            #plt.imshow(im[i])
            #plt.show()

            
            count = count + 1 
            
            
            if xamp[i]>amp_threshold_x and yamp[i]>amp_threshold_y and xrms[i]>min_sigma and yrms[i]>min_sigma and xrms[i]<max_sigma and yrms[i]<max_sigma:
                # if conditions are met, stop resampling this image
                repeat = False
            elif count==3:
                # if still bad after 3 tries, return nan
                xrms[i], yrms[i], xrms_err[i], yrms_err[i], xamp[i], yamp[i] = np.nan,np.nan,np.nan,np.nan,np.nan,np.nan
                repeat = False
    
    #average images before taking fits
    if avg == True:
        

        im = np.mean(im, axis=0)
        
        im = Image(im, ncol, nrow, bg_image = bg_image)
        
        im.reshape_im()
        if subtract_bg:
            im.subtract_bg()
        im.get_im_projection()
        
        #plt.imshow(im.proc_image)
        
        mean_xrms, mean_yrms, mean_xrms_err, mean_yrms_err, mean_xamp, mean_yamp = im.get_sizes(show_plots=False)
        
        if mean_xamp<amp_threshold_x or mean_yamp<amp_threshold_y:
            return [np.nan, np.nan, np.nan, np.nan, np.nan, np.nan]
                
        #save_beam = list(np.array(beamsizes[0:4])*resolution/1e-6)
        
        timestamp = (datetime.datetime.now()).strftime("%Y-%m-%d_%H-%M-%S-%f") 
        saveimage(im.proc_image, ncol, nrow, timestamp, avg_img=True) # pass beamsizes in um
        
        return [mean_xrms, mean_yrms, mean_xrms_err, mean_yrms_err, mean_xamp, mean_yamp, im.proc_image]
        
        
        
    #average individual rms fits
    else: 
        idx = ~np.isnan(xrms)
        
        if True not in idx:
            return np.nan, np.nan, np.nan, np.nan, np.nan, np.nan

        mean_xrms = np.mean(np.array(xrms)[idx])
        mean_yrms = np.mean(np.array(yrms)[idx])
        mean_xrms_err = np.std(np.array(xrms)[idx])/np.sqrt(len(idx))
        mean_yrms_err = np.std(np.array(yrms)[idx])/np.sqrt(len(idx))
    #     mean_xrms_err = np.sqrt(np.mean(np.array(xrms_err)[idx]**2))
    #     mean_yrms_err = np.sqrt(np.mean(np.array(yrms_err)[idx]**2))
        mean_xamp = np.mean(np.array(xamp)[idx])
        mean_yamp = np.mean(np.array(yamp)[idx]) 
        
        im = Image(im[0], ncol, nrow, bg_image = bg_image)
        
        im.reshape_im()
        im.get_im_projection()
        
    
        return [mean_xrms, mean_yrms, mean_xrms_err, mean_yrms_err, mean_xamp, mean_yamp, im.proc_image]
        


def get_beamsizes(use_profMon=False, reject_bad_beam=True, save_summary = True, post = None):
    """Returns xrms, yrms, xrms_err, yrms_err, with options to reject bad beams and use either profmon or image processing"""
        
    xrms = np.nan
    yrms =  np.nan
    xrms_err =  np.nan
    yrms_err =  np.nan
    xamp =  np.nan
    yamp =  np.nan
    beamsizes = [np.nan,np.nan,np.nan,np.nan,np.nan,np.nan]
    im = None


    #remove hardcoding
    min_sigma=0.00000002
    max_sigma=0.1

    if reject_bad_beam:

        count = 0
        
        while xrms<=min_sigma or yrms<=min_sigma or xrms>max_sigma or yrms>max_sigma or np.isnan(np.array(beamsizes[0:6])).any():

            if count > 1:
                print("Low beam intensity/noisy or beam too small/large.")
                print("Waiting 1 sec and repeating measurement...")
                time.sleep(1)

            #make sure stats is checked on profmon gui
            if use_profMon:
                beamsizes=[]
                xrms, xrms_err = x_size_pv.get()*1e-6, 0 # in meters
                yrms, yrms_err = y_size_pv.get()*1e-6, 0 # in meters
                #print('bz',xrms,yrms)
                #print(x_size_pv)
                #print(y_size_pv)
                count = count + 1

            if not use_profMon:

                #if post:
                #    beamsizes = getbeamsizes_from_img(post = post)
                #:
                beamsizes = getbeamsizes_from_img()

                xrms = beamsizes[0]
                yrms = beamsizes[1]
                xrms_err = beamsizes[2]
                yrms_err = beamsizes[3]
                xamp = beamsizes[4]
                yamp = beamsizes[5] 
                #im = beamsizes[6]
                
                #print('raw px',xrms,yrms)              
                # convert to meters
                xrms = xrms*resolution 
                yrms = yrms*resolution 
                xrms_err = xrms_err*resolution
                yrms_err = yrms_err*resolution 
                
                print('bzs', beamsizes)
                

                #print('res',resolution)
                print('after', xrms, min_sigma, yrms, min_sigma, xrms, max_sigma, yrms, max_sigma)
                

                
                if count == 1:
                    # resample beamsize only 3 times
                    return np.nan, np.nan, np.nan, np.nan
                
                count = count + 1
                
                
    else:
            
            #make sure stats is checked on profmon gui
            if use_profMon:
                xrms, xrms_err = x_size_pv.get()*1e-6, 0 # in meters
                yrms, yrms_err = y_size_pv.get()*1e-6, 0 # in meters
                print("USING PROFMONNNNNNN!!!!")

            else:
                if post:
                    beamsizes = getbeamsizes_from_img(post = post)
                else:
                    beamsizes = getbeamsizes_from_img()
                #print(beamsizes)

                xrms = beamsizes[0]
                yrms = beamsizes[1]
                xrms_err = beamsizes[2]
                yrms_err = beamsizes[3]
                xamp = beamsizes[4]
                yamp = beamsizes[5]  
                #im = beamsizes[6]
                

                # convert to meters
                xrms = xrms*resolution 
                yrms = yrms*resolution 
                xrms_err = xrms_err*resolution
                yrms_err = yrms_err*resolution 


    if save_summary:
        timestamp=(datetime.datetime.now()).strftime("%Y-%m-%d_%H-%M-%S-%f")

        save_config(xrms,yrms,xrms_err,yrms_err,timestamp=timestamp)#,#im)
        numpy_save(xrms, yrms, xrms_err, yrms_err,timestamp)

    return xrms, yrms, xrms_err, yrms_err

def setinjector(set_list):
    sol_cntrl_pv.put(set_list[0])
    cq_cntrl_pv.put(set_list[1])
    sq_cntrl_pv.put(set_list[2])
    
def setquad(value):
    """Sets Q525 to new scan value"""
    meas_cntrl_pv.put(value)
    
def get_beamsize_inj(set_list, quad=meas_read_pv.get(), use_profMon=False):
    """Get beamsize fn that changes upstream cu injector and returns xrms and yrms in [m]"""   
    setinjector(set_list)
    setquad(quad)
    time.sleep(3)
    beamsize = get_beamsizes(use_profMon=use_profMon)
    return np.array([beamsize[0], beamsize[1]])

def quad_control(val=None, action="get"):
    if action == "get":
        return meas_read_pv.get()
    elif action == "set":
        setquad(val)
        return None
    else:
        raise ValueError("Invalid quad function")

def get_twiss0(filepath= DATA_PATH+'/beamline_info.json'):
    """Import Twiss0 from config file"""

    beamline_info = json.load(open(filepath))
    twiss0 = beamline_info['Twiss0']
    # emit, beta, alpha
    twiss0_by_dim = {'x': [twiss0[0], twiss0[2], twiss0[4]],
                     'y': [twiss0[1], twiss0[3], twiss0[5]]
                     }

    return twiss0_by_dim

def get_rmat(filepath=DATA_PATH+'/beamline_info.json'):
    """Import r-matrix from config file"""

    beamline_info = json.load(open(filepath))
    # if only separated by a drift:
    # rMat will be [1, L, 0, 1]
    rmat = beamline_info['rMat']
    # m11, m12, m21, m22
    rmat_array = np.array([ [rmat[0], rmat[1]],
                            [rmat[2], rmat[3]]
                            ])

    return rmat_array