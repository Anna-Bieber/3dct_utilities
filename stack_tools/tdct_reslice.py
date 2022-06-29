# -*- coding: utf-8 -*-
"""
Created on Sun Feb 16 13:08:47 2020

@author: Anna
"""
import os
import sys
import time
import logging
import numpy as np
from scipy import interpolate

import tifffile as tf


log = logging.getLogger("Reslice")

def reslice(img, step_in, step_out, interpolationmethod='linear', 
         fname_out='resliced.tif'):
    """Main function handling the file type and parsing of filenames/directories"""

    if len(img.shape) < 3:
        log.error( "ERROR: This seems to be a 2D image with the shape {0}. Please select a stack image file.".format(img.shape))
        return


    ## Start Processing
    log.debug( "Interpolating...")
    # Actual interpolation
    img_int = interpol(img, step_in, step_out, interpolationmethod)

    if type(img_int) == str:
        log.debug( img_int)
        return
    if img_int is not None:
        log.debug( "Saving interpolated stack as: {}".format(fname_out) )
        tf.imsave(fname_out, img_int)            
        log.debug( "        ...done.")



def interpol(img, step_in, step_out, interpolationmethod):
    """Main function for interpolating image stacks via polyfit"""
    ## Depending on tiff format the file can have different shapes; e.g. z,y,x or c,z,y,x
    if len(img.shape) == 4 and img.shape[0] == 1:
        img = np.squeeze(img, axis=0)
    elif len(img.shape) == 4 and img.shape[0] > 1:
        return "ERROR: I'm sorry, I cannot handle multichannel files: "+str(img.shape)

    if len(img.shape) == 3:
        ## Number of slices in original stack
        sl_in = img.shape[0]
        ## Number of slices in interpolated stack
        # Discarding last data point. e.g. 56 in i.e.
        # 55 steps * (309 nm original spacing / 161.25 nm new spacing) = 105.39 -> int() = 105 + 1 = 106
        sl_out = int((sl_in-1)*(step_in/step_out)) + 1
        ## Interpolate image stack shape
        img_int_shape = (sl_out, img.shape[1], img.shape[2])
    else:
        return "ERROR: I only know tiff stack image formats in z,y,x or c,z,y,x with one channel"

    if interpolationmethod == 'none':
        return None
    elif interpolationmethod == 'linear':
        log.debug( "Nr. of slices: {} in, {} out): ".format(sl_in, sl_out) )
        return linear(img, img_int_shape, step_in, step_out, sl_in, sl_out)
    elif interpolationmethod == 'spline':
        log.debug( "Nr. of slices: {} in, {} out): ".format(sl_in, sl_out) )
        return spline(img, img_int_shape, step_in, step_out, sl_in, sl_out)
    else:
        return "Please specify the interpolation method ('linear', 'spline', 'none')."

def spline(img, img_int_shape, step_in, step_out, sl_in, sl_out):
    """
    Spline interpolation

    # possible depricated due to changes in code -> marked for futur code changes
    step_in : step size input stack
    step_out : step size output stack
    sl_in : slices input stack
    sl_out : slices output stack
    """
    ## Known x values in interpolated stack size.
    zx = np.arange(0,sl_out,step_in/step_out)
    zxnew = np.arange(0, (sl_in-1)*step_in/step_out, 1)  # First slice of original and interpolated are both 0. n-1 to discard last slice
    if step_in/step_out < 1.0:
        zx_mod = []
        for i in range(img.shape[0]):
            zx_mod.append(zx[i])
        zx = zx_mod

    ## Create new numpy array for the interpolated image stack
    img_int = np.zeros(img_int_shape,img.dtype)
    log.debug( "Interpolated stack shape: {}".format(img_int.shape) )

    r_sl_out = list(range(sl_out))

    ping = time.time()
    for px in range(img.shape[-1]):
        for py in range(img.shape[-2]):
            spl = interpolate.InterpolatedUnivariateSpline(zx, img[:,py,px])
            np.put(img_int[:,py,px], r_sl_out, spl(zxnew))
        sys.stdout.write("\r%d%%" % int(px*100/img.shape[-1]))
        sys.stdout.flush()
    pong = time.time()
    log.debug( "This interpolation took {0} seconds".format(pong - ping))
    return img_int


def linear(img, img_int_shape, step_in, step_out, sl_in, sl_out):
    """Linear interpolation"""
    ##  Determine interpolated slice positions
    sl_int = np.arange(0,sl_in-1,step_out/step_in)  # sl_in-1 because last slice is discarded (no extrapolation)

    ## Create new numpy array for the interpolated image stack
    img_int = np.zeros(img_int_shape,img.dtype)
    log.debug( "Interpolated stack shape: {} ".format(img_int.shape) )

    ## Calculate distances from every interpolated image to its next original image
    sl_counter = 0
    ping = time.time()
    for i in sl_int:
        int_i = int(i)
        lower = i-int_i
        upper = 1-(lower)
        img_int[sl_counter,:,:] = img[int_i,:,:]*upper + img[int_i+1,:,:]*lower
        sl_counter += 1
    pong = time.time()
    log.debug( "This interpolation took {0} seconds".format(pong - ping))
    return img_int