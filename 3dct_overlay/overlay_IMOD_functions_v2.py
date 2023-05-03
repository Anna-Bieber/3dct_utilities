# -*- coding: utf-8 -*-
"""
Created on Thu May  7 18:15:06 2020
Changelog

2023/02: Fix to make functions work for 16bit fluorescence images. 
        Main changes in plot_overlay_IMOD and make_overlay_data_multifluo
@author: Anna
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, rgb2hex

import os
import tifffile as tf
import mrcfile
import subprocess
from scipy.interpolate import interpn
from scipy.spatial.transform import Rotation

from itertools import cycle
from pathlib import Path



def cmap_fading(base_color):
    """Creates a fading colormap from base_color. Base_color is given as list or tuple of rgb values,
    with values between 0 and 1. Used for plotting with matplotlib."""
    # scale base color to 1 
    if np.max(base_color) > 1:
        c = [val/255 for val in base_color] # adjust range
    else:
        c = list(base_color)
    colors = np.zeros((100000, 4), np.float64)
    for i, val in enumerate(c):
        colors[:,i] = val # set RGB to base color
    colors[:,3] = np.linspace(0, 1, colors.shape[0]) # make alpha gradient
    name = 'fading_{}'.format(rgb2hex(c))
    cmap = ListedColormap(colors, name=name)
    plt.cm.register_cmap(name=name, cmap=cmap)
    
    return cmap

def lut_imagej(base_color, dtype=np.uint8):
    """Create a look-up table for saving files for ImageJ. base_color is given as 
    [r,g,b] with values ranging from 0 to 1. Changelog 2023/02: account for dtype. However,
    in practice it doesn't seem to matter if a 256-color LUT is given to uint16 data."""
    # scale base color to 1. base color should be (1,1,1) for gray LUT!
    if np.max(base_color) > 1:
        c = [val/255 for val in base_color] # adjust range
    else:
        c = list(base_color) 
    # Intensity value range
    if dtype==np.uint8:
        max_val = 255
    else:
        max_val = 65535
    
    val_range = np.arange(max_val+1, dtype=dtype)
    # set up LUT
    lut = np.zeros((3, max_val+1), dtype=dtype)
    for i, val in enumerate(c): # go through RGB
        lut[i,:] = val*val_range
    
    return lut

# Function for reading 3DCT results    
def read_3dct_results(results_file_name, verbose=True):
    """
    Function that reads parameters from 3dct correlation txt file
    saves parameters into a dictionary with parameter names as keys
    """
    d={}
    with open(results_file_name, 'r') as results_fd:
        for line in results_fd:
            if line.startswith('#   - rotation'):
                d['phi'], d['psi'], d['theta'] = np.fromstring(line.split(': [')[1].split(']')[0], sep=',')
                continue
            if line.startswith('#   - scale'):
                d['scale'] = float(line.split('= ')[1])
                continue
            if line.startswith('#   - translation for rotation around [0,0,0]'):
                d['trans'] = np.fromstring(line.split('= [')[1].split(']')[0], sep=',')
                continue
            if line.startswith('#   - translation for rotation around'):
                d['origin'] = np.fromstring(line.split('[')[1].split(']')[0], sep=',')
                d['trans_or'] = np.fromstring(line.split('= [')[1].split(']')[0], sep=',')
                break
    if verbose:
        print("Please check that the following parameters are the same as those in the 3DCT results file:")
        print("phi, theta, psi = {}".format((d['phi'], d['theta'], d['psi'])))
        print("scale = {}".format(d['scale']))
        print("translation around [0, 0, 0] = {}".format(d['trans']))
        print("translation around {} = {}".format(d['origin'], d['trans_or']))
    
    return d

# Function for defining extent of fluo image after rotation
def extent_after_rotation(orig_image_size, rot_image_size, 
                          rotation_matrix, scale, translation, 
                          round_values=True, verbose=True):
    """
    Input:
            Sizes of original stack and rotated stack in order x,y,z (flip if read by tifffile)
            rotation matrix, scale and translation from 3DCT output
            round_values: should output values be rounded to integers?
    """
    fluo_center_orig = np.dot(orig_image_size, 0.5) # Original image center
    # Location of image center after rotation around 0,0,0 & scaling & translation
    fluo_center_rot0 = (scale * np.dot(rotation_matrix,fluo_center_orig.T)+translation) # do normal rotate_array translation on these coordinates
    # Image center after rotatevol
    fluo_center_rotatevol = np.dot(rot_image_size, 0.5) 
    
    if verbose:
        print("Original image stack center: {}".format(fluo_center_orig))
        print("Image stack center after rotation around 0,0,0: {}".format(fluo_center_rot0))
        print("Image stack center after rotatevol: {}".format(fluo_center_rotatevol))
        
    # Calculate extent for plotting scaled MIP on IB image
    # fluo_center_rotatevol is at the same time the distance to the edges in x and y
    xmin = fluo_center_rot0[0] - scale*fluo_center_rotatevol[0]
    xmax = fluo_center_rot0[0] + scale*fluo_center_rotatevol[0]
    
    ymin = fluo_center_rot0[1] - scale*fluo_center_rotatevol[1]
    ymax = fluo_center_rot0[1] + scale*fluo_center_rotatevol[1]
    
    if round_values:
        xmin = np.rint(xmin)
        xmax = np.rint(xmax)
        ymin = np.rint(ymin)
        ymax = np.rint(ymax)
    
    if verbose:
        print("Extent in x: {} to {}".format(xmin,xmax))
        print("Extent in y: {} to {}".format(ymin,ymax))
        
    return xmin, xmax, ymin, ymax



def make_overlay_after_rotatevol(ib_tif, cf_orig_tif, cf_data_rot, 
                                 rotation_matrix, scale, translation,
                                 verbose=True):
    """
    Create overlay after rotating fluo image with IMOD's rotatevol function.
    Note: make_overlay_data_multifluo is more general (allows multiple fluo channels) and can be used instead.

    Parameters
    ----------
    ib_tif : str
        Filename of target image, e.g. ion beam / SEM / TEM image. Should be a tif file.
    cf_orig_tif : str
        Filename of original fluo image (cf stands for confocal, widefield can also be used). Should be a tif file.
    cf_data_rot : np.ndarray
        Array of rotated fluo data.
    rotation_matrix : np.ndarray
        Rotation matrix, needed to calculate position of fluo image relative to target image. Calculated from correlation output euler angles.
    scale : float
        Scaling factor of fluo relative to target image.
    translation : np.ndarray
        (x,y,z) translation array.
    verbose : bool, optional
        The default is True.

    Returns
    -------
    composite_data : np.ndarray
        Array containing target and rotated fluo image.

    """
    if verbose:
        print("Loading image data..")
    # load IB image
    ib_data = tf.imread(ib_tif)
    if len(ib_data.shape) == 3:
        ib_data = ib_data[:,:,0]
    # get size of orig cf data
    with tf.TiffFile(cf_orig_tif) as tif:
        cf_orig_size = tif.asarray().shape
        cf_orig_dtype = tif.pages[0].dtype
    
    # Make MIP of rotatevol data
    cf_rot_MIP = cf_data_rot.max(axis=0)
    
    if verbose:
        print("Preparing interpolation..")
    # Get extent of rotated fluo data after rotation, scaling and translation
    (xmin, xmax, ymin, ymax) = extent_after_rotation(np.flip(cf_orig_size), np.flip(cf_data_rot.shape),
                                                 rotation_matrix, scale, translation,
                                                 round_values=True, verbose=verbose)
    # Range of rotatevol data
    x_range_rotatevol = np.linspace(xmin, xmax, cf_rot_MIP.shape[1])
    y_range_rotatevol = np.linspace(ymin, ymax, cf_rot_MIP.shape[0])
    
    # Target grid
    grid_x, grid_y = np.mgrid[xmin:xmax:(xmax-xmin+1)*1j, 
                              ymin:ymax:(ymax-ymin+1)*1j]
    # Make array of grid coordinates
    grid_coordinates = np.hstack( ( np.reshape(grid_x,(-1,1)) , np.reshape(grid_y, (-1,1)) ) )
    
    # Interpolation with scipy.interpolate.interpn
    MIP_interpol = interpn((x_range_rotatevol, y_range_rotatevol),
                             cf_rot_MIP.T,
                             grid_coordinates,
                             method='linear',
                             bounds_error=False,
                             fill_value=None)
    
    if verbose:
        print("Processing interpolation results..")
    # reshape result corresponding to grid_coordinates back to target dimensions
    x_size = int(xmax-xmin+1)
    MIP_interpol_reshaped = np.reshape(MIP_interpol, (x_size, -1))
    
    # MIP interpol values need to be transposed (x,y->y,x), flipped along y and value range adjusted to 0,255
    MIP_interpol_transformed = np.flip(MIP_interpol_reshaped.T, 0) - MIP_interpol_reshaped.min()
    # Embed reshaped data into array of same size as IB image
    MIP_large = np.zeros(ib_data.shape)
    # Define ranges
    L_xmin, L_ymin, L_xmax, L_ymax = int(xmin), int(ymin), int(xmax+1), int(ymax+1)
    I_xmin = I_ymin = 0
    (I_ymax, I_xmax) = MIP_interpol_transformed.shape
    # If MIP_interpol has values outside range of MIP_large:
    if xmin<0:
        L_xmin = 0
        I_xmin = int(-1*xmin)
    if ymin<0:
        L_ymin = 0
        I_ymin = int(-1*ymin)
    if L_xmax > MIP_large.shape[1]:
        L_xmax = MIP_large.shape[1]
        I_xmax = I_xmin + MIP_large.shape[1]
    if L_ymax > MIP_large.shape[0]:
        L_ymax = MIP_large.shape[0]
        I_ymax = I_ymin + MIP_large.shape[0]

    MIP_large[L_ymin:L_ymax, L_xmin:L_xmax] = MIP_interpol_transformed[I_ymin:I_ymax, I_xmin:I_xmax]
    # Convert fluo to same dtype as orig
    if MIP_large.dtype != cf_orig_dtype:
        MIP_large = MIP_large.astype(cf_orig_dtype)
    # Create composite array
    composite_data = np.array([ib_data, MIP_large])
    
    if verbose:
        print("Done!")
    return composite_data
    
# Same as make_overlay_after_rotatevol, but for multiple fluo channels    
def make_overlay_data_multifluo(ib_tif, cf_orig_tif, list_cf_data_rot, 
                                rotation_matrix, scale, translation,
                                cf_data_as_stack=True, verbose=True):
    """
    Create overlay after rotating fluo images with IMOD's rotatevol function.

    Parameters
    ----------
    ib_tif : str
        Filename of target image, e.g. ion beam / SEM / TEM image. Should be a tif file.
    cf_orig_tif : str
        Filename of original fluo image (cf stands for confocal, widefield can also be used). Should be a tif file.
    list_cf_data_rot : list of np.ndarrays
        List of arrays of rotated fluo data.
    rotation_matrix : np.ndarray
        Rotation matrix, needed to calculate position of fluo image relative to target image. Calculated from correlation output euler angles.
    scale : float
        Scaling factor of fluo relative to target image.
    translation : np.ndarray
        (x,y,z) translation array.
    cf_data_as_stack : bool, optional
        If True, assumes cf_data are given as stacks. If False, assumes they're already MIPs.
    verbose : bool, optional
        The default is True.

    Returns
    -------
    composite_data : np.ndarray
        Array containing target and rotated fluo images.

    """
    if verbose:
        print("Loading image data..")
    # load IB image
    ib_data = tf.imread(ib_tif)
    if len(ib_data.shape) == 3:
        ib_data = ib_data[:,:,0]
    # get size of orig cf data
    with tf.TiffFile(cf_orig_tif) as tif:
        cf_orig_size = tif.asarray().shape
        cf_orig_dtype = tif.pages[0].dtype
        
    # Get list of MIPs
    if cf_data_as_stack:
        list_cf_MIP = [data.max(axis=0) for data in list_cf_data_rot]
    else:
        list_cf_MIP = list_cf_data_rot
    # Use first fluo image to determine parameters
    cf_rot_MIP = list_cf_MIP[0]
    
    if verbose:
        print("Preparing interpolation..")
    # Get extent of rotated fluo data after rotation, scaling and translation
    (xmin, xmax, ymin, ymax) = extent_after_rotation(np.flip(cf_orig_size), np.flip(cf_rot_MIP.shape),
                                                 rotation_matrix, scale, translation,
                                                 round_values=True, verbose=verbose)
    # Range of rotatevol data
    x_range_rotatevol = np.linspace(xmin, xmax, cf_rot_MIP.shape[1])
    y_range_rotatevol = np.linspace(ymin, ymax, cf_rot_MIP.shape[0])
    
    # Target grid
    grid_x, grid_y = np.mgrid[xmin:xmax:(xmax-xmin+1)*1j, 
                              ymin:ymax:(ymax-ymin+1)*1j]
    # Make array of grid coordinates
    grid_coordinates = np.hstack( ( np.reshape(grid_x,(-1,1)) , np.reshape(grid_y, (-1,1)) ) )
    
    # Start composite data np array 
    composite_data = np.array([ib_data])
    for i, channel_rot_MIP in enumerate(list_cf_MIP):
        if verbose:
            print("Starting interpolation for channel {}".format(i))

        # Interpolation with scipy.interpolate.interpn
        MIP_interpol = interpn((x_range_rotatevol, y_range_rotatevol),
                                 channel_rot_MIP.T,
                                 grid_coordinates,
                                 method='linear',
                                 bounds_error=False,
                                 fill_value=None)
        
        if verbose:
            print("Processing interpolation results..")
        # reshape result corresponding to grid_coordinates back to target dimensions
        x_size = int(xmax-xmin+1)
        MIP_interpol_reshaped = np.reshape(MIP_interpol, (x_size, -1))
        
        # MIP interpol values need to be transposed (x,y->y,x), flipped along y and value range adjusted to 0,255
        MIP_interpol_transformed = np.flip(MIP_interpol_reshaped.T, 0) - MIP_interpol_reshaped.min()
        # Embed reshaped data into array of same size as IB image
        MIP_large = np.zeros(ib_data.shape)
        # Define ranges
        L_xmin, L_ymin, L_xmax, L_ymax = int(xmin), int(ymin), int(xmax+1), int(ymax+1)
        I_xmin = I_ymin = 0
        (I_ymax, I_xmax) = MIP_interpol_transformed.shape
        # If MIP_interpol has values outside range of MIP_large:
        if xmin<0:
            L_xmin = 0
            I_xmin = int(-1*xmin)
        if ymin<0:
            L_ymin = 0
            I_ymin = int(-1*ymin)
        if L_xmax > MIP_large.shape[1]:
            L_xmax = MIP_large.shape[1]
            I_xmax = I_xmin + MIP_large.shape[1]
        if L_ymax > MIP_large.shape[0]:
            L_ymax = MIP_large.shape[0]
            I_ymax = I_ymin + MIP_large.shape[0]
    
        MIP_large[L_ymin:L_ymax, L_xmin:L_xmax] = MIP_interpol_transformed[I_ymin:I_ymax, I_xmin:I_xmax]

        # Convert fluo to same dtype as original
        if MIP_large.dtype != cf_orig_dtype:
            MIP_large = MIP_large.astype(cf_orig_dtype)
        # Create composite array
        composite_data = np.concatenate([ composite_data, np.array([MIP_large]) ], axis=0)
    
    if verbose:
        print("Done!")
    return composite_data    
    
def run_IMOD_rotatevol(in_file, out_file, angles, verbose=True):
    """
    Rotate input file by the given angles, using the function rotatevol from IMOD.
    Dependencies: IMOD. 
    rotatevol is run using subprocess. Output is directly saved as out_file.

    Parameters
    ----------
    in_file : str
        Filename of input file, usually tif.
    out_file : str
        Filename of output file, usually mrc.
    angles : list
        List of euler angles for rotatevol.
    verbose : bool, optional
        The default is True.

    Returns
    -------
    None.

    """
    # Shell commands are run using subprocess
    # useful package to find out how to split commands into a list: import shlex > shlex.split(command_str)
    
    rot_query_list = ['rotatevol', '-i', in_file, 
                      '-angles', '{} {} {}'.format(angles[0], angles[1], angles[2]),
                      '-query']
    # run rotatevol query in system shell
    size = subprocess.check_output(rot_query_list, text=True) # text=True to get string and not bytes output
    size = size.strip() # orig output has \n after numbers
    # parameters for actual rotatevol
    rot_query_list = ['rotatevol', 
                  '-i', in_file, 
                  '-angles', '{} {} {}'.format(angles[0], angles[1], angles[2]),
                  '-size', size,
                  '-fill', '0',
                  '-output', out_file]

    rotatevol_out = subprocess.run(rot_query_list, capture_output=True, text=True) 
    if verbose:
        print(rotatevol_out.stdout)

def plot_overlay_IMOD(fname_corr_txt, fname_target, fluo_file_list, list_colors=None, ax=None,
                      return_fnames_rot=False):
    """
    Given input image files and the 3DCT correlation result, use IMOD to rotate the fluorescence data and generate an overlay.

    Parameters
    ----------
    fname_corr_txt : str
        Filename of correlation text file.
    fname_target : str
        Filename of target image, e.g. ion beam / SEM / TEM image. Should be a tif file.
    fluo_file_list : list of str
        List of resliced fluorescence images to overlay on target image. Should be tif files.
    list_colors : list of tuples or lists, optional
        List of colors corresponding to the input fluorescence images. 
        If None, the default color cycle is red - green - blue. The default is None.
    ax : matplotlib.pyplot.axis, optional
        Figure axis to plot output into. If None, a new figure is generated. The default is None.
    return_fnames_rot : bool, optional
        If True, returns a list of filenames of rotated fluo mrc files.

    Returns
    -------
    composite_data : np.ndarray
        Array containing target image and fluorescence overlays.
    im_list : list
        List of plotted images.
    list_fnames_fluo_rot : list of str, optional, if return_fnames_rot==True
        List of file names of rotated fluo mrc files, can be used for cleanup.

    """
    print("Started plot_overlay_IMOD")
    
    # Take care of ax
    if ax is None:
        fig, ax = plt.subplots(dpi=300)
    
    # Take care of colors
    if list_colors is None:
        list_colors = [(1,0,0),(0,1,0),(0,0,1)]
    elif not isinstance(list_colors[0],(list,tuple)):
        list_colors = [list_colors] # makes a list
    list_colors = cycle(list_colors)
    
    # Read parameters from 3DCT results
    p = read_3dct_results(fname_corr_txt, verbose=True)
    
    # Rotation matrix
    # gives same results as first making matrix with make_r_euler or pyto and then converting into zyx angles
    # rotation matrix can be generated by r.as_matrix()
    
    r = Rotation.from_euler('zxz', [p['phi'], p['theta'], p['psi'] ], degrees=True)
    zyx_angles = r.as_euler('zyx', degrees=True)
    rotatevol_angles = [-zyx_angles[0], zyx_angles[1], -zyx_angles[2]]
    
    list_fluo_rot = []
    list_fnames_fluo_rot = []
    # Rotate and load fluo data
    for i, fname_fluo in enumerate(fluo_file_list):
        print("Rotating channel {} with IMOD rotatevol".format(i))
        # Build rotatevol output name
        fstem_target = Path(fname_target).stem # NEW: add target image name to rot output name, in case several different overlays are made (target name seems more readable than correlation txt name)
        fname_rot_mrc = fname_fluo.replace('.tif', f'_rot_{fstem_target}.mrc')
        list_fnames_fluo_rot.append(fname_rot_mrc)
        
        # Only run rotatevol if rotated stack doesn't exist already
        if not os.path.exists(fname_rot_mrc):
            # Run rotatevol
            run_IMOD_rotatevol(fname_fluo, fname_rot_mrc, rotatevol_angles, verbose=True)
        # Load mrc output of rotatevol
        print("Loading channel {} after rotation".format(i))
        with mrcfile.mmap(fname_rot_mrc) as mrc: # mmap might be faster than open
            mrc_data = mrc.data.max(axis=0) # Load MIP directly
        
        # Rotatevol makes 8bit files into int8, but leaves intensities of 16 bit files
        if mrc_data.dtype == 'int8':
            mrc_data = (mrc_data + 128).astype(np.uint8)
        else:
            # Otherwise, simply make rotated data same dtype as input
            with tf.TiffFile(fname_fluo) as tif:
                dtype_fluo = tif.pages[0].dtype
            mrc_data = mrc_data.astype(dtype_fluo)
                
        list_fluo_rot.append( mrc_data ) # directly convert to uint8
    # Generate overlay data
    print("Making overlay data")
    composite_data = make_overlay_data_multifluo(fname_target, 
                                                 fluo_file_list[0],
                                                 list_fluo_rot,
                                                 r.as_matrix(), p['scale'], p['trans'],
                                                 cf_data_as_stack=False,
                                                 verbose=True)
    # Make plot
    print("Plotting overlay data")
    im = ax.imshow(composite_data[0], cmap='Greys_r', origin='upper')
    ax.set_xlim(ax.get_xlim()) # fix x,y lim
    ax.set_ylim(ax.get_ylim())   

    im_list = [im]
    for i, (data, color) in enumerate(zip(composite_data[1:], list_colors)):
        print("Plotting fluo channel {}".format(i))
        cmap = cmap_fading(color)
        im = ax.imshow(data, cmap=cmap)
        im.set_label('fluo_{}'.format(i))        
        im_list.append(im)
    
    if return_fnames_rot:
        return composite_data, im_list, list_fnames_fluo_rot
    
    return composite_data, im_list

        
def save_composite_tif_fromarray(composite_data, fname, list_colors=None):
    """
    Save overlay image from composite data.

    Parameters
    ----------
    composite_data : np.ndarray
        Output array of plot_overlay_IMOD. Array of target image and rotated overlay fluo images
    fname : str
        Output file name.
    list_colors : list of tuples/lists, optional
        List of colors for the output image, should contain black (-> gray cmap) for the target image. 
        If None, the color cycle is: black - red - green - blue. The default is None.

    Returns
    -------
    None.

    """
    print("Directly saving composite data!")
    # Take care of colors
    if list_colors is None:
        list_colors = [(1,1,1),(1,0,0),(0,1,0),(0,0,1)]
    elif not isinstance(list_colors[0],(list,tuple)):
        list_colors = [list_colors] # makes a list
    list_colors = cycle(list_colors)
    lut_list = []
    for color, data in zip(list_colors, composite_data):
        lut_list.append(lut_imagej(color))
        
    print('Saving overlay image as {}'.format(fname))

    tf.imsave(fname,
              composite_data,
              imagej=True,
              metadata={'mode': 'composite'},
              ijmetadata={'LUTs': lut_list}) 
    
    
    
    
    
    
    
    
    

