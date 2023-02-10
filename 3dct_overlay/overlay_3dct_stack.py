#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  9 15:48:22 2021

Script for creating overlays of fluorescent stacks on FIB/SEM/TEM images after 3DCT
uses IMOD to perform stack rotation

@author: anbieber
"""

from pathlib import Path
from overlay_IMOD_functions_v1 import plot_overlay_IMOD, save_composite_tif_fromarray

#%% INPUT

# Paths to directories with files. Directory with fluorescence stacks can be set as different path
in_dir = Path('./test_data/')
fluo_dir = in_dir


# File names
fname_target = 'FOV8_IB_before.tif' # Target image (usually ion beam / SEM / TEM image). Should be in in_dir.
fname_correlation_txt = '2020-01-16_15-29-23_correlation.txt' # Correlation txt output. Should be in in_dir.
# List of fluo file names. Should be resliced images as used for 3dct. Should be in fluo_dir.
list_fnames_fluo = ['20200113_EA006_5_FOV8_7_cmle_ch00_resliced.tif'] 

list_fluo_colors = [(1,0,1)] # List of fluo colors, given as (r,g,b) tuples.

# Output filename
fname_overlay = fname_target.replace('.tif', '_overlay.tif')

clean_rot_mrc = True # If True, delete the intermediate rotated fluo mrc files.


#%% Processing

# Generate full file paths from file names

fpath_target = (in_dir / fname_target).as_posix()
fpath_correlation = (in_dir / fname_correlation_txt).as_posix()

list_fpaths_fluo = [(fluo_dir / fname).as_posix() for fname in list_fnames_fluo]

fpath_overlay = (in_dir / fname_overlay).as_posix()

# Process data & plot

composite_data, im_list, list_fpaths_fluo_rot = plot_overlay_IMOD(fpath_correlation,
                                                                  fpath_target,                                              
                                                                  list_fpaths_fluo,
                                                                  list_fluo_colors,
                                                                  return_fnames_rot=True)

# Save overlay image
overlay_colors = [(1,1,1), *list_fluo_colors]
save_composite_tif_fromarray(composite_data, fpath_overlay, overlay_colors)


#%% Clean up

if clean_rot_mrc:
    for fname in list_fpaths_fluo_rot:
        Path(fname).unlink(missing_ok=True) # unlink removes the files

