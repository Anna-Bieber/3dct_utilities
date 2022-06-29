# -*- coding: utf-8 -*-
"""
Created on Fri Feb 14 18:19:29 2020

Batch script for reslicing stacks for 3DCT and for creating MIPs
Reslicing script adapted from original 3DCT reslicing function

Dependencies: 
    tifffile for loading stacks (pip install tifffile)
    stack_tools module (folder) including __init__.py, MIP_functions.py and tdct_reslice.py

@author: Anna
"""

import os
import re
import matplotlib.pyplot as plt
import logging
from pathlib import Path
import tifffile as tf 

from stack_tools import tdct_reslice, MIP_functions

#%% INPUT

logging.basicConfig(level=logging.DEBUG) # set logging levels, e.g. logging.DEBUG or logging.INFO

reslice = True
make_MIP = True
MIP_mode = 'sum' # 'sum' for Z projection, 'max' for Maximum Intensity Projection

step_in = 300.595 # Input focus step (z step size)
step_out = 89.018 # Output focus step (xy step size, same unit as step_in)

in_dir = Path(r'/fs/pool/pool-erdmann/Cristina/Yeast_autophagy/Ede1/EA_GFP_Ede1_RFP_Atg8/EA017_6/20220202_EA017/batch_decon_green')
out_subdir = 'resliced'
MIP_subdir = 'MIP_decon'

out_dir = in_dir / out_subdir
MIP_dir = in_dir / MIP_subdir

channels_to_process = ['cmle']#['ch00', 'ch01'] # If only one channel, no channel number is given in tif file

color_list = [(0,1,0)] # color list with normalized (r,g,b) values for channels in right order
hist_clip_factor = 0.9 # histogram clipping factor for enhancing signal

#%%

log = logging.getLogger("Main")

# Make output directory
if not os.path.exists(out_dir):
    log.debug("Making output directory {}".format(out_dir))
    os.mkdir(out_dir)

if not os.path.exists(MIP_dir):
    log.debug("Making MIP output directory {}".format(MIP_dir))
    os.mkdir(MIP_dir)

# find and sort files by channels and FOVs
d = {}

log.debug("Sorting files to process...")
for channel in channels_to_process:
    files = list( in_dir.glob('*{}*.tif'.format(channel)) )
    files.sort()
    d[channel] = {}
    
    for file in files:
        file = file.as_posix()
        if file.endswith('_resliced.tif'):
            continue
        FOV = re.search('(FOV\d*)_', file)
        FOV = FOV.group(1)
        
        d[channel][FOV] = file

# Make list of FOVs
k = channels_to_process[0]
FOV_list = list(d[k].keys())       

log.debug("...done.")        
#%%

for FOV in FOV_list:    # loop over FOVs
    log.debug("Processing {}".format(FOV))
    stack_list = []
    for channel in channels_to_process:
        log.debug("Reading {} of {}".format(channel, FOV))
        fname_in = d[channel][FOV]
        stack = tf.imread(fname_in)
        # sometimes output of tf.imread has 4 dimensions
        if len(stack.shape) > 3 and stack.shape[0] == 1:
            stack = stack[0]
        log.debug("...done.")
        # reslice stack
        if reslice:
            fname_out = out_dir / os.path.basename(fname_in).replace('.tif', '_resliced.tif')
            fname_out = fname_out.as_posix()
            log.debug('Resliced stack file: {}'.format(fname_out))
            # Reslice:
            tdct_reslice.reslice(stack, step_in, step_out, interpolationmethod='linear', 
               fname_out=fname_out)
        # add to stack list
        stack_list.append(stack) # stack list contains all stacks before reslicing (for MIP)
    
    if make_MIP:
        fig, ax = plt.subplots() 
        
        # Make plot
        im_list = MIP_functions.plot_MIP(stack_list, color_list, axis=0, ax=ax, mode=MIP_mode, clim_factor=hist_clip_factor)
        fig.set_facecolor("black")
           
        MIP_file = MIP_dir / '{}_MIP_decon.png'.format(FOV)
        MIP_file = MIP_file.as_posix()
       
        # Save MIP as png
        fig.savefig(MIP_file, dpi=300, bbox_inches='tight')
        
        plt.close()

