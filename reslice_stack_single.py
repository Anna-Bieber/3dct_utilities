#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 29 12:24:16 2020.

@author: anbieber
"""
import tifffile as tf
import os
from pathlib import Path
from stack_tools import tdct_reslice


step_in = 300.595 #.0000001 # Input focus step (z step size)
step_out = 89.018 # Output focus step (xy step size, same unit as step_in)

fname_in = '/fs/pool/pool-erdmann/Cristina/Yeast_autophagy/Ede1/EA_GFP_Ede1_RFP_Atg8/EA017_6/20220202_EA017/Stacks_GFP_FOV2_z00_ch00.tif'
out_dir = Path(fname_in).parent #Path('/fs/pool/pool-erdmann/Cristina/Meteor_development/20210608_EI002_1/FOV1/')

stack = tf.imread(fname_in)

 # sometimes output of tf.imread has 4 dimensions
if len(stack.shape) > 3 and stack.shape[0] == 1:
    stack = stack[0]
    
fname_out = out_dir / os.path.basename(fname_in).replace('.tif', '_resliced.tif')
fname_out = fname_out.as_posix()

tdct_reslice.reslice(stack, step_in, step_out, interpolationmethod='linear', fname_out=fname_out)


