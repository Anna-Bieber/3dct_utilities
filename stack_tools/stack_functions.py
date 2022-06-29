#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 16 09:35:33 2021

@author: anbieber
"""
import numpy as np


def bin_stack_xy(stack, bin_factor=2):
    """Bin a zyx stack in yx. Code adapted from (11.02.2021): 
        scipython.com/blog/binning-a-2d-array-in-numpy/"""
    bin_factor = int(bin_factor)
    # If x or y is not divisible by bin factor, get rid of values in the back
    for i in [1,2]:
        i = int(i)
        division_rest = stack.shape[i] % bin_factor
        if division_rest > 0:
            id_range = np.arange(stack.shape[i] - division_rest)
            stack = np.take(stack, id_range, axis=i)
            
            
    # Determine original and tmp shape
    orig_shape = stack.shape    
    tmp_shape = (orig_shape[0], int(orig_shape[1]/bin_factor), bin_factor,  int(orig_shape[2]/bin_factor), bin_factor)
    
    # Reshape and average along relevant axes (the ones with length of bin_factor)
    stack_reshaped = stack.reshape(tmp_shape)    
    stack_binned = stack_reshaped.mean(-1).mean(2)
    
    return stack_binned

def bin_stack_xyz(stack, bin_factor=2):
    """Bin a zyx stack in yx. Code adapted from (11.02.2021): 
        scipython.com/blog/binning-a-2d-array-in-numpy/"""
    bin_factor = int(bin_factor)
    # If z, x or y is not divisible by bin factor, get rid of values in the back
    for i in [0,1,2]:
        i = int(i)
        division_rest = stack.shape[i] % bin_factor
        if division_rest > 0:
            id_range = np.arange(stack.shape[i] - division_rest)
            stack = np.take(stack, id_range, axis=i)
            
            
    # Determine original and tmp shape
    orig_shape = stack.shape    
    tmp_shape = (int(orig_shape[0]/bin_factor), bin_factor, 
                 int(orig_shape[1]/bin_factor), bin_factor,  
                 int(orig_shape[2]/bin_factor), bin_factor)
    
    # Reshape and average along relevant axes (the ones with length of bin_factor)
    stack_reshaped = stack.reshape(tmp_shape)    
    stack_binned = stack_reshaped.mean(5).mean(3).mean(1)
    
    return stack_binned, stack_reshaped


# Binning tests
# test_data = np.arange(1,49,1).reshape((3,4,4))

# shape_tmp = (3,2,2,2,2)
# test_tmp = test_data.reshape(shape_tmp)
# test_tmp2 = test_tmp.mean(-1).mean(2)