#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 10 10:38:08 2021

@author: anbieber
"""
import numpy as np
from pathlib import Path
from read_corr_txt import read_corr_txt

def write_points_3dct(p, fname):
    """Write out point coordinates in 3DCT style."""
    # Check that xyz are given, if not, add column of zeros
    if p.shape[1] < 3:
        p = np.concatenate([p, np.expand_dims(np.zeros(p.shape[0]), -1)], axis=1)
    np.savetxt(fname, p, fmt='%.2f', delimiter='\t')
    
def write_points_fiji(p, fname):
    """Write out point coordinates that can be loaded in FIJI."""
    np.savetxt(fname, p[:,0:2], fmt='%.2f', delimiter=' ', 
               header='point\n{:d}'.format(p.shape[0]), comments='')
    

def extract_3dct_coordinates(fname_txt, points_LR=['L','R'], output_style='3DCT', beads_only=True,
                             fname_out_labels=None):
    """
    Write out coordinates from a 3DCT results file.

    Parameters
    ----------
    fname_txt : str
        Path to the correlation results txt file.
    points_LR : list, optional
        List indicating which points to write out, 'L' for left, 'R' for right. The default is ['L','R'].
    output_style : str, optional
        '3DCT' or 'FIJI'. The default is '3DCT'.
    beads_only : True or False, optional
        Write out coordinates of beads only or spots as well. The default is True.
    fname_out_labels : None or list of str, optional
        List of labels for generating output filenames. If None, labels 'L' and 'R' will be used.. The default is None.

    Returns
    -------
    None.

    """
    # Read results file
    results = read_corr_txt(fname_txt)
    # Extract relevant coordinates
    points_all = {'L': {'beads': results['transformation_markers']['Initial (3D) markers'],
                        'spots': results['correlation_3d']['Spots (3D)']},
                  'R': {'beads': results['transformation_markers']['Final (2D) markers'],
                        'spots': results['correlation_3d']['Correlated spots']} }
    # Prepare output names
    fname_out_ending = '.csv' if output_style=='3DCT' else '.txt'
    if fname_out_labels == None:
        fname_out_list = ['corr_points_{}_{}{}'.format(output_style.lower(), pos, fname_out_ending) for pos in points_LR]
    else:
        fname_out_list = ['corr_points_{}_{}{}'.format(output_style.lower(), label, fname_out_ending) for label in fname_out_labels]

    # Make paths from output names
    out_dir = Path(fname_txt).parent
    fpaths_out = [(out_dir / f).as_posix() for f in fname_out_list]
            
    for key, filename in zip(points_LR, fpaths_out):
        # Choose which points to use: beads only or beads+spots
        if beads_only:
            p = points_all[key]['beads']
        else:
            p0 = points_all[key]['beads']
            if p0.shape[1] < 3:
                p0 = np.concatenate([p0, np.expand_dims(np.zeros(p0.shape[0]), -1)], axis=1)
            p = np.concatenate([p0, points_all[key]['spots']], axis=0)
        # Write out files
        if output_style.upper() == '3DCT':
            write_points_3dct(p, filename)
        elif output_style.upper() == 'FIJI':
            write_points_fiji(p, filename)