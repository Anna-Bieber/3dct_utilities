# -*- coding: utf-8 -*-
"""
Created on Wed Jan 29 15:59:32 2020
Parser written by Micha

"""
import numpy as np

def read_corr_txt(fname, verbose=False):
    """
    Logfile parser for the '*_correlated.txt' files
    
    Parameters
    ----------
    fname : str
        path to the logifle
    verbose : bool, optional
        switch for additional information during read in.
        Default is False.
    
    Returns
    -------
    dict 
        Dictionary with informations from the logfile
    """
    def extract_list_from_line(line, start=None, end=None):
        line_split = line.split()
        alist = np.safe_eval(' '.join(line_split[start:end]))
        return alist
    
    results = {}
    
    if verbose: print('LogFileParser : read logfile : {}'.format(fname))
    with open(fname) as fp:
        text = fp.read().split('\n\n')  # HACK because we got two \n symbols...
        
        while len(text) > 0:
            line = text.pop(0)
            line_strip = line.strip()
            line_split = line_strip.split()
            if len(line_strip) == 0:  # empty line
                continue 
            elif line_strip == '# Transformation parameters':
                if verbose: print('LogFileParser : read Transformation parameters')
                trans_para = {}
                text.pop(0) # empty line
                # get the three angles from the line
                line = text.pop(0)
                trans_para['rotation'] = (' '.join(line.split()[-6:-3])[:-2],
                                          extract_list_from_line(line, start=-3))
                # get the float from the line
                trans_para['scale'] = float(text.pop(0).split()[-1])
                # get the two lists from the line
                line = text.pop(0)

                trans_para['translation_0'] = (
                    extract_list_from_line(line, 6, -4), # first 
                    extract_list_from_line(line, -3),    # second field
                )    

                line = text.pop(0)
                trans_para['translation_1'] = (
                    extract_list_from_line(line, 6, -4), # first 
                    extract_list_from_line(line, -3),    # second field
                )   

                trans_para['rms_error'] = float(text.pop(0).split()[-1])
                trans_para['optimization_status'] = text.pop(0).split()[-1]
                text.pop(0) # empty line
                results['transformation_parameters'] = trans_para
            elif line_strip == '# Transformation of initial (3D) markers':
                if verbose: print('LogFileParser : read Transformation of initial (3D) markers')
                text.pop(0) # empty line
                trans_markers = {
                    'Initial (3D) markers' : [],
                    'Transformed initial' : [],
                    'Final (2D) markers' : [],
                    'Transformed-Final' : [],
                }
                line = text.pop(0)
                header = [h for h in line.split('\t')[1:] if h != '']  # discard '#'
                if verbose: print(' -> header : {}'.format(header))
                # CHECK if everything is OK
                assert all([h in trans_markers.keys() for h in header])
                assert all([k in header for k in trans_markers.keys()])

                # parse lines
                while len(text) > 0:
                    line = text.pop(0)
                    line_strip = line.strip()
                    if len(line_strip) == 0 : continue
                    if line_strip[0] == '#' : break # if section ended
                    line_split = line_strip.split()
                    trans_markers['Initial (3D) markers'].append( tuple(map(float, line_split[0:3])) )
                    trans_markers['Transformed initial'].append(  tuple(map(float, line_split[3:6])) )
                    trans_markers['Final (2D) markers'].append(   tuple(map(float, line_split[6:8])) )
                    trans_markers['Transformed-Final'].append(    tuple(map(float, line_split[8:10])) )
                for k in list(trans_markers.keys()):
                    trans_markers[k] = np.array(trans_markers[k])
                results['transformation_markers'] = trans_markers
            elif line_strip == '# Correlation of 3D spots (POIs) to 2D':
                if verbose: print('LogFileParser : read Correlation of 3D spots (POIs) to 2D')
                text.pop(0) # empty line
                correlation_3d = {
                    'Spots (3D)' : [],
                    'Correlated spots' : [],
                }
                line = text.pop(0)
                header = [h for h in line.split('\t')[1:] if h != '']  # discard '#'
                if verbose: print(' -> header : {}'.format(header))
                # parse lines
                while len(text) > 0:
                    line = text.pop(0)
                    line_strip = line.strip()
                    if len(line_strip) == 0 : continue
                    if line_strip[0] == '#' : break # if section ended
                    line_split = line_strip.split()
                    correlation_3d['Spots (3D)'].append( tuple(map(float, line_split[0:3])) )
                    correlation_3d['Correlated spots'].append( tuple(map(float, line_split[3:6])) )
                for k in list(correlation_3d.keys()):
                    correlation_3d[k] = np.array(correlation_3d[k])    
                results['correlation_3d'] = correlation_3d
    return results