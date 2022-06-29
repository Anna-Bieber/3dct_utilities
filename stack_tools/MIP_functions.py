# -*- coding: utf-8 -*-
"""
Created on Sun Feb 16 17:57:51 2020

@author: Anna
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, rgb2hex, to_rgb
from itertools import cycle

#NOT DONE YET
def plot_MIP(stack_list, color_list, axis=0, ax=None, mode='max', clim_factor=1):
    "Makes maximum intensity projection along `axis`, overlaying channels"
    if color_list is None:
        color_list = ['red', 'green', 'blue']
    # Check stack and color list
    if not isinstance(stack_list, (list, tuple)):
        stack_list = [stack_list,]
        color_list = [color_list,]
    # Check color list
    color_list = [to_rgb(c) for c in color_list] # turn any type of color input into rgb
    color_list = cycle(color_list)
    
    if ax is None:
        fig, ax = plt.subplots()
    im_list = []
    im = ax.imshow(np.ones([s for i,s in enumerate(stack_list[0].shape) if i != axis]),
                   cmap=plt.cm.Greys_r)
    for stack, color in zip(stack_list, color_list):
        assert stack.ndim == 3, "Wrong stack dimension"
        if mode == 'max':
            img = np.max(stack, axis=axis)
        elif mode == 'sum':
            img = np.sum(stack, axis=axis)
        cmax = img.max()
        climit = cmax*clim_factor
        #img = np.ma.array(img, mask=img<1)
        cmap = cmap_fading(color)
        im = ax.imshow(img, cmap=cmap, clim = [0,climit])
        im_list.append(im)
        
    ax.axis('off')
        
    return im_list
        
    
def cmap_fading(base_color):
    # scale base color to 1 
    if np.max(base_color) > 1:
        c = [val/255 for val in base_color] # adjust range
    else:
        c = list(base_color)
    colors = np.zeros((1000, 4), np.float64)
    for i, val in enumerate(c):
        colors[:,i] = val # set RGB to base color
    colors[:,3] = np.linspace(0, 1, colors.shape[0]) # make alpha gradient
    name = 'fading_{}'.format(rgb2hex(c))
    cmap = ListedColormap(colors, name=name)
    plt.cm.register_cmap(name=name, cmap=cmap)
    
    return cmap