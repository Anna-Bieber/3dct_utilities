# 3dct_overlay

Create overlays of fluorescence data on target images after correlating with 3DCT. Scripts are written in Python 3 but use IMOD for fast rotation of the fluorescence stacks.

## Scripts

* overlay_3dct_stack.py: Main script for creating overlays.
* overlay_IMOD_functions_v1.py: Functions needed for creating overlays.

## Dependencies:

* Python 3 and standard packages including numpy, matplotlib, scipy, pathlib, itertools, subprocess
* Extra python packages: tifffile, mrcfile
* IMOD
