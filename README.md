# 3dct_utilities
Scripts related to 3dct functions.

## 3dct_overlay

Scripts to create an overlay of fluorescence data on top of ion beam / SEM / TEM data. IMOD is used for fast rotation of fluorescence stacks.

## stack tools

Functions for manipulating fluorescence stacks, e.g., making a maximum intensity projection, binning or reslicing. Reslicing functions are mostly taken from the 3DCT.

## Other scripts

- **Shift_corr_points_v2.ipynb**: 
	- Used for shifting a set of correlation points between similar FOVs, e.g. SEM overviews of the same grid square before and after milling. 
	- Usage: Do correlation of first image and save the correlation points. Then load the second image, shift 3-5 points to their correct new position and save them. Load both the original points and the partially shifted points as in_file and shift_file, respectively, and run the script to get a new file with all points shifted.
- **extract_3dct_coordinates**: Functions to extract point coordinates from the .txt output of 3DCT, either to reload points in 3DCT or to load them in FIJI. The _terminal.py script allows running this from a terminal.

- **read_corr_txt.py**: Function for reading the correlation txt file into python. This is a rather complicated version, there's an easier one in 3dct_overlay/overlay_IMOD_functions_v1.py

- **reslice_stack_single** and **reslice_stacks_batch**: For reslicing fluorescence data and optionally creating MIPs. Reslicing is done using the 3DCT functions. Especially the batch script is useful in preparation for a correlative FIB session, to avoid waiting times during reslicing and to get a quick overview of all available fields of view.
