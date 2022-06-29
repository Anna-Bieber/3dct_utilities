#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 10 10:40:28 2021

@author: anbieber
"""
import sys, os
import argparse



sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)))

from extract_3dct_coordinates import extract_3dct_coordinates

parser = argparse.ArgumentParser()
parser.add_argument('--input', help="Correlation results file (*_correlation.txt)")
parser.add_argument('--LR', default='LR', help="Which points to write out? L for left, R for right. (default: LR)")
parser.add_argument('--output_style', default='3DCT', help="3DCT or FIJI. default: 3DCT.")
parser.add_argument('--beads_only', default=True, help="If True, write only bead coordinates, otherwise also write out spot coordinates. default: True.")
parser.add_argument('--output_labels', default=None, nargs='+', help="Labels for output filenames (separate with space). If None, filenames are generated automatically. default is None.")

args = parser.parse_args()

fname_txt = args.input
points_LR = [s for s in args.LR]
output_style = args.output_style
beads_only = args.beads_only
if type(beads_only) == str:
    beads_only = True if beads_only.lower()=='true' else False
fname_out_labels = args.output_labels
if not all((fname_txt, points_LR, output_style)):
    print("ERROR: Not all arguments specified")
    parser.print_help()
    sys.exit(1)

print("""\
Extracting 3DCT coordinates with input:
- input : {fname_txt}
- LR   : {points_LR}
- output_style  : {output_style}
- beads_only : {beads_only}
- output_labels : {fname_out_labels}
""".format(**locals()))

# Run function
extract_3dct_coordinates(fname_txt, points_LR, output_style, beads_only, fname_out_labels)


