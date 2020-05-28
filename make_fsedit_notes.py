#!/Users/zeynepenkavi/anaconda/envs/py37/bin/python

"""
    Creates report of manual edits made on Freesurfer recon outputs

    Parameters:
        subnum = subject number (e.g. CC0058_core1)
        typ = type of edits to extract reports for
        basedir = directory where to find both the directory of edited images and compressed unedited images
        outdir = director where edit report will be saved into

    Returns:
        sub-[subnum].csv
"""

from argparse import ArgumentParser
import os
from os import path
import pandas as pd
import sys

parser = ArgumentParser()
parser.add_argument("--subnum")
parser.add_argument("--typ")
parser.add_argument("--basedir", default='/Users/zeynepenkavi/Downloads/ConteQC')
parser.add_argument("--outdir", default='/Users/zeynepenkavi/Downloads/ConteQC')
args = parser.parse_args()
subnum = args.subnum
typ = args.typ
basedir = args.basedir
outdir = args.outdir

# Check if original zip exists. If not print error and quit
exists = path.isfile(path.join(basedir, 'sub-', subnum, '.zip'))
if not exists:
    print("Unedited images not found in directory. Quitting...")
    try:
        sys.exit()
    except SystemExit:
        print("Quit because difference images cannot be calculated without unedited images.")
else:
    #Unzip unedited images


#Output:
#1. raw report; sort by coronal slice number in the end
#Case | Slice (Sag) | Slice (Axe) | Slice (Cor) | Action | Vol

#2. rolled up report (don't roll up just by distance; action must match as well)
#Case | Slice (Sag) | Slice (Axe) | Slice (Cor) | Action | Vol | Num_voxels

import os

os.system()

def get_bm_edits(subnum = subnum):
    return out

def get_wm_edits(subnum = subnum):
    return out

def get_cp_edits(subnum = subnum):
    return out

def get_bfm_edits(subnum = subnum):
    return out

fn_dict = {'bm': get_bm_edits, 'wm': get_wm_edits, 'cp': get_cp_edits, 'bfm': get_bfm_edits}

if vol == 'all':
    out = pd.DataFrame()
    for k,v in fn_dict.items():
        out = out.append(fn_dict[k]())
else:
    out = fn_dict[typ]() #apply corresponding function as specified in the dicitonary of functions
    out.to_csv(path.join(outdir, subnum+typ+'_edits.csv'))

#Steps:

#Unzip the original files again to be used for difference images
    #Make sure to save the new unzip with new suffix e.g. sub-SUBNUM_core1_unedited
#if vol != cp (control points)
#Convert the .mgz into .nii.gz using freesurfer > mri_convert
#Calculate difference image between the original and edited nii.gz

#
