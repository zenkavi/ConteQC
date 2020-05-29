#!/Users/zeynepenkavi/anaconda/envs/py2/bin/python

"""
    Creates report of manual edits made on Freesurfer recon outputs

    Parameters:
        subnum = subject number (e.g. CC0058_core1)
        typ = type of edits to extract reports for
        basedir = directory where to find both the directory of edited images and compressed unedited images
        outdir = director where edit report will be saved into

    Returns:
        sub-[subnum]-[typ]-edits.csv
"""

from argparse import ArgumentParser
import os
from os import path
import nibabel as nib
import numpy as np
import pandas as pd
import sys
from zipfile import ZipFile

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
exists = path.isfile(path.join(basedir, 'sub-'+subnum+'.zip'))
if not exists:
    print("Unedited images not found in directory. Quitting...")
    try:
        sys.exit()
    except SystemExit:
        print("Quit because difference images cannot be calculated without unedited images.")
else:
    # Unzip unedited images of potential interest
    toExtract = ['brainmask.mgz', 'wm.mgz', 'brain.finalsurfs.mgz']
    suffix = path.join('sub-'+subnum, 'mri')
    toExtract = [path.join(suffix, s) for s in toExtract]
    with ZipFile(path.join(basedir, 'sub-'+subnum+'.zip'), 'r') as zipObj:
        # Get a list of all archived file names from the zip
        listOfFileNames = zipObj.namelist()
        # Iterate over the file names
        for fileName in listOfFileNames:
            # Check filename endswith csv
            if fileName in toExtract:
                # Extract a single file from zip
                zipObj.extract(fileName, path.join(basedir, 'sub-'+subnum+'_unedited'))

# Store paths in vars for faster string comprehension
editp = path.join(basedir, 'sub-'+subnum)
uneditp = path.join(basedir, 'sub-'+subnum+'_unedited')
typ_dict = {'bm': 'brainmask', 'wm': 'wm', 'bfs': 'brain.finalsurfs'}
vol = typ_dict[typ]

# Calculate difference image
bm_edit_img = nib.load(path.join(editp,'mri/%s.mgz')%(vol))
bm_edit_data = bm_edit_img.get_fdata()
bm_unedit_img = nib.load(path.join(uneditp,'sub-%s/mri/%s.mgz')%(subnum, vol))
bm_unedit_data = bm_unedit_img.get_fdata()
diff_data = bm_unedit_data - bm_edit_data

# Extract non-zero values from difference data and arrange in df
out = pd.DataFrame(np.asarray(np.asarray(diff_data != 0).nonzero()).T).rename(columns={0:"Sag", 1:"Axe", 2:"Cor"})
out['diff_val'] = diff_data[np.where(diff_data != 0)]
out['Action'] = np.where(out.diff_val>0, "delete voxel", "add voxel")
out = out.drop(columns="diff_val")
out['Vol'] = vol
out = out.sort_values(by=['Cor'])
out.reset_index(drop=True)

def get_bm_edits(subnum = subnum):

#Output:
#1. raw report; sort by coronal slice number in the end
#Case | Slice (Sag) | Slice (Axe) | Slice (Cor) | Action | Vol

#2. rolled up report (don't roll up just by distance; action must match as well)
#Case | Slice (Sag) | Slice (Axe) | Slice (Cor) | Action | Vol | Num_voxels

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
    out = fn_dict[typ]() # apply corresponding function as specified in the dicitonary of functions

out.to_csv(path.join(outdir, subnum+'_'+typ+'_edits.csv'))

# Clean up: delete unzipped unedited directory
os.rmdir(uneditp)
