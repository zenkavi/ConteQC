#!/Users/zeynepenkavi/anaconda/envs/py37/bin/python

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
import shutil
import sys
from zipfile import ZipFile
from fsedit_notes_helpers import get_diff_data, summarize_edits, get_bm_edits, get_wm_edits, get_cp_edits, get_bfs_edits

parser = ArgumentParser()
parser.add_argument("--subnum")
parser.add_argument("--typ")
parser.add_argument("--basedir", default='/Users/zeynepenkavi/Downloads/ConteQC')
parser.add_argument("--outdir", default='/Users/zeynepenkavi/Downloads/ConteQC')
parser.add_argument("--savediff", default=True)
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
    print("Extracting unedited files complete.")

# Store paths in vars for faster string comprehension
editp = path.join(basedir, 'sub-'+subnum)
uneditp = path.join(basedir, 'sub-'+subnum+'_unedited')
if savediff:
    diffp = path.join(basedir, 'sub-'+subnum+'_diff')
    if not os.path.exists(diffp):
        os.makedirs(diffp)
typ_dict = {'bm': 'brainmask', 'wm': 'wm', 'cp': 'control.dat', 'bfs': 'brain.finalsurfs'}
fn_dict = {'bm': get_bm_edits, 'wm': get_wm_edits, 'cp': get_cp_edits, 'bfs': get_bfs_edits}

if typ == 'all':
    out = pd.DataFrame()
    # loop through function lookup dictionary and run all edit extraction functions
    for k,v in fn_dict.items():
        print("Extracting %s edits"%(k))
        vol = typ_dict[k]
        out = out.append(v(editp=editp, uneditp=uneditp, vol=vol, subnum=subnum, savediff=savediff, diffp=diffp))
else:
    print("Extracting %s edits"%(typ))
    # apply corresponding function as specified in the dicitonary of functions
    vol = typ_dict[typ]
    out = fn_dict[vol](editp=editp, uneditp=uneditp, vol=vol, subnum=subnum, savediff=savediff, diffp=diffp) 

out.to_csv(path.join(outdir, subnum+'_'+typ+'_edits.csv'))

# Clean up: delete unzipped unedited directory
shutil.rmtree(uneditp)