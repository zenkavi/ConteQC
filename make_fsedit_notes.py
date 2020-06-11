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
typ_dict = {'bm': 'brainmask', 'wm': 'wm', 'cp': 'control.dat', 'bfs': 'brain.finalsurfs'}
fn_dict = {'bm': get_bm_edits, 'wm': get_wm_edits, 'cp': get_cp_edits, 'bfs': get_bfs_edits}

if typ == 'all':
    out = pd.DataFrame()
    # loop through function lookup dictionary and run all edit extraction functions
    for k,v in fn_dict.items():
        vol = typ_dict[k]
        out = out.append(v(editp=editp, uneditp=uneditp, vol=vol, subnum=subnum))
else:
    # apply corresponding function as specified in the dicitonary of functions
    vol = typ_dict[typ]
    out = fn_dict[vol](editp=editp, uneditp=uneditp, vol=vol, subnum=subnum) 

out.to_csv(path.join(outdir, subnum+'_'+typ+'_edits.csv'))

# Clean up: delete unzipped unedited directory
os.rmdir(uneditp)
