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
typ_dict = {'bm': 'brainmask', 'wm': 'wm', 'bfs': 'brain.finalsurfs', 'cp': 'control.dat'}
vol = typ_dict[typ]


def get_diff_data(editp=editp, uneditp=uneditp, vol=vol, subnum=subnum):
    # Calculate difference image
    edit_img = nib.load(path.join(editp,'mri/%s.mgz')%(vol))
    edit_data = edit_img.get_fdata()
    unedit_img = nib.load(path.join(uneditp,'sub-%s/mri/%s.mgz')%(subnum, vol))
    unedit_data = unedit_img.get_fdata()
    diff_data = unedit_data - edit_data
    return(diff_data)

def summarize_edits(diff_data):

    # Extract non-zero values from difference data and arrange in df
    out = pd.DataFrame(np.asarray(np.asarray(diff_data != 0).nonzero()).T).rename(columns={0:"Sag", 1:"Axe", 2:"Cor"})
    out['diff_val'] = diff_data[np.where(diff_data != 0)]
    out['Action'] = np.where(out.diff_val>0, "delete voxel", "add voxel")
    out = out.drop(columns="diff_val")
    out['Vol'] = vol
    out = out.sort_values(by=['Cor'])
    out.reset_index(drop=True)

    # Make distance matrix
    dm = distance_matrix(out[['Sag', 'Axe', 'Cor']], out[['Sag', 'Axe', 'Cor']])

    # Apply clustering on the distance matrix
    # eps: The maximum distance between two samples for one to be considered as in the neighborhood of the other. This is not a maximum bound on the distances of points within a cluster. This is the most important DBSCAN parameter to choose appropriately for your data set and distance function.
    # min_samples: The number of samples (or total weight) in a neighborhood for a point to be considered as a core point. This includes the point itself.
    clustering = DBSCAN(eps=50, min_samples=5).fit(dm)

    # Assign cluster labels to use as edit groups later
    out['cluster'] = clustering.labels_

    # Roll up by edit group (i.e. cluster) summarizing the beginning and the end of each cluster in each orientation
    out_grouped = out.groupby('cluster')
    out_grouped = out_grouped.agg(min_sag=('Sag', min),
                   max_sag=('Sag', max),
                   min_axe=('Axe', min),
                   max_axe=('Axe', max),
                   min_cor=('Cor', min),
                   max_cor=('Cor', max),
                   num_vox=('Cor', 'count')).reset_index()
    out_grouped['vol'] = out.Vol.unique()[0]
    out_grouped['action'] = out.Action.unique()[0]
    return(out_grouped)

def get_bm_edits(editp=editp, uneditp=uneditp, vol=vol, subnum=subnum):
    
    bm_diff_data = get_diff_data(editp=editp, uneditp=uneditp, vol=vol, subnum=subnum)
    bm_edits = bm_diff_data.groupby('action').apply(summarize_edits).reset_index()

    return bm_edits

def get_wm_edits(editp=editp, uneditp=uneditp, vol=vol, subnum=subnum):
    
    wm_diff_data = get_diff_data(editp=editp, uneditp=uneditp, vol=vol, subnum=subnum)
    wm_edits = wm_diff_data.groupby('action').apply(summarize_edits).reset_index()

    return wm_edits

def get_cp_edits(subnum = subnum):
    
    # Check if control points exist
    if path.isfile(path.join(editp,'tmp/%s')%(vol)):
        # Read in control.dat
    
    
        # Translating to voxel space/slice numbers comparable to bm and wm edits
        # http://www.grahamwideman.com/gw/brain/fs/coords/fscoords.htm
        # Right = 128 - Byte
        # Ant = Slice - 128
        # Sup = 128 - Row

        # control.dat format TkReg RAS
        # -36 -12 24
        # Translation to [Sag, Axe, Cor]
        # [164, 104, 116]
        # -36 = 128 - Sag; Sag = 128 - (-36) = 164
        # -12 = Axe - 128; Axe = 128 + (-12) = 116
        # 24 = 128 - Cor; Cor = 128 - (24) = 104
    
    else:
        print("No control points found for %s"%(subnum))
        out = pd.DataFrame()
    
    # Always returns something to avoid exceptions when looping through and appending all types of edits below
    # (might not be necessary)
    return out

def get_bfm_edits(subnum = subnum):
    return out

fn_dict = {'bm': get_bm_edits, 'wm': get_wm_edits, 'cp': get_cp_edits, 'bfm': get_bfm_edits}

if vol == 'all':
    out = pd.DataFrame()
    # loop through function lookup dictionary and run all edit extraction functions
    for k,v in fn_dict.items():
        out = out.append(fn_dict[k]())
else:
    # apply corresponding function as specified in the dicitonary of functions
    out = fn_dict[typ]() 

out.to_csv(path.join(outdir, subnum+'_'+typ+'_edits.csv'))

# Clean up: delete unzipped unedited directory
os.rmdir(uneditp)
