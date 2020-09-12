import os
from os import path
import nibabel as nib
import numpy as np
import pandas as pd
from scipy.spatial import distance_matrix
from sklearn.cluster import DBSCAN
from copy import copy

def get_diff_data(editp, uneditp, vol, subnum, savediff, diffp):
    # Calculate difference image
    
    if vol == "brain.finalsurfs":
        edit_img = nib.load(path.join(editp,'mri/%s.manedit.mgz')%(vol))
    else:
        edit_img = nib.load(path.join(editp,'mri/%s.mgz')%(vol))
    
    edit_data = edit_img.get_fdata()
    unedit_img = nib.load(path.join(uneditp,'sub-%s/mri/%s.mgz')%(subnum, vol))
    unedit_data = unedit_img.get_fdata()
    diff_data = unedit_data - edit_data
    
    if savediff:
        diff_img = nib.MGHImage(diff_data.astype(np.int32), edit_img.affine)
        nib.save(diff_img, os.path.join(diffp, 'diff_%s.mgz')%(vol))
        
    # Extract non-zero values from difference data and arrange in df
    out = pd.DataFrame(np.asarray(np.asarray(diff_data != 0).nonzero()).T).rename(columns={0:"Sag", 1:"Axe", 2:"Cor"})
    out['diff_val'] = diff_data[np.where(diff_data != 0)]
    out['Action'] = np.where(out.diff_val>0, "delete voxel", "add voxel")
    
    # When WM edits are not perfect it can look like voxel additions
    # For the clarity of the summary drop these
    if vol == "wm":
        out = out.query("diff_val != -1")
    
    out = out.drop(columns="diff_val")
    out['Vol'] = vol
    out = out.sort_values(by=['Cor'])
    out.reset_index(drop=True)
    
    return(out)

def summarize_edits(diff_data, vol):
    
    out = pd.DataFrame()
    
    for cur_action in diff_data['Action'].unique():
        cur_diff_data = diff_data.query("Action==@cur_action")
        
        # Make distance matrix
        #dm = distance_matrix(cur_diff_data[['Sag', 'Axe', 'Cor']], cur_diff_data[['Sag', 'Axe', 'Cor']])
        dm = distance_matrix(cur_diff_data.loc[:,('Sag', 'Axe', 'Cor')], cur_diff_data.loc[:,('Sag', 'Axe', 'Cor')])

        # Apply clustering on the distance matrix
        # eps: The maximum distance between two samples for one to be considered as in the neighborhood of the other. This is not a maximum bound on the distances of points within a cluster. This is the most important DBSCAN parameter to choose appropriately for your data set and distance function.
        # min_samples: The number of samples (or total weight) in a neighborhood for a point to be considered as a core point. This includes the point itself.
        clustering = DBSCAN(eps=50, min_samples=5).fit(dm)

        # Assign cluster labels to use as edit groups later
        # cur_diff_data.loc[:,('cluster')] = list(clustering.labels_)
        cur_diff_data.insert(0, "cluster", list(clustering.labels_), True)

        # Roll up by edit group (i.e. cluster) summarizing the beginning and the end of each cluster in each orientation
        
        if len(cur_diff_data['cluster'].unique())>1:
            out_grouped = cur_diff_data.groupby('cluster')
            out_grouped = out_grouped.agg(min_sag=('Sag', min),
                       max_sag=('Sag', max),
                       min_axe=('Axe', min),
                       max_axe=('Axe', max),
                       min_cor=('Cor', min),
                       max_cor=('Cor', max),
                       num_vox=('Cor', 'count')).reset_index()
        else:
            out_grouped = pd.DataFrame(data = {'cluster':cur_diff_data.loc[:, 'cluster'].unique()[0],
                                               'min_sag': min(cur_diff_data.Sag),
                                              'max_sag': max(cur_diff_data.Sag),
                                              'min_axe': min(cur_diff_data.Axe),
                                              'max_axe': max(cur_diff_data.Axe),
                                              'min_cor': min(cur_diff_data.Cor),
                                              'max_cor': max(cur_diff_data.Cor),
                                              'num_vox': cur_diff_data.shape[0]},
                                      index=[0])
        
        
        out_grouped.loc[:,('vol')] = [vol]*out_grouped.shape[0]
        out_grouped.loc[:,('action')] = [cur_action]*out_grouped.shape[0]
        out = out.append(out_grouped)
    
    return(out)

def get_bm_edits(editp, uneditp, vol, subnum, savediff, diffp):
    
    if vol == "brainmask":
        bm_diff_data = get_diff_data(editp=editp, uneditp=uneditp, vol=vol, subnum=subnum, savediff=savediff, diffp=diffp)
        bm_edits = summarize_edits(bm_diff_data, vol)
    else:
        print("Incorrect vol specification!")
        bm_edits = pd.DataFrame()

    return bm_edits

def get_wm_edits(editp, uneditp, vol, subnum, savediff, diffp):
    
    if vol == "wm":
        wm_diff_data = get_diff_data(editp=editp, uneditp=uneditp, vol=vol, subnum=subnum, savediff=savediff, diffp=diffp)
        wm_edits = summarize_edits(wm_diff_data, vol)
    else:
        print("Incorrect vol specification!")
        wm_edits = pd.DataFrame()

    return wm_edits

def get_cp_edits(editp, uneditp, vol, subnum, savediff=False, diffp=None):
    
    # Check if control points exist
    cp_fname = path.join(editp,'tmp/%s')%(vol)
    if path.isfile(cp_fname):
        # Read in control.dat
        cps = []
        with open(cp_fname, 'r') as fd:
            for line in fd:
                line = line.strip()
                vals = line.split(' ')
                if len(vals) == 3:
                    cps.append(vals)
        cps = np.array(cps, dtype=np.float)
        cps = pd.DataFrame(cps).rename(columns={0: "Sag", 1: "Cor", 2: "Axe"})

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
        # -12 = Cor - 128; Cor = 128 + (-12) = 116
        # 24 = 128 - Axe; Axe = 128 - (24) = 104
        
        cps['Sag'] = 128 - cps['Sag']
        cps['Cor'] = 128 + cps['Cor']
        cps['Axe'] = 128 - cps['Axe']
        
        # Reshape to make it have the same columns as the bm/wm edits for appending later
        # cluster	min_sag	max_sag	min_axe	max_axe	min_cor	max_cor	num_vox	vol	action
        out = pd.DataFrame(data={"cluster": np.nan,
                                "min_sag": cps['Sag'],
                                "max_sag": cps['Sag'],
                                "min_axe": cps['Axe'],
                                "max_axe": cps['Axe'],
                                "min_cor": cps['Cor'],
                                "max_cor": cps['Cor'],
                                "num_vox": 1,
                                "vol": "cp",
                                "action": "added cp"})
    
    else:
        print("No control points found for %s using vol = %s"%(subnum, vol))
        out = pd.DataFrame()
    
    # Always returns something to avoid exceptions when looping through and appending all types of edits below
    # (might not be necessary)
    return out

def get_bfs_edits(editp, uneditp, vol, subnum, savediff, diffp):
    
    if path.isfile(path.join(editp,'mri/%s.manedit.mgz')%(vol)):
        bfs_diff_data = get_diff_data(editp=editp, uneditp=uneditp, vol=vol, subnum=subnum, savediff=savediff, diffp=diffp)
        bfs_edits = summarize_edits(bfs_diff_data, vol)
    else: 
        print("No brain.finalsurfs.manedit.mgz file found for %s"%(subnum))
        bfs_edits = pd.DataFrame()

    return bfs_edits