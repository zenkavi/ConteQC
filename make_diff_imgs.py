#!/Users/zeynepenkavi/anaconda/envs/py37/bin/python

from argparse import ArgumentParser
import os
from os import path
import nibabel as nib
import numpy as np

parser = ArgumentParser()
parser.add_argument("--subnum")
#Conte-Freesurfer/FS7/pass0_edits
parser.add_argument("--editp", default='/Users/zeynepenkavi/Downloads/ConteQC/pass0_edits')
#Conte-Freesurfer/FS7/pass0_defaced
parser.add_argument("--uneditp", default='/Users/zeynepenkavi/Downloads/ConteQC/pass0')

args = parser.parse_args()
subnum = args.subnum
editp = args.editp
uneditp = args.uneditp

vols = [i for i in os.listdir(path.join(editp, subnum)) if i.endswith('.mgz')]

for vol in vols:

    edit_img = nib.load(path.join(editp,'%s/%s')%(subnum, vol))
    edit_data = edit_img.get_fdata()
    unedit_img = nib.load(path.join(uneditp,'%s/mri/%s')%(subnum, vol))
    unedit_data = unedit_img.get_fdata()
    diff_data = unedit_data - edit_data


    diff_img = nib.MGHImage(diff_data.astype(np.int32), edit_img.affine)
    nib.save(diff_img, os.path.join(editp, subnum, 'diff_%s')%(vol))