import os
import nibabel as ni
from scipy.ndimage.measurements import center_of_mass
import numpy as np
import pickle

def voxel_centerofmass(dat, label_dat, label):
    """calculate center of mass (in voxel coords)
    given inputs
    
    Parameters
    ----------
    dat : np.array 
        "weight" of data, all ones for binary image
    label_dat : nd.array
        same shape as dat, has labelled regions
    label : int
        value of label in label_dat to use to calc
        center of mass
    """
    if not dat.shape == label_dat.shape:
        raise IOError('Shape mismatch dat: %s, label: %s'%(dat.shape,
                                                           label_dat.shape))
    if not label in set(label_dat.flatten()):
        raise IOError('label %d not in label_dat'%label)

    return center_of_mass(dat, label_dat, label)

def from_label_img(infile):
    """ parse infile data and create binary dat and labels to be used
    to calc center of mass
    Returns
    -------
    aff : np.array
        affine transform mapping voxels to world space
    dat : np.array
        array of label data where each label represents a unique region
    mask : np.array
        array of ones the same shape as dat
    """
    img = ni.load(infile)
    mask = np.ones(img.get_shape())
    return img.get_affine(), img.get_data(), mask

def find_labels(labeldat):
    """given an array with labels (ints)
    return unique labels, (omits 0)omits"""
    unique = set(labeldat.flatten())
    return [x for x in unique if not x == 0]



def world_coords(affine, voxel_coords):
    """calcs world coords given affine and voxel coords"""
    coords = np.ones(len(voxel_coords) + 1)
    coords[:-1] = voxel_coords
    world_coords = np.dot(affine, coords)
    return world_coords[:-1]


if __name__ == '__main__':

    infile  ='/home/jagust/cindeem/TEMPLATES/aal.nii.nii.nii.gz'
    aff, dat, mask = from_label_img(infile)
    labels = find_labels(dat)
    region_centers = {}
    for label in labels:
        vox_coords = voxel_centerofmass(mask, dat, label)
        w_coords = world_coords(aff, vox_coords)
        region_centers[label] = w_coords
    with open('aal_coords' , 'w+') as fid:
        pickle.dump(region_centers, fid)






