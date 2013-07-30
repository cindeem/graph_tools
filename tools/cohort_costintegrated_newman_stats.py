
""" for a cohort of subjects
Across all:
create cohort array, check for NA data and drop subjects or nodes
Save out new files if necessary (as npy files)

in parallel:
1. for a range of costs
2. calc newman partition and binary matrix graph stats
3. save to subject specific directories
"""
import os, sys
import re
from glob import glob
import numpy as np
import pandas


def get_txt_files(datadir, pattern = 'PIDN*.txt'):
    """ given data directory and patternt o match, globs files and 
    returns a sorted list"""
    globstr = os.path.join(datadir, pattern)
    result = sorted(glob(globstr))
    return result

def get_subid(instr, pattern):
    """given a string and pattern , return re result, else raise 
    an error"""
    m = re.search(pattern, instr)
    try:
        return m.group()
    except:
        raise ValueError('%s not found in %s'%(pattern, instr))

    
def files_to_panel(file_list):
    """ given a list of txt files, create a 3D pandas Data Panel
    """
    subdict = {}
    for infile in file_list:
        subid = get_subid(infile, 'PIDN_[0-9]{5}')
        tmpdat = np.loadtxt(infile)
        subdict.update({subid: tmpdat})
    panel = pandas.Panel(subdict)
    return subdict

def qa_data(panel):
    """ check 3D data in panel to find nan values across cohort"""
    
    nsub, nnode, _ = panel.shape
    allmask = np.zeros((nsub, nnode))
    for val, subid in enumerate(panel):
        print subid
        tmp = panel[subid]
        ind = util.triu_indices_from(tmp)
        mask = ~np.isnan(tmp.values[ind])
        allmask[val] = mask.astype(int)

