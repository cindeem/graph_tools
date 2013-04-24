import sys, os
import numpy as np
import statsmodels.stats.api as sms
import scipy.stats as ss
sys.path.insert(0,'/home/jagust/graph/scripts/brainx')
from brainx import (util, metrics)
from brainx import modularity as md
from brainx import nxplot
import networkx as nx
import matplotlib.pyplot as plt
import pickle
import json

#sys.path.insert(0, '/home/jagust/graph/scripts/graph_tools/graph_tools')
from .. import cohort_graph_tools as cgt

import numpy.testing as npt


def get_data(which=None):
    pth, _ = os.path.split(__file__)
    if which is None:
        datadir = os.path.join(pth, 'data')
    else:
        datadir = os.path.join(pth, 'data', which)
    return datadir

def test_sparsity():
    jnk = np.ones((10), dtype=bool)
    yield npt.assert_equal, cgt.calc_sparsity(jnk), 1
    jnk[:5] = 0
    yield npt.assert_equal, cgt.calc_sparsity(jnk), 0.5


def test_bonferroni():
    datadir = get_data('test_fisher_data.npy')
    mask, tvals, sparse = cgt.gen_bonferroni_corrected_graph(datadir)
    yield npt.assert_equal, np.diag(tvals), np.zeros((10))
    yield npt.assert_equal, np.diag(mask), np.zeros((10), dtype=bool)
    yield npt.assert_equal, tvals[0, :5].sum(),0.0
    yield npt.assert_equal, np.all(mask[0,:5]==False), True
    yield npt.assert_almost_equal, tvals.max(), 13.3856011

def test_parse_aal():
    aal_txt = get_data('aal.txt')
    aal_pickle = get_data('aal_coords')
    aal_90, labels = cgt.parse_aal(aal_pickle, aal_txt)
    yield npt.assert_almost_equal, aal_90[0], np.array([-38.6495705,  -5.6833250])
    yield npt.assert_equal, labels[0], 'Precentral_L'
    aal_90, labels = cgt.parse_aal(aal_pickle, aal_txt, dim='yz')
    yield npt.assert_almost_equal, aal_90[0], np.array([ -5.6833250,  50.9442038])
    yield npt.assert_equal, labels[0], 'Precentral_L'
    npt.assert_raises(IOError,  cgt.parse_aal, aal_pickle, aal_txt, dim='mm')


def test_gen_node_colors():
    """docstring for test_gen_node_colors"""
    modules = {0:[1,2], 1:[3,4]}
    node_colors = cgt.gen_node_colors(modules)
    yield npt.assert_equal, node_colors[1], plt.cm.jet(0)
    yield npt.assert_equal, node_colors.has_key(0), False

def test_boot_corrected_graph():
    datadir = get_data('test_fisher_data.npy')
    mask, boot_correl = cgt.boot_corrected_graph(datadir, nreps=10)
    yield npt.assert_equal, mask[0,:5].mean(), 0.0
    yield npt.assert_equal, mask[0, 5:].mean(), 1.0
    sparse = cgt.calc_sparsity(mask)
    yield npt.assert_equal, sparse, .75
