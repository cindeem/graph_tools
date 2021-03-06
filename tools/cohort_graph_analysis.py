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

"""
citation: He Y, Wang J, Wang L, Chen ZJ, Yan C, et al. (2009) Uncovering
Intrinsic Modular Organization of Spontaneous Brain Activity in Humans. PLoS
ONE 4(4): e5226. doi:10.1371/journal.pone.0005226


A random-effect one-sample t test was further performed on these correlation matrices in an element-by-element manner to
obtain the significance level (i.e. P value) of each inter-regional correlation across the subjects. Finally, the
P-value matrix was thresholded by using a conservative Bonferroni-corrected P value (P = 0.001) to reduce the chance of
false positives, which resulted in a binarized matrix (sparsity = 8.41%) that captured the functional connectivity
backbone underlying the topological organization of spontaneous human brain activity at a time domain (Figure 1B).
"""

def gen_bonferroni_corrected_graph(indata, alpha = .01):
    """ given indata (4D array subjects X nnodes X nnodes)
    and alpha (rejection level for bonferroni correction)
    calculate significant connections across cohort
    return thresholded t-value map and binary map
    """

    data = np.load(indata)
    data = data.squeeze() # remove singular dims (eg we only have 1 block)
    ind = util.triu_indices(data.shape[-1], 1)# ind of upper tri minus diag
    lowind = util.tril_indices(data.shape[-1], -1) # ind of lower triag minusdiag
    tvals, pvals = ss.ttest_1samp(data, 0, axis = 0) ## rfx ttest across subjects
    ## for each region
    (reject, pvals_corrected, 
    alphacSidak, alphacBonf) = sms.multipletests(pvals[ind], 
                                                 alpha = alpha, 
                                                 method='bonferroni')
    (lreject, lpvals_corrected, 
     lalphacSidak, lalphacBonf) = sms.multipletests(pvals[lowind], 
                                                    alpha = alpha, 
                                                    method='bonferroni')
    print 'sparsity', reject.astype(float).sum() / reject.shape[0]

    mask = np.zeros(tvals.shape, dtype = np.bool)
    mask[ind] = reject
    mask[lowind] = lreject
    mask[tvals < 0] = False
    tvals[mask == False] = 0
    return mask, tvals


def boot_corrected_graph(indata, nreps = 1000, alpha = .01):
    """Uses bootstrapping to generate mask, correct for
    multiple comparisons and give robust estimate for 
    cohort"""
    mat = np.load(indata)
    np.random.seed()
    nb, nsub, nnode, _ = mat.shape
    resampled = np.zeros((nreps, nnode, nnode))
    resamp_pvals = np.zeros((nreps, nnode, nnode))
    for i in np.arange(nreps):
        indicies = np.random.randint(0,nsub-1, nsub)
        mean = mat[0,indicies,:,:].mean(axis=0)
        tval, pval = ss.ttest_1samp(mat[0,indicies,:,:], 0, axis=0)
        resampled[i, : :] = mean
        resamp_pvals[i,:,:] = pval

    boot_correl = resampled.mean(0)
    # significance is at .001, correction at alpha

    mask = resamp_pvals <= .001
    mask = mask.sum(0) / np.float(nreps)
    mask[mask < 1-alpha] = 0
    mask[mask >= 1-alpha] = 1
    return mask, boot_correl




def plot_cohort_map(mask, outdir, other = None):
    """ plots N X N connectivity map in outdir
    if other is another N X N map, will plot side by side"""
    
    
    if other is not None:
        plt.subplot(122)
        
        plt.imshow(other, cmap= plt.cm.binary_r, interpolation='none')
        plt.title('Connectivity')
        plt.colorbar()
        plt.subplot(121)

    plt.imshow(mask, cmap= plt.cm.binary_r, interpolation='none')
    plt.title('Connectivity')
    plt.colorbar()
    outf = os.path.join(outdir, 'connectivity_plot.png')
    plt.savefig(outf)


def calc_modularity(inmat,  ideal_cost = 0.1):
    """
    Once we have a thresholded matrix, we can use simulated annealing to 
    calculate the modularity of the graph, 
    datgrap.index shows us the components that make up the 
    distinct modules in the network (note ids start at 0, aal count from 1)
    """

    G = nx.Graph(weighted=False)
    G = nx.from_numpy_matrix(inmat, G)
    ## fixed config parameters for simulated annealing
    temperature = 0.1
    temp_scaling = 0.9995
    tmin = 1e-4
    ideal_cost = 0.1
    (datgraph, 
     modularity_value) = md.simulated_annealing(G,
                                                temperature = temperature,
                                                temp_scaling = temp_scaling,
                                                tmin = tmin,
                                                extra_info = False,
                                                debug = False)
    print modularity_value
    return datgraph, modularity_value

def save_modularity(datgraph, mod_val, outdir):
    """ saves index values form datgraph and moulatiry value to txt file"""
    outf = os.path.join(outdir, 'cohort_modularity_index.pkl')
    pickle.dump(datgraph.index, open(outf, 'w+'))
    print 'wrote', outf
    outf = outf.replace('index.pkl', 'value.txt')
    with open (outf, 'w+') as fid:
        fid.write('%f'%mod_val)
    print 'wrote', outf


def load_modularity(datdir, fname = 'cohort_modularity_index.pkl'):
    """ load pre computed module indicies from pickled file"""
    infile = os.path.join(datdir, fname)
    mod_index = pickle.load(open(infile))
    return mod_index


def parse_aal(dim):
    """ load aal coords from pickled aal_coords file
    parse, and return coords in xy, or yz plane
    parse aal.txt to get label names
    dim can be 'xy' or 'yz' """

    aal_coords = pickle.load(open('aal_coords'))
    aal_90 = {}
    for label, val in sorted(aal_coords.items()):
        if label < 91:
            if dim == 'xy':
                aal_90.update({label-1:val[:-1]})
            elif dim == 'yz':
                aal_90.update({label-1: val[1:]})
            else:
                raise IOError('%s is not \'xy\' or \'yz\', invalid dim'%dim)
    tmp = np.loadtxt('aal.txt', dtype = str)
    labels = tmp[:90,1]
    labels = [x for x in labels]
    return aal_90, labels

def gen_node_colors(modules):
    """gen_node_colors given a dict of nodes calulated but
    running cal_modularity, assign a color to each member of a 
    node group"""
    newlist = {}
    n_mods = len(modules) # number of distinct modules in network
    for val, node_set in enumerate(modules.values()):
        color = plt.cm.jet(1. * val / n_mods)
        tmp = dict.fromkeys(node_set, color)
        newlist.update(tmp)
    return newlist



def plot_weighted_graph(G, outdir, dim = 'xy', node_colors = None, title='graph'):
    """docstring for plot_weighted_graph"""
    aal_90, labels = parse_aal(dim)
    if node_colors is None:
        nxplot.draw_graph(G, layout = aal_90, labels = labels)
    else:
        nxplot.draw_graph(G, layout = aal_90, labels = labels,
                          node_colors = node_colors,
                          title = title)
    plt.show()
    #plt.gca()
    #plt.savefig(os.path.join(outdir, 'cohort_graph_%s.png'%(dim)), dpi=600)


def save_graph_metrics(G, outdir):
    summary = metrics.graph_summary(G)
    outf = os.path.join(outdir, 'cohort_graph_summary.json')
    json.dump(summary, open(outf, 'w+'))
    print summary
    print 'saved',  outf

if __name__ == '__main__':
             
    try:
        indata = sys.argv[1]
    except:
        indata = '/home/jagust/graph/data/spm_189/All_unmasked_data_fisher.npy'



