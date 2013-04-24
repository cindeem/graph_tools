


""" Use boot strapping to take a cohort of subjects with individual network adjacent matricies, and with resampling,
estimate a cohort mean, and significance for connection between each set of regions

will return a estimate to use for a weighted graph, along with a thresholded mask for calculating modularity.

"""
import os
import numpy as np
from scipy.stats import ttest_1samp

datdir = '/home/jagust/UCSF/Manja_Lehmann/ICN/rsfMRI_in_AD/data/extended_cohort_30controls_49AD/graph_theory/Controls'
correlations = os.path.join(datdir, 'Controls_unmasked_data_fisher.npy')
mat = np.load(correlations)

np.random.seed()
resampled = np.zeros((1000, 90, 90))
resamp_pvals = np.zeros((1000, 90, 90))
for i in np.arange(1000):
    indicies = np.random.randint(0,29,30)
    mean = mat[0,indicies,:,:].mean(axis=0)
    tval, pval = ttest_1samp(mat[0,indicies,:,:], 0, axis=0)
    resampled[i, : :] = mean
    resamp_pvals[i,:,:] = pval

boot_correl = resampled.mean(0)
boot_pvals = resamp_pvals.mean(0)

mask = np.zeros(boot_correl.shape)
mask[boot_pvals <=.001] = 1
mask[boot_correl < 0] = 0



