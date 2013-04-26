#!/usr/bin/env python

""" Use boot strapping to take a cohort of subjects with individual network adjacent matricies, and with resampling,
estimate a cohort mean, and significance for connection between each set of regions

will return a estimate to use for a weighted graph, along with a thresholded mask for calculating modularity.

"""
import os, sys
import time
import numpy as np
import argparse
sys.path.insert(0, '/home/jagust/graph/scripts/graph_tools/graph_tools')
import cohort_graph_tools as cgt

def boot_main(indata, outd, nperms=1000, pval=.001, alpha=.01):
    """docstring for boot_main"""
    mask, boot_correl = cgt.boot_corrected_graph(indata,
                                                 nreps=nperms,
                                                 pthr = pval,
                                                 alpha = alpha)
    
    currtime = time.strftime('%Y-%m-%d-%H-%M')
    fname = '_'.join(['bootstrap_thresholded', '%04d'%nperms, 
                      'pval%s'%pval, 'alpha%s'%alpha,
                      currtime, '.npy'])
    outmask = os.path.join(outd, 
                           fname.replace('thresholded', 'thresholded_mask'))
    np.save(outmask, mask)
    outdat = os.path.join(outd, 
                          fname.replace('thresholded', 'thresholded_correl'))
    np.save(outdat, boot_correl)
    cgt.plot_chohort_map(mask, other = boot_correl)
    print 'saved', outmask, outdat

    

if __name__ == '__main__':

    parser = argparse.ArgumentParser(
            description = 'Generate thresholded cohort map')

    parser.add_argument('infile', type=str, nargs = 1,
                        help = "numpy file holding cohort connectivity\
                                matricies")

    parser.add_argument('-o', type=str, dest='outd',
                        help='Directory to save new matrix\
                              (default same as input)')
    parser.add_argument('-nperm', type = int, dest = 'nperm', default = 1000,
                        help = 'Number of bootstrap permutations\
                                (default 1000)')
    parser.add_argument('-pval', type = float, dest = 'pval', default = .001,
                        help='pval considered significant\
                             (dafault .001)')
    parser.add_argument('-alpha', type=float, dest='alpha', default = .01,
                        help='alpha level for multiple corrections (default .01)')

    if len(sys.argv) == 1:
        parser.print_help()
    else:
        args = parser.parse_args()
        if args.outd is None:
            args.outd, _  = os.path.split(args.infile)

        print args
        boot_main(args.infile[0], args.outd, 
                  nperms=args.nperm, 
                  pval= args.pval, 
                  alpha= args.alpha)


