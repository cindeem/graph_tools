#!/usr/bin/env python

import os, sys
import time
import numpy as np
import argparse
sys.path.insert(0, '/home/jagust/graph/scripts/graph_tools/graph_tools')
import cohort_graph_tools as cgt


def plot_graph(inmat, inmod, dim, title):
    
    inmat = np.load(inmat)
    mod_index = cgt.load_modularity(inmod)
    # make graph
    G = cgt.util.mat2graph(inmat)
    node_colors = cgt.gen_node_colors(mod_index)
    cgt.plot_weighted_graph(G, dim=dim, node_colors = node_colors,
                            title = title)


if __name__ == '__main__':

    parser = argparse.ArgumentParser(
            description = 'Generate thresholded cohort map')

    parser.add_argument('inmat', type=str, nargs = 1,
                        help = "numpy file holding cohort connectivity\
                        matricies adjacency matrix\
                        eg bootstrap_thresholded_correl_1000_pval0.001_\
                        alpha0.01_2013-04-26-10-17_.npy")
    parser.add_argument('inpkl', type = str, nargs = 1,
                        help = 'pickle file holding network modules')

    parser.add_argument('-dim', type = str, default = 'xy', dest = 'dim', 
                        help =  'plane to draw network, xy or yz')
    parser.add_argument('-title', type = str, default = 'graph', 
                        help = 'optional title for figure')

    if len(sys.argv) == 1:
        parser.print_help()
    else:
        args = parser.parse_args()
        plot_graph(args.inmat[0], args.inpkl[0], args.dim, args.title)

