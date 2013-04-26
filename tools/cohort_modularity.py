import os, sys
import argparse
sys.path.insert(0, '/home/jagust/graph/scripts/graph_tools/graph_tools')
import cohort_graph_tools as cgt


def run_modularity(infile, outdir, ideal_cost = .01):
    datgraph, modval = cgt.calc_modularity(infile, ideal_cost = ideal_cost)
    pkl, txt = cgt.save_modularity(datgraph, modval, outdir)




if __name__ == '__main__':

    parser = argparse.ArgumentParser(
            description = 'Calculates modularity from binary\
                    adjacency matrix')
    parser.add_argument('infile', type=str, nargs = 1,
                        help= 'numpy file holding binary adjacency matrix\
                        eg: bootstrap_thresholded_mask_1000_pval0.001_alpha0.01_2013-04-26-10-15_.npy')
    parser.add_argument('-cost', type=float, dest = 'cost', default = .01,
                        help = 'ideal cost for running simulated annealing')
    parser.add_argument('-outd', type=str, dest = 'outd', 
                        help = 'out directoryto save files')

    if len(sys.argv) == 1:
        parser.print_help()
    else:
        args = parser.parse_args()
        if args.outd is None:
            args.outd, _ = os.path.split(args.infile[0])
        run_modularity(args.infile[0], args.outd, args.cost) 
