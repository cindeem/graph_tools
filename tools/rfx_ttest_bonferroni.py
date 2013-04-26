import sys, os
import argparse
import time
import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0, '/home/jagust/graph/scripts/graph_tools/graph_tools')
import cohort_graph_tools as cgt

"""
citation: He Y, Wang J, Wang L, Chen ZJ, Yan C, et al. (2009) Uncovering
Intrinsic Modular Organization of Spontaneous Brain Activity in Humans. PLoS
ONE 4(4): e5226. doi:10.1371/journal.pone.0005226


A random-effect one-sample t test was further performed on these correlation matrices in an element-by-element manner to
obtain the significance level (i.e. P value) of each inter-regional correlation across the subjects. Finally, the
P-value matrix was thresholded by using a conservative Bonferroni-corrected P value (P = 0.001) to reduce the chance of
false positives, which resulted in a binarized matrix (sparsity = 8.41%) that captured the functional connectivity
backbone underlying the topological organization of spontaneous human brain activity at a time domain (Figure 1B).

try:
    indata = sys.argv[1]
except:
    indata = '/home/jagust/graph/data/spm_189/All_unmasked_data_fisher.npy'
alpha = .01
"""
def main(indata, alpha, save = False):
    try:
        data = np.load(indata)
    except:
        raise IOError('unable to open %s, is this a *.npy file?'%indata)

    mask, tvals, sparsity = cgt.gen_bonferroni_corrected_graph(indata, alpha=alpha)
    if save:
        pth, fname = os.path.split(indata)
        nme,ext = os.path.splitext(fname)
        curr_time = time.strftime('%Y-%m-%d-%H-%M')
        newf = os.path.join(pth, 'bonferroni_cohort_mask_' + nme + curr_time + ext)
        np.save(newf, mask)
        datf = newf.replace('cohort_mask', 'cohort_correl')
        avg_data = data.squeeze().mean(0)
        avg_data[mask == False] = 0
        np.save(datf, avg_data)
        print 'wrote %s'%(newf)

    plt.subplot(121)
    plt.imshow(mask , cmap= plt.cm.binary_r, interpolation='none')
    plt.subplot(122)
    plt.imshow(avg_data, cmap = 'jet', interpolation = 'none')
    plt.colorbar()
    plt.title('Group Network corrected at %2s (sparsity = %2.2f)'%(alpha,
                                                                   sparsity * 100))
    plt.show()


if __name__ == '__main__':
    
    # create the parser
    parser = argparse.ArgumentParser(
                        description='Generate thresholded cohort network map')

    # add the arguments
    parser.add_argument(
                        'infile', type=str, nargs=1,
                        help="""npy file containing cohort 
                                connectivity matricies""")
    parser.add_argument('--alpha', type=float,
                        default=.001, dest='alpha',
                        help='Threshold for Bonferroni correction (default .001)')

    parser.add_argument('--save_thr_matrix', default = False, 
                        action = 'store_true',
                        help = 'save bonferroni corrected binary matrix')
    if len(sys.argv) == 1:
        parser.print_help()
    else:
        args = parser.parse_args()
        main(args.infile[0], args.alpha, save = args.save_thr_matrix)

