
from glob import glob
import os
import sys
import argparse

def write_mods(indir, gstr='*cost*.txt'):
    """collates mods written to individual files in a directory
    into one csv spreadsheet"""
    
    globstr = os.path.join(indir, gstr)
    allmods = sorted(glob(globstr))

    outfile = os.path.join(indir, 'all_modularity.csv')

    with open(outfile, 'w+') as fid:
        fid.write('SUBID,modularity,\n')
        for item in allmods:
            jnk = open(item).read()
            id, mod, _ = jnk.split(',')
            fid.write('%s,%s,\n'%(id, mod))
    print 'wrote', outfile


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Collate mod files to one csv')
    parser.add_argument('indir',
                        help='directory holding mod files')
    parser.add_argument('--gstr', dest='gstr',
                        default = '*cost*.txt',
                        help='file specific glob (default *cost*.txt)')

    if len(sys.argv) < 2:
        parser.print_help()
    else:
        args = parser.parse_args()
        print args
        write_mods(args.indir, gstr = args.gstr)
