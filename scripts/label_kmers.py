"""
Script for labeling k-mer matches.
@author: Rachael Cox <rachaelcox@utexas.edu>
"""

__author__ = "Rachael Cox (rachaelcox@utexas.edu)"
__version__ = "0.1.0"
__license__ = "MIT"

import argparse
import pandas as pd
import time
from Bio.Seq import Seq

''' Functions '''
def get_oligo_set(truth_file):
    df = pd.read_csv(truth_file)
    df.columns = map(str.lower, df.columns)
    df['edge'] = df['node1']+df['node2'].str.rstrip().str[-1]
    truth_set = df['edge'].unique()
    truth_set = set(truth_set)
    return(truth_set)

def get_library_set(lib_file):
    df = pd.read_csv(lib_file)
    df.columns = map(str.lower, df.columns)
    df['edge'] = df['node1']+df['node2'].str.rstrip().str[-1]
    lib_set = df['edge'].unique()
    lib_set = set(lib_set)
    return(lib_set)

''' Main '''
def main():

    t0 = time.time()
    print('Reading in library/bait files ...')
    lib_kmers = get_library_set(args.lib_file)
    comp_kmers = get_oligo_set(args.bait_file)
    direct_kmers = []
    for i in comp_kmers:
        i = Seq(i)
        oligo_kmer = i.complement()
        direct_kmers.append(str(oligo_kmer))
    direct_kmers = set(direct_kmers)
    print(f'Reading in {args.results_file} ...')
    df = pd.read_csv(args.results_file)
    print('Labeling k-mers ...')
    df['label'] = (['complement oligo match' if edge in comp_kmers \
                    else 'direct oligo match' if edge in direct_kmers \
                    else 'library match' if edge in lib_kmers \
                    else 'no match' for edge in df['edge']])

    if args.remove_lib:
        print('Removing k-mers matching the library static regions ...')
        df = df[df['label'] != 'library match']
        if args.drop_label_col:
            print('Dropping label column ...')
            df = df.drop(['label'], axis=1)

    print('Writing out results ...')
    df.to_csv(args.outfile, index=False)
    print('Done!')
    print(f'Total run time: {round((time.time()-t0)/60, 2)} minutes.')

if __name__ == "__main__":
    """ This is executed when run from the command line """
    parser = argparse.ArgumentParser()

    parser.add_argument("--results_file", action="store", required=False,
                        help="(Required) Path to file with experimental k-mer counts/enrichment.")

    parser.add_argument("--lib_file", action="store", required=False,
                        help="(Required) Path to file with library k-mers.")

    parser.add_argument("--bait_file", action="store", required=False,
                        help="(Required) Path to file with oligonucleotide bait k-mers (complementary).")

    parser.add_argument("--remove_lib", action="store_true", required=False, default=False, help="(Optional) Remove library matches from output.")

    parser.add_argument("--drop_label_col", action="store_true", required=False, default=False, help="(Optional) Drop label column from output, i.e., as a follow up to '--remove_lib'")

    parser.add_argument("--outfile", action="store", required=False,
                        help="(Required) Outfile path/name (.csv)")

    # optional verbosity counter (eg. -v, -vv, -vvv, etc.)
    parser.add_argument(
        "-v",
        "--verbose",
        action="count",
        default=0,
        help="Verbosity (-v, -vv, etc)")

    # Specify output of "--version"
    parser.add_argument(
        "--version",
        action="version",
        version="%(prog)s (version {version})".format(version=__version__))

    args = parser.parse_args()
    main()