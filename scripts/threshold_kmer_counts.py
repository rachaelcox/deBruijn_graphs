"""
Script for thresholding de Bruijn graph k-mer edges by count.
@author: Rachael Cox <rachaelcox@utexas.edu>
"""

__author__ = "Rachael Cox (rachaelcox@utexas.edu)"
__version__ = "0.1.0"
__license__ = "MIT"

import argparse
import pandas as pd
import time

def read_data(infile, sep=','):
    print(f'Reading in {infile} ...')
    df = pd.read_csv(infile, sep=sep)
    print(f'Number of unique k-mers = {len(df)}')
    return(df)

def calc_quartiles(df):
    print(f'Calculating quartiles ...')
    Q1 = df['count'].quantile(0.25)
    Q3 = df['count'].quantile(0.75)
    IQR = Q3-Q1
    print(f'Q1={Q1}\tQ3={Q3}\tIQR={IQR}')

def calc_cutoff(df, cut):
    print(f'Calculating quantile cutoff for q={cut} ...')
    q = df['count'].quantile(cut)
    df_cut = df[df['count'] > q]
    n_removed = len(df) - len(df_cut)
    print(f'Count={q} for q={cut} ({n_removed} rows removed) ...')
    print(f'Final table:')
    print(df_cut)
    return(df_cut)

def main(args):
    
    t0=time.time()
    df = read_data(args.infile, args.sep)
    calc_quartiles(df)
    df_out = calc_cutoff(df, args.cutoff)
    print(f'Writing to {args.outfile} ...')
    df_out.to_csv(args.outfile, index=False)
    print(f'Total run time: {round((time.time()-t0)/60, 2)} minutes.')

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser()

    # file in
    parser.add_argument("-i", "--infile", help="(Required) Path to and name of infile.")

    # delimiter for infile
    parser.add_argument("-s", "--sep", default=',', help="(Optional) Infile delimiter; default=','.")

    # desired threshold
    parser.add_argument("-c", "--cutoff", default=0.05, type=float, help="(Optional) Quantile threshold; default=0.05.")

    # file out
    parser.add_argument("-o", "--outfile", help="(Required) Path to and name of outfile.")

    # Optional verbosity counter (eg. -v, -vv, -vvv, etc.)
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
    main(args)