"""
Script for calculating k-mer enrichment.
@author: Rachael Cox <rachaelcox@utexas.edu>
"""

__author__ = "Rachael Cox (rachaelcox@utexas.edu)"
__version__ = "0.1.0"
__license__ = "MIT"

import argparse
import math
import pandas as pd
import time
import re
from scipy.stats import norm
from statsmodels.stats.multitest import multipletests

def merge_data(exp_file, bg_file):
    exp_df = pd.read_csv(exp_file).rename(columns={'count': 'exp_count'})
    bg_df = pd.read_csv(bg_file).rename(columns={'count': 'bg_count'})
    exp_label = re.search('(?<=_)\w\w(?=_)', exp_file)[0]
    join_cols = ['node1','node2']
    join_df = exp_df.merge(bg_df, how='outer', left_on=join_cols, right_on=join_cols).reset_index(drop=True)
    min_exp = join_df['exp_count'].min(skipna=True)
    min_bg = join_df['bg_count'].min(skipna=True)
    join_df['exp_count'] = join_df['exp_count'].fillna(min_exp).astype(int)
    join_df['bg_count'] = join_df['bg_count'].fillna(min_bg).astype(int)
    join_df['exp_norm'] = join_df['exp_count']/join_df['exp_count'].sum()
    join_df['bg_norm'] = join_df['bg_count']/join_df['bg_count'].sum()
    join_df['exp'] = exp_label
    join_df['node_edge'] = join_df['node1']+join_df['node2'].str.rstrip().str[-1]
    return(join_df)

def get_oligo_set(truth_file):
    df = pd.read_csv(truth_file)
    df.columns = map(str.lower, df.columns)
    df['node_edge'] = df['node1']+df['node2'].str.rstrip().str[-1]
    truth_set = df['node_edge'].unique()
    truth_set = set(truth_set)
    return(truth_set)

def get_library_set(lib_file):
    df = pd.read_csv(lib_file)
    df.columns = map(str.lower, df.columns)
    df['node_edge'] = df['node1']+df['node2'].str.rstrip().str[-1]
    lib_set = df['node_edge'].unique()
    lib_set = set(lib_set)
    return(lib_set)

def calc_fc(exp, bg):
    if bg > exp:
        fc = -bg/exp
        log2fc = -(math.log2(abs(fc)))
    elif bg <= exp:
        fc = exp/bg
        log2fc = math.log2(fc)
    return(fc, log2fc)

def calc_pseudo_fc(df):
    # pseudocounts
    df['exp_pc'] = df['exp_norm']
    df['bg_pc'] = df['bg_norm']
    df['f0_exp_pc'] = df['exp_pc']/df['exp_pc'].sum()
    df['f0_bg_pc'] = df['bg_pc']/df['bg_pc'].sum()
    # fold change cols
    df[['fc','log2fc']] = [calc_fc(i,j) for i,j in zip(df['f0_exp_pc'], df['f0_bg_pc'])]
    # drop extra cols
    df.reset_index(drop=True, inplace=True)
    df.drop(['exp_pc','bg_pc','f0_exp_pc','f0_bg_pc'], axis=1, inplace=True)
    df.sort_values(by='log2fc', ascending=False, inplace=True)
    return(df)

def calc_apex_zscore(f0_exp, f0_bg, f1, exp_sum, bg_sum):
    apex_zscore = (
        (f0_exp - f0_bg)/math.sqrt(
        (((f1*(1-f1))/exp_sum) + ((f1*(1-f1))/bg_sum)))
                  )
    return(apex_zscore)

def calc_pval(zscore, one_side=True):
    p_value = norm.cdf(zscore)
    # for one-sided test, subtract the obtained p-value from 1
    if one_side:
        p_val = 1 - p_value
    return(p_val)

def main(args):

    t0 = time.time()
    print('Reading in files ...')
    lib_kmers = get_library_set(args.lib_file)
    oligo_kmers = get_oligo_set(args.oligo_file)
    print('Normalizing data ...')
    df = merge_data(args.exp_file, args.bg_file)
    exp_sum = df['exp_count'].sum()
    bg_sum = df['bg_count'].sum()
    df['label'] = (['oligo match' if edge in oligo_kmers \
                    else 'library match' if edge in lib_kmers \
                    else 'no match' for edge in df['node_edge']])
    print('Calculating enrichment ...')
    df = calc_pseudo_fc(df)
    df['f0_exp'] = df['exp_count']/exp_sum
    df['f0_bg'] = df['bg_count']/bg_sum
    # apex zscore calculation
    df['f1'] = (df['exp_count']+df['bg_count'])/(df['exp_count'].sum()+df['bg_count'].sum())
    df['apex_zscore'] = [calc_apex_zscore(i, j, k, exp_sum, bg_sum) for i, j, k in zip(df.f0_exp, df.f0_bg, df.f1)]
    df['apex_pval'] = [calc_pval(z) for z in df.apex_zscore]
    bonf_pvals = multipletests(df['apex_pval'], method='bonferroni')[1]
    bh_pvals = multipletests(df['apex_pval'], method='fdr_bh')[1]
    df['apex_pval_bonferroni'] = bonf_pvals
    df['apex_pval_bh'] = bh_pvals
    df.sort_values('apex_pval_bh', ascending=True, inplace=True)
    print(df)
    print('Writing out results ...')
    df.to_csv(args.outfile, index=False)
    print('Done!')
    print(f'Total run time: {round((time.time()-t0)/60, 2)} minutes.')

if __name__ == "__main__":
    """ This is executed when run from the command line """
    parser = argparse.ArgumentParser()

    parser.add_argument("--exp_file", action="store", required=False,
                        help="(Required) Path to file with experimental k-mer counts.")

    parser.add_argument("--bg_file", action="store", required=False,
                        help="(Required) Path to file with background k-mer counts.")

    parser.add_argument("--lib_file", action="store", required=False,
                        help="(Required) Path to file with library k-mers.")

    parser.add_argument("--oligo_file", action="store", required=False,
                        help="(Required) Path to file with oligonucleotide bait k-mers (complementary).")

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
    main(args)