"""
Script for combining k-mer data across conditions and replicates.
@author: Rachael Cox <rachaelcox@utexas.edu>
"""

__author__ = "Rachael Cox (rachaelcox@utexas.edu)"
__version__ = "0.1.0"
__license__ = "MIT"

import argparse
import pandas as pd
import time
import re
import glob 
import pickle

''' Functions '''
def get_data(data_dir, k):
    
    files = glob.glob(f'{data_dir}/*{k}mer*')
    print(f'Files detected:')
    print(files)
    
    print(f'Reading in data ...')
    df_list = []
    for i in files:
        name = re.search('([^\/]+$)', i)[0]
        meta_data = name.split('_')
        if len(meta_data) != 6:
            substrate = meta_data[0]
            condition = meta_data[1]
            rep = meta_data[2].replace('rep','')
            k = meta_data[3].replace('mers','')
        else:
            substrate = '_'.join([meta_data[0], meta_data[1]])
            condition = meta_data[2]
            rep = meta_data[3].replace('rep','')
            k = meta_data[4].replace('mers','')
        df = pd.read_csv(i)
        df[['substrate','condition','rep','k']] = [substrate, condition, rep, k]
        df_list.append(df)

    return(df_list)

def format_long(df_list, cutoff=None):

    df_long = pd.concat(df_list).reset_index(drop=True)
    if cutoff:
        df_long = df_long[df_long['count'] >= cutoff]
    return(df_long)

def format_wide(df_long):
    
    print('Converting data to wide format ... ')
    # pivot count data to long format
    df_wide = df_long.pivot(index=['node1','node2'], columns=['condition','rep'], values='count').reset_index()
    # flatten multi-index columns
    df_wide.columns = df_wide.columns.map('_'.join).str.strip('_')
    # reorder columns lexicographically
    df_wide = df_wide.reindex(sorted(df_wide.columns), axis=1)
    # move node columns back to front
    df_wide = df_wide.set_index(['node1','node2']).reset_index()
    # convert floats back to int
    num_cols = df_wide.columns[df_wide.dtypes.eq('float64')]
    df_wide[num_cols] = df_wide[num_cols].fillna(0).astype(int)

    return(df_wide)

''' Main '''

def main():
    
    data_list = get_data(args.data_dir, args.k)
    df_long = format_long(data_list, args.cutoff)
    df_wide = format_wide(df_long)

    # write out results
    print(f'Writing long table to {args.out_prefix}_long.csv ...')
    df_long.to_csv(f'{args.out_prefix}_long.csv', index=False)
    print(f'Writing wide table to {args.out_prefix}_wide.csv ...')
    df_wide.to_csv(f'{args.out_prefix}_wide.csv', index=False)
    

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Script for combining k-mer data across conditions and replicates.")
    parser.add_argument("--data_dir", action="store", required=True,
                                        help="Directory of k-mer count files")
    parser.add_argument("--k", action="store", required=True,
                                        help="K-mer size for target files")
    parser.add_argument("--out_prefix", action="store", required=True,
                                        help="Path to result directory + prefix of outfile name")
    parser.add_argument("--cutoff", action="store", required=False, default=None, type=int,
                                        help="(Optional) Minimum k-mer count to keep")
    args = parser.parse_args()
    main()