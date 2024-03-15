import csv
import pandas as pd
import argparse
import pickle
from Bio import SeqIO
from itertools import groupby
from itertools import filterfalse
import numpy as np
from datetime import datetime as dt

# generate unique edge list with aggregate protein IDs

def map_kmer(kmer, dct):
    if kmer in dct.keys():
        return(dct[kmer])
    return(np.nan)

def unique(input_file, output_file, pickle_file, return_ids=False):
    
    # read input to dataframe & create new file name
    if not output_file:
        output_file = input_file.replace(".csv","")
        output_file = "{}_unique.csv".format(output_file)

    print(f"[{dt.now()}] Reading in de Bruijn graph from {input_file} ...")
    try:
        df = pd.read_csv(input_file)
    except:
        with open(input_file, 'rb') as handle:
            df = pickle.load(handle)  
    df.columns = map(str.lower, df.columns)


    # ----------REWORK AGAIN----------------
    print(f"[{dt.now()}] Counting unique de Bruijn kmers ...")
    df['node_edge'] = df['node1']+df['node2'].str.rstrip().str[-1]
    counts = {}
    for kmer_edge in df['node_edge']:
        counts[kmer_edge] = counts.get(kmer_edge, 0) + 1
    out_df = pd.DataFrame.from_dict(counts, orient='index', columns=['count'])
    out_df = out_df.reset_index().rename(columns={'index':'node_edge'})
    out_df['node1'] = out_df.node_edge.str.rstrip().str[:-1]
    out_df['node2'] = out_df.node_edge.str.rstrip().str[1:]
    print(out_df)
    cols = out_df.columns.tolist()
    col_order = ['node1', 'node2', 'count']
    if return_ids:
        print(f"[{dt.now()}] Summarizing IDs per node ...")
        id_df = df.groupby('node_edge').agg(lambda x: ', '.join(set(x)))
        print(id_df)
        id_df.reset_index(inplace=True)
        print(id_df.head())
        ids = id_df.set_index('node_edge')['proteinid'].to_dict()
        print({k: ids[k] for k in list(ids)[:10]})
        out_df['proteinid'] = [map_kmer(i, ids) for i in out_df['node_edge']]
        col_order = ['node1', 'node2', 'count', 'proteinid']
    out_df = out_df[col_order]
    print(out_df)

    # ----------rework----------------
    # print(f"[{dt.now()}] Counting unique de Bruijn kmers ...")
    # df['tuple'] = df[['node1','node2']].apply(tuple, axis=1)
    # kmer_pairs = df.tuple.to_list()
    # count_dct = {n:kmer_pairs.count(n) for n in unique_everseen(kmer_pairs)}
    
    # print(f"[{dt.now()}] Converting de Bruijn kmer counts to table ... ")
    # out_df = df.drop('proteinid',axis=1).drop_duplicates()
    # out_df['count'] = [map_kmer(kmer, count_dct) for kmer in out_df.tuple]
    # if return_ids:
    #     id_df = df.groupby('tuple').agg(lambda x: ','.join(set(x))).drop(columns=['node1','node2'])
    #     id_dct = id_df.to_dict()['proteinid']
    #     out_df['ids'] = [map_kmer(kmer, id_dct) for kmer in out_df.tuple]
    # out_df = out_df.drop('tuple', axis=1)
    
    # --------old code below here-----
    # # unique the kmer edges & aggregate the protein IDs in column 3             
    # uniques_df = kmer_df.groupby(['Node1','Node2'],sort=False).agg(lambda x: set(x)).reset_index()
    
    # # change the protein ID aggregate lists to a less annoying format
    # # and also create a protein count column
    # count_list = []
    # prot_list = []

    # for entry in uniques_df['ProteinID']:
    #     count_list.append(len(entry))
    #     prot_list.append(list(entry))

    # uniques_df['ProteinID'] = prot_list
    # uniques_df.insert(3, 'ProteinCount', count_list)      
    # ----------------------------------
    
    # write the output to a new file
    print(f"[{dt.now()}] Writing results to {output_file} ...")
    out_df = out_df.sort_values('count', ascending=False)
    print(out_df)
    out_df.to_csv(output_file,index=False)
    if pickle_file:
        pkl_out = output_file.replace(".csv", ".pkl")
        print(f"[{dt.now()}] Writing results to {pkl_out} ...")
        out_df.to_pickle(pkl_out)

def main():

    unique_kmer_file = unique(args.input_file, args.outfile, args.pickle, args.return_ids)

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description="Collapses an annotated deBruijn kmer file into uniques while aggregating IDs")
    parser.add_argument("--input_file", action="store", required=True,
                                        help="Filename for kmer edge list (.csv)")
    parser.add_argument("--sep", action="store", required=False, default=',',
                                        help="Column separator for input file, default=,")
    parser.add_argument("--outfile", action="store", required=False,
                                        help="(Optional) Filename for unique weighted kmer edge list (.csv)")
    parser.add_argument("--pickle", action="store_true", default=False, help="(Optional) Serialize results in a pickle file.")

    parser.add_argument("--return_ids", action="store_true", default=False, help="(Optional) Return IDs associated with each unique edge; not recommended for huge data sets.")
    
    args = parser.parse_args()
    main()
