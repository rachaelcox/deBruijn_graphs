import csv
import pandas as pd
import argparse
import pickle
from Bio import SeqIO
from itertools import groupby
from datetime import datetime as dt

# generate unique edge list with aggregate protein IDs
def unique(input_file, output_file, pickle_file, return_ids):
    
    # read input to dataframe & create new file name
    if not output_file:
        output_file = input_file.replace(".csv","")
        output_file = "{}_unique.csv".format(output_file)

    print(f"[{dt.now()}] Reading in de Bruijn graph from {input_file} ...")
    try:
        n_df = pd.read_csv(input_file)
    except:
        with open(input_file, 'rb') as handle:
            n_df = pickle.load(handle)
        
    n_df.columns = ['Node1','Node2','ProteinID']    
    kmer_df = n_df[['Node1','Node2','ProteinID']]   
    
    # ----------rework----------------
    print(f"[{dt.now()}] Counting unique de Bruijn kmers ...")
    edge_n_dict = dict()
    edge_id_dict = dict()
    for i in range(len(kmer_df)):
        n1 = kmer_df['Node1'][i]
        n2 = kmer_df['Node2'][i]
        pid = kmer_df['ProteinID'][i]
        edge = frozenset({n1,n2})
        if edge not in edge_n_dict.keys():
            edge_n_dict[edge] = 1
            edge_id_dict[edge] = set()
            edge_id_dict[edge].add(pid)
        elif set(pid) not in edge_id_dict[edge]:
            edge_n_dict[edge] += 1
            edge_id_dict[edge].add(pid)

    print(f"[{dt.now()}] Converting de Bruijn kmer counts to table ... ")
    edge_df = pd.DataFrame.from_dict(edge_n_dict, orient='index').reset_index()
    edge_df = edge_df.rename(columns={'index':'edge', 0:'count'})

    if return_ids:
        edge_df['ids'] = edge_df.edge.map(edge_id_dict)

    print(f"[{dt.now()}] Formatting de Bruijn kmer tabe ... ")
    node_df = pd.DataFrame(edge_df['edge'].values.tolist())
    node_df = node_df.rename(columns = lambda x: 'node{}'.format(x+1))
    edge_df.drop(['edge'], axis=True, inplace=True)
    uniques_df = pd.concat([node_df, edge_df], axis = 1)
    uniques_df = uniques_df.dropna()
    uniques_df.reset_index(inplace=True, drop=True)

    print(f"[{dt.now()}] Redirecting weighted de Bruijn graph with {len(uniques_df)} edges ... ")
    redirected = []
    for i in range(len(uniques_df)):
        node1 = uniques_df.iloc[i, 0]
        node2 = uniques_df.iloc[i, 1]
        if node1[:-1] == node2[1:]:
            node1 = uniques_df.iloc[i, 1]
            node2 = uniques_df.iloc[i, 0]
        redirected.append([node1, node2])

    redirected_df = pd.DataFrame(redirected, columns=['node1','node2'])
    redirected_df.reset_index(inplace=True, drop=True)
    redirected_df['count'] = uniques_df['count']
    if return_ids:
        redirected_df['ids'] = uniques_df['ids']
    
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
    redirected_df = redirected_df.sort_values('count', ascending=False)
    print(redirected_df)
    redirected_df.to_csv(output_file,index=False)
    if pickle_file:
        pkl_out = output_file.replace(".csv", ".pkl")
        print(f"[{dt.now()}] Writing results to {pkl_out} ...")
        redirected_df.to_pickle(pkl_out)
        
    #print(uniques_df)

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
    parser.add_argument("--pickle", action="store_true", default=None, help="(Optional) Serialize results in a pickle file.")

    parser.add_argument("--return_ids", action="store_true", default=None, help="(Optional) Return IDs associated with each unique edge; not recommended for huge data sets.")
    
    args = parser.parse_args()
    main()
