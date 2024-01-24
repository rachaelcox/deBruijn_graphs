import csv
import pandas as pd
import argparse
import re
import os
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
np.random.seed(666)

""" recursive function for growing kmer sequence in N-terminal/5' direction """
def grow_seq_n(seed_kmer, input_df, nterm_kmers=None, length=None, cons_seq_n=None):
    
    # if a list has not been passed as argument create an empty one
    if(nterm_kmers == None):
        nterm_kmers = []
    # if an length has not been passed as argument, initialize one from 0
    if (length == None):
        length = 0  
    # choose next kmer in the N-terminal direction
    if input_df['Node2'].str.contains(seed_kmer).any()\
        and seed_kmer not in nterm_kmers:  # greedy path only, no cycles
        length += 1
        nterm_kmers.append(seed_kmer)
        poss_paths_nterm = input_df[input_df['Node2'] == seed_kmer]
        #print("possible future directions =",poss_paths_nterm)

        # find edge with highest counts 
        most_common_kmer_n = poss_paths_nterm['Weight'] == poss_paths_nterm['Weight'].max()
        # if more than one most common edge, randomly pick one
        next_n_edge_df = poss_paths_nterm[most_common_kmer_n].sample(1)
        next_n_kmer = str(next_n_edge_df.iloc[0]['Node1'])
        # recursively add kmers to the consensus sequence
        grow_seq_n(next_n_kmer, input_df, nterm_kmers, length)
        
    # after the kmer finding terminates
    else:
        #nterm_kmers.append(seed_kmer) #for locating where the first repeat is found

        print(f"Found beginning of sequence; total downstream k-mers traveled = {len(nterm_kmers)}")
        return nterm_kmers        

    # calculate the consensus sequence from the most common kmer list
    cons_seq_n = nterm_kmers[0]
    for kmer in nterm_kmers[1:]:
        cons_seq_n = kmer[0]+cons_seq_n       

    # exit function and return consensus n-terminal kmers
    return cons_seq_n

""" function for recursively growing kmer sequence in C-terminal/3' direction """

def grow_seq_c(seed_kmer, input_df, cterm_kmers=None, length=None):

    # if a list has not been passed as argument create an empty one
    if(cterm_kmers == None):
        cterm_kmers = []
    # if an length has not been passed as argument, initialize one from 0
    if (length == None):
        length = 0
    # choose next kmer in the N-terminal direction
    if input_df['Node1'].str.contains(seed_kmer).any()\
        and seed_kmer not in cterm_kmers:  # greedy path only, no cycles
        length += 1
        cterm_kmers.append(seed_kmer)
        poss_paths_cterm = input_df[input_df['Node1'] == seed_kmer]
        # find edge with highest counts
        most_common_kmer_c = poss_paths_cterm['Weight'] == poss_paths_cterm['Weight'].max()
        # if more than one most common edge, randomly pick one
        next_c_edge_df = poss_paths_cterm[most_common_kmer_c].sample(1)
        next_c_kmer = str(next_c_edge_df.iloc[0]['Node2'])
        # recursively add kmers to the consensus sequence
        grow_seq_c(next_c_kmer, input_df, cterm_kmers, length)

    # after the kmer finding terminate
    else:
        #cterm_kmers.append(seed_kmer) #for locating where the first repeat is found                                     
        print(f"Found end of sequence; total upstream k-mers traveled = {len(cterm_kmers)}")
        return cterm_kmers
    
    # calculate the consensus sequence from the most common kmer list
    cons_seq_c = cterm_kmers[0]
    for kmer in cterm_kmers[1:]:
        cons_seq_c = cons_seq_c + kmer[-1]
    
    # exit function and return consensus n-terminal kmers
    return cons_seq_c

""" functions for optional parameters """

def format_out(infile_name, outfile_name=None):
    # format outfile paths/names
    if outfile_name:
        path, basename = os.path.split(os.path.realpath(outfile_name))
        outpath = path+'/'
    else:
        path, basename = os.path.split(os.path.realpath(infile_name))
        outpath = path+'/'
        basename = basename.replace(".csv","")
    return(outpath, basename)
    
def seed_diversity(input_df):
    input_df['edge'] = [i[1:] for i in input_df.Node1]
    input_df['div'] = [len(set(i)) for i in input_df.edge]
    div_df = input_df.sort_values(['div','Weight'], ascending=False)
    return div_df

def grow_cons_seq(seed_kmer_1, seed_kmer_2, input_df):
    # grow sequence recursively in both directions
    nterm_seq = grow_seq_n(seed_kmer_1, input_df)
    print(f"Upstream sequence [length={len(nterm_seq)}]:\n{nterm_seq}")
    cterm_seq = grow_seq_c(seed_kmer_2, input_df)
    print(f"Downstream sequence [length={len(cterm_seq)}]:\n{cterm_seq}")
    # format output in the event connecting edges are not found
    try:
        consensus_seq = str(nterm_seq)+str(cterm_seq)
        if len(nterm_seq)==0 and len(cterm_seq)==0:
            print('WARNING: No connecting paths found.')
            consensus_seq = None
    except TypeError:
        if len(cterm_seq) < len(seed_kmer_1) and len(nterm_seq) > len(seed_kmer_1):
            consensus_seq = nterm_seq+seed_kmer_2[-1]
        elif len(nterm_seq) < len(seed_kmer_1) and len(cterm_seq) > len(seed_kmer_1):
            consensus_seq = seed_kmer_1[0]+cterm_seq
    if consensus_seq:
        consensus_seq = consensus_seq.replace('[]','')
        print(f"Consensus [length={len(consensus_seq)}]: {consensus_seq}")
    return consensus_seq

""" wrappers """

## WORK IN PROGRESS
def calc_n_seqs(input_df, n, seed_div=False):
    
    # get n starting seeds with highest weights
    if seed_div:
        print('Prioritizing lexically diverse k-mers ...')
        div_df = seed_diversity(input_df)
        print(div_df)
        ranked_edges = div_df.head(n).reset_index(drop=True)
    else:
        ranked_edges = input_df.nlargest(n, 'Weight').reset_index(drop=True)

    # loop through each of the top kmers as a starting point
    rows = []
    records = []
    print()
    for idx, row in ranked_edges.iterrows():
        rank = idx+1
        seed_kmer_1, seed_kmer_2, score = row[['Node1','Node2','Weight']]
        print(f'{seed_kmer_1} {seed_kmer_2}\trank:{rank}\tscore:{score}')
        seed_kmer_seq = seed_kmer_1[:-1]+seed_kmer_2[-1]
        header = f'[seed={seed_kmer_seq} | rank={rank} | score={score}]'
        cons_seq = grow_cons_seq(seed_kmer_1, seed_kmer_2, input_df)
        record = SeqRecord(Seq(str(cons_seq)), id=header)
        records.append(record)
        rows.append([seed_kmer_seq, rank, score, cons_seq])
        print()
    df_out = pd.DataFrame(rows, columns=['seed_kmer','rank','score','consensus_seq'])
    print(f'Consensus seqs from top {len(ranked_edges)} kmers:')
    print(df_out)
    return records, df_out
            
def calc_consensus_seq(input_df, seed_div=False):

    # find the most common edge for starting seeds
    most_common_edge = input_df['Weight'] == input_df['Weight'].max()

    # if more than one most common edge, randomly pick one
    seed_df = input_df[most_common_edge].sample(1)

    # assign kmers as seeds & initialize the consensus list
    if seed_div:
        div_df = seed_diversity(input_df)
        seed_kmer_1 = str(div_df.iloc[0]['Node1'])
        seed_kmer_2 = str(div_df.iloc[0]['Node2'])
        print(seed_kmer_1, '\t', seed_kmer_2)
        consensus_seq = grow_cons_seq(seed_kmer_1, seed_kmer_2, input_df)
        # --- try again w/o diversity if bad seed
        if len(consensus_seq) == len(seed_kmer_1+1):
            print('Trying again w/o diversity seed ...')
            seed_kmer_1 = str(seed_df.iloc[0]['Node1'])
            seed_kmer_2 = str(seed_df.iloc[0]['Node2'])
            consensus_seq = grow_cons_seq(seed_kmer_1, seed_kmer_2, input_df)
    else:
        seed_kmer_1 = str(seed_df.iloc[0]['Node1'])
        seed_kmer_2 = str(seed_df.iloc[0]['Node2'])
        consensus_seq = grow_cons_seq(seed_kmer_1, seed_kmer_2, input_df)

    # print output
    print("\nthe consensus sequence (length {} aa) is: \n{}".\
          format(len(consensus_seq), consensus_seq))
    return consensus_seq, seed_kmer_1, seed_kmer_2

""" main """

def main():

    # read in and format data
    dBg_unique = args.input_file
    input_df = pd.read_csv(dBg_unique)
    cols=['Node1','Node2','Weight']
    input_df.columns = cols
    input_df.sort_values('Weight', ascending=False, inplace=True)

    # format outfile name
    outpath, basename = format_out(args.input_file, args.outfile_prefix)

    if not args.calc_n:
        # for one consensus seq
        cons_seq, sk1, sk2 = calc_consensus_seq(input_df)
        # format new file name & write consensus sequence to output
        header = f'consensus sequence [{basename}] [seed={sk1[:-1]+sk2[-1]} | rank=1]'
        fasta = SeqRecord(Seq(str(cons_seq)), id=header)
    else:
        # for many consensus seqs
        fasta, consensus_df = calc_n_seqs(input_df, args.calc_n, args.seed_diversity)
        consensus_df.to_csv(f'{outpath+basename}.csv', index=False)
        for entry in fasta:
            entry.id = f'consensus sequence [{basename}] {entry.id}'
    
    with open(f'{outpath+basename}.fasta', "w") as f:
        SeqIO.write(fasta, f, "fasta")
    

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Calculates a consensus sequence for a deBruijn kmer file with unique kmers")
    parser.add_argument("--input_file", action="store", required=True,
                                        help="Filename for unique kmer edge list (.csv)")
    parser.add_argument("--seed_diversity", action="store_true", required=False, default=False, help="(Optional) Use the (1) most lexically diverse and (2) highest weight k-mer as a seed.")
    parser.add_argument("--calc_n", action="store", required=False, default=None, type=int,
                                        help="(Optional) Number of consensus sequences to generate.")
    parser.add_argument("-o", "--outfile_prefix", action="store", default=None, help="(Optional) Specify the outfile path/name. Default='{input_name}.consensus.fasta'")
    args = parser.parse_args()
    main()
    

