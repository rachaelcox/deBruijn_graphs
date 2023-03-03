import csv
import pandas as pd
import argparse
from Bio import SeqIO
from itertools import groupby

 # generate unique edge list with aggregate protein IDs

def alph_unique(input_file):
    
    output_file = input_file.replace(".csv","")
    n_df = pd.read_csv(input_file)
    n_df.columns = ['Node1','Node2','ProteinID']   
   
    columns2alphabetize = ['Node1','Node2']
    print(n_df[columns2alphabetize].head())
    intermediate_df = n_df[columns2alphabetize].apply(sorted,axis=1)
    n_df[columns2alphabetize] = intermediate_df
  #     intermediate_df = n_df[n_df.columns[columns2alphabetize]].apply(sorted,axis=1)
   # n_df[n_df.columns[columns2alphabetize]] = intermediate_df
    
    print(n_df)              
    uniques_df = n_df.groupby(['Node1','Node2'],sort=False).agg(lambda x: set(x)).reset_index()
    
    count_list = []
    prot_list = []
    for entry in uniques_df['ProteinID']:
        count_list.append(len(entry))
        prot_list.append(list(entry))

    uniques_df['ProteinID'] = prot_list
    uniques_df.insert(3, 'ProteinCount', count_list)

    #uniques["ProteinID"] = pd.to_list(uniques["ProteinID"] # going to be some type of apply )

    uniques_df.to_csv("{}_alph_uniques.csv".format(output_file),index=False)
        
    print(uniques_df)

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description="Takes a file with annonated de Bruijn kmers, alphabetizes and collapes into uniques")
    parser.add_argument("--input_file", action="store", required=True,
                                        help="Filename for kmer edge list (.csv)")
    parser.add_argument("--sep", action="store", required=False, default=',',
                                        help="Column separator for input file, default=,")
    
    inputs = parser.parse_args()
    unique_kmer_file = alph_unique(inputs.input_file)
