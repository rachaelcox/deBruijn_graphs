import csv
import pandas as pd
import argparse
from Bio import SeqIO
from itertools import groupby

 # generate unique edge list with aggregate protein IDs

def unique(input_file):
    
    output_file = input_file.replace(".csv","")
    n_df = pd.read_csv(input_file)
    n_df.columns = ['Node1','Node2','ProteinID']    

    kmer_df = n_df[['Node1','Node2','ProteinID']]   
    #print(kmer_df)               
    uniques_df = kmer_df.groupby(['Node1','Node2'],sort=False).agg(lambda x: set(x)).reset_index()

    count_list = []
    prot_list = []
    for entry in uniques_df['ProteinID']:
        count_list.append(len(entry))
        prot_list.append(list(entry))

    uniques_df['ProteinID'] = prot_list
    uniques_df.insert(3, 'ProteinCount', count_list)      

    uniques_df.to_csv("{}_uniques.csv".format(output_file),index=False)
        
    print(uniques_df)

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description="Collapses an annotated deBruijn kmer file into uniques while aggregating IDs")
    parser.add_argument("--input_file", action="store", required=True,
                                        help="Filename for kmer edge list (.csv)")
    parser.add_argument("--sep", action="store", required=False, default=',',
                                        help="Column separator for input file, default=,")
    
    inputs = parser.parse_args()
    unique_kmer_file = unique(inputs.input_file)
