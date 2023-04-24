import csv
import pandas as pd
import argparse
from Bio import SeqIO
from itertools import groupby

# generate unique edge list with aggregate protein IDs
def unique(input_file):
    
    # read input to dataframe & create new file name
    output_file = input_file.replace(".csv","")
    n_df = pd.read_csv(input_file)
    n_df.columns = ['Node1','Node2','ProteinID']    
    kmer_df = n_df[['Node1','Node2','ProteinID']]   
    
    # unique the kmer edges & aggregate the protein IDs in column 3             
    uniques_df = kmer_df.groupby(['Node1','Node2'],sort=False).agg(lambda x: set(x)).reset_index()
    
    # change the protein ID aggregate lists to a less annoying format
    # and also create a protein count column
    count_list = []
    prot_list = []

    for entry in uniques_df['ProteinID']:
        count_list.append(len(entry))
        prot_list.append(list(entry))

    uniques_df['ProteinID'] = prot_list
    uniques_df.insert(3, 'ProteinCount', count_list)      

    # write the output to a new file
    uniques_df = uniques_df.sort_values('ProteinCount', ascending=False)
    uniques_df.to_csv("{}_unique.csv".format(output_file),index=False)
        
    #print(uniques_df)

def main():

    unique_kmer_file = unique(inputs.input_file)

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description="Collapses an annotated deBruijn kmer file into uniques while aggregating IDs")
    parser.add_argument("--input_file", action="store", required=True,
                                        help="Filename for kmer edge list (.csv)")
    parser.add_argument("--sep", action="store", required=False, default=',',
                                        help="Column separator for input file, default=,")
    
    inputs = parser.parse_args()
    main()
