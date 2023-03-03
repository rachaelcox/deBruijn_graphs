import csv
import pandas as pd
import re
import argparse
from Bio import SeqIO
from itertools import groupby

# in a de Bruijn graph, edges=kmers and nodes=(k-1)mers
# specify FASTA file, desired kmer and number of proteins
def edgelist(input_fasta,kmer_size,num_prot):
    # define new file name, create empty list(s) and data frame(s)
    writefile = input_fasta.replace(".fasta",'')
    nodes = []
    n_df = pd.DataFrame()
    
    # extract lists of nodes and edges from FASTA file, then write edge list to csv
    with open("{}_nodes_{}mer_{}proteins.csv".format(writefile,kmer_size,num_prot),"w") as f:
        for index,record in enumerate((SeqIO.parse(open(input_fasta,"r"), "fasta"))):
            index += 1
            seq_str = record.seq
            ProteinID = record.id
            f.write("%s,%s,%s\n"%("Node1","Node2","ProteinID"))
            for i in range(len(seq_str)):
                n1 = seq_str[i:i+(kmer_size-1)]
                n2 = seq_str[i+1:i+kmer_size]
                if len(n1) and len(n2) == kmer_size-1:
                    nodes.append([str(n1),str(n2),ProteinID])
                    f.write("%s,%s,%s\n"%(n1,n2,ProteinID))
            if index == num_prot:
                break
    
    # generate unique edge list with aggregate protein IDs
    #n_df = pd.DataFrame(nodes)
    #n_df.columns = ['Node1','Node2','ProteinID']
    #uniques = n_df.groupby(['Node1','Node2'],sort=False).agg(lambda x: set(x)).reset_index()
    #uniques.to_csv("{}_nodes_{}mer_{}prots_uniques.csv".format(writefile,kmer_size,num_prot))
    #print(uniques)

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description="Takes a FASTA file and generates de Bruijn kmers")
    parser.add_argument("--input_fasta", action="store", required=True,
                                        help="Filename of .fasta input")
    parser.add_argument("--kmer_size", action="store", type=int, required=True,
                                        help="Size of kmer edges")
    parser.add_argument("--num_prot", action="store", required=False, default='all',
                                        help="Number of proteins to process")
    parser.add_argument("--sep", action="store", required=False, default=',',
                                        help="Column separator for input file, default=,")
    
    inputs = parser.parse_args()
    print(inputs)
    kmer_file = edgelist(inputs.input_fasta, inputs.kmer_size, inputs.num_prot)
