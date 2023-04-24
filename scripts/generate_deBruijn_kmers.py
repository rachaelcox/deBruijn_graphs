import csv
import re
import argparse
from Bio import SeqIO

# in a de Bruijn graph, edges=kmers and nodes=(k-1)mers
# specify FASTA file, desired kmer and number of proteins
def edgelist(input_fasta, kmer_size):
    
    # define new file name, create empty list(s) and data frame(s)
    writefile = input_fasta.replace(".fasta",'')
    nodes = []
    
    # extract lists of nodes and edges from FASTA file, then write edge list to csv
    with open("{}_dBg_{}mers.csv".format(writefile,kmer_size),"w") as f:
        
        f.write("%s,%s,%s\n"%("Node1","Node2","ProteinID"))
        for record in SeqIO.parse(open(input_fasta,"r"), "fasta"):
            
            seq_str = record.seq

            pattern = r'(tr|sp)\|(.*)\|\S*'
            match = re.search(pattern, str(record.id))
            if match:
                ProteinID = match.group(2)
            else:
                ProteinID = record.id

            for i in range(len(seq_str)):
                n1 = seq_str[i:i+(kmer_size-1)]
                n2 = seq_str[i+1:i+kmer_size]
                
                if len(n1) and len(n2) == kmer_size - 1:
                    f.write("%s,%s,%s\n"%(n1,n2,ProteinID))
                    print("{} edge #{} = {}-{}".format(ProteinID,i,n1,n2))

def main():

    kmer_file = edgelist(args.input_fasta, args.kmer_size)


if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description="Takes a FASTA file and generates de Bruijn kmers")
    parser.add_argument("--input_fasta", action="store", required=True,
                                        help="Filename of .fasta input")
    parser.add_argument("--kmer_size", action="store", type=int, required=True,
                                        help="Size of kmer edges") 
    args = parser.parse_args()
    main()