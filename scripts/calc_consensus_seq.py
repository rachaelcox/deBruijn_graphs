import csv
import pandas as pd
import argparse
import re
import numpy as np
np.random.seed(13)

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Calculates a consensus sequence for a deBruijn kmer file with unique kmers")
    parser.add_argument("--input_file", action="store", required=True,
                                        help="Filename for unique kmer edge list (.csv)")
    inputs = parser.parse_args()
    dBg_unique = inputs.input_file

########################################################################
# recursive function for growing kmer sequence in N-terminal direction
########################################################################
def grow_seq_n(seed_kmer, input_df, nterm_kmers=None, length=None, cons_seq_n=None):
    
    # if a list has not been passed as argument create an empty one
    if(nterm_kmers == None):
        nterm_kmers = []
    
    # if an length has not been passed as argument, initialize one from 0
    if (length == None):
        length = 0
            
    # choose next kmer in the N-terminal direction
    if input_df['Node2'].str.contains(seed_kmer).any()\
        and seed_kmer not in nterm_kmers:
        
        length += 1
        #print("calculating kmer#{}...".format(length))
        nterm_kmers.append(seed_kmer)
        
        poss_paths_nterm = input_df[input_df['Node2'] == seed_kmer]
        ##print("possible future directions =",poss_paths_nterm)

        # find edge with highest counts 
        most_common_kmer_n = poss_paths_nterm['ProteinCount'] == poss_paths_nterm['ProteinCount'].max()
        
        # if more than one most common edge, randomly pick one
        next_n_edge_df = poss_paths_nterm[most_common_kmer_n].sample(1)
        next_n_kmer = str(next_n_edge_df.iloc[0]['Node1'])

        #print("\nnext consensus kmer =",next_n_kmer)
        
        # recursively add kmers to the consensus sequence
        grow_seq_n(next_n_kmer, input_df, nterm_kmers, length)
        
    # after the kmer finding terminates
    else:
        #nterm_kmers.append(seed_kmer) #for locating where the first repeat is found

        #print("found end of sequence")
        #print("\nlength = {}; kmers = {}".format(length,nterm_kmers))
        return nterm_kmers        

    # calculate the consensus sequence from the most common kmer list
    cons_seq_n = nterm_kmers[0]
    for kmer in nterm_kmers[1:]:
        cons_seq_n = kmer[0]+cons_seq_n       

    # exit function and return consensus n-terminal kmers
    return cons_seq_n

# function for recursively growing kmer sequence in N-terminal direction

def grow_seq_c(seed_kmer, input_df, cterm_kmers=None, length=None):

    # if a list has not been passed as argument create an empty one
    if(cterm_kmers == None):
        cterm_kmers = []

    # if an length has not been passed as argument, initialize one from 0
    if (length == None):
        length = 0

    # choose next kmer in the N-terminal direction
    if input_df['Node1'].str.contains(seed_kmer).any()\
        and seed_kmer not in cterm_kmers:

        length += 1
        #print("calculating kmer#{}...".format(length))
        cterm_kmers.append(seed_kmer)

        poss_paths_cterm = input_df[input_df['Node1'] == seed_kmer]
        ##print("possible future directions =",poss_paths_cterm)

        # find edge with highest counts
        most_common_kmer_c = poss_paths_cterm['ProteinCount'] == poss_paths_cterm['ProteinCount'].max()

        # if more than one most common edge, randomly pick one
        next_c_edge_df = poss_paths_cterm[most_common_kmer_c].sample(1)
        next_c_kmer = str(next_c_edge_df.iloc[0]['Node2'])

        #print("\nnext consensus kmer =", next_c_kmer)

        cons_seq_c = cterm_kmers[0]
        for kmer in cterm_kmers[1:]:
            cons_seq_c = cons_seq_c + kmer[-1] 
        #print(cons_seq_c)

        # recursively add kmers to the consensus sequence
        grow_seq_c(next_c_kmer, input_df, cterm_kmers, length)

       

    # after the kmer finding terminate
    else:
        #nterm_kmers.append(seed_kmer) #for locating where the first repeat is found                                                                         
        #print("found end of sequence")
        #print("\nlength = {}; kmers = {}".format(length,cterm_kmers))
        return cterm_kmers
    
    # calculate the consensus sequence from the most common kmer list
    cons_seq_c = cterm_kmers[0]
    for kmer in cterm_kmers[1:]:
        cons_seq_c = cons_seq_c + kmer[-1]
    
    # exit function and return consensus n-terminal kmers
    return cons_seq_c

########################################################################
# primary script
########################################################################

# read in file
input_df = pd.read_csv(dBg_unique)

# find the most common edge for starting seeds
most_common_edge = input_df['ProteinCount'] == input_df['ProteinCount'].max()

# if more than one most common edge, randomly pick one
seed_df = input_df[most_common_edge].sample(1)

# assign kmers as seeds & initialize the consensus list
seed_kmer_1 = str(seed_df.iloc[0]['Node1'])
seed_kmer_2 = str(seed_df.iloc[0]['Node2'])

# grow sequence recursively
nterm_seq = grow_seq_n(seed_kmer_1, input_df)
#print("nterm seq = {}\nlength = {}".format(nterm_seq,len(nterm_seq)))
cterm_seq = grow_seq_c(seed_kmer_2, input_df)
#print("cterm seq = {}\nlength = {}".format(cterm_seq,len(cterm_seq)))

consensus_seq = nterm_seq+cterm_seq
#print("\nthe consensus sequence (length {} aa) is: \n{}".\
      #format(len(consensus_seq), consensus_seq))

# format new file name & write consensus sequence to output
writefile = dBg_unique.replace(".csv","")
writename = re.match("(.*)(?=_)",dBg_unique).group()

with open("{}.consensus.fasta".format(writefile),"w") as f:
    f.write(">{} consensus sequence; length = {}\n".format(writename, len(consensus_seq)))
    f.write(consensus_seq+"\n")

