import csv
import argparse
from Bio import SeqIO

def extract_prot_lengths(input_fasta):

    # set up renamed output file
    database_file = input_fasta
    database = database_file.replace('.fasta','')

    # open file to write and specify headers
    with open("{}_lengths.txt".format(database),"w") as f:
        f.write("{}\t{}\n".format("ProteinID","ProteinLength"))

        # loop through the fasta file and extract protein ID(s) and lengths
        for record in SeqIO.parse(open(database_file,"r"), "fasta"):
            prot_id = record.id
            prot_len = len(record.seq)

            # write the ID(s) and lengths to the open csv
            f.write("%s\t%i\n"%(prot_id, prot_len))

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Takes a FASTA file and generates tab-delimited                                                        file with protein IDs and lengths")
    parser.add_argument("--input_fasta", action="store", required=True,
                                        help="Filename of .fasta input")
    inputs = parser.parse_args()
    print(inputs)
    outputs = extract_prot_lengths(inputs.input_fasta)
