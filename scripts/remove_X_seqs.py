import csv
import argparse
from Bio import SeqIO

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Remove sequences containing X")
    parser.add_argument("--input_fasta", action="store", required=True,
                                        help="Filename for FASTA with sequences containing X")
    inputs = parser.parse_args()
    input_fasta = inputs.input_fasta

writefile = input_fasta.replace(".fasta",'')

with open("{}_noX.fasta".format(writefile), "w") as f:
    for record in SeqIO.parse(open(input_fasta,"r"), "fasta"):
        if 'X' not in record.seq:
            print(record.id)
            print(record.seq)
            f.write('>'+record.id+'\n')
            f.write(str(record.seq)+'\n')

