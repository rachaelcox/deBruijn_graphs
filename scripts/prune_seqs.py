"""
Script for pruning NGS reads to retain those with exact or close matches to the 5' and 3' flanking regions of a random sequence library.
@author: Rachael Cox <rachaelcox@utexas.edu>
"""

__author__ = "Rachael Cox (rachaelcox@utexas.edu)"
__version__ = "0.1.0"
__license__ = "MIT"

import argparse
from Bio import SeqIO
from Bio import Align
import re
import statistics
import time
from datetime import datetime as dt

""" Functions """

def mismatch_count(s1, s2):
    aligner = Align.PairwiseAligner(mode='global', match_score=2, mismatch_score=-1)
    alignment = aligner.align(s1, s2)
    score = alignment.score
    num_mm = sum(1 for a, b in zip(alignment[0][0], alignment[0][1]) if a != b)
    return num_mm

def print_stats(records):
    sizes = [len(rec) for rec in records]
    print("Mean read length:", statistics.mean(sizes))
    print("Median:", statistics.median(sizes))
    print("Mode:", statistics.mode(sizes))
    print("Max:", max(sizes))
    print("Min:", min(sizes))

""" Main """
def main(args):

    t0 = time.time()
    
    # 5' and 3' library flanking seqs
    f1 = 'TATTGCGATAGCTGAGAGAGAAGACGCGAGGG'
    f2 = 'GCGAAAACAAAAAACAAAAATAAGAATCCAAGCAGCAGCAACA'

    # read in NGS seqs
    print(f"[{dt.now()}] Reading in {args.fasta_file} ...")
    records = list(SeqIO.parse(args.fasta_file, "fasta"))

    # prune NGS seqs
    n = args.num_mismatch
    if args.exact:
        print(f"[{dt.now()}] Getting exact matches ...")
        matches = [rec for rec in records if rec.seq.startswith(f1) and rec.seq.endswith(f2)]
    else:
        print(f"[{dt.now()}] Getting reads where max # mismatches = {n} ...")
        matches = [rec for rec in records if mismatch_count(f1, rec.seq[0:len(f1)]) <= n and mismatch_count(f2, rec.seq[-len(f2):]) <= n]
        
    # print stats
    print(f"[{dt.now()}] ---------------------------------------------------------")
    print("Total reads in:", len(records))
    print_stats(records)
    print(f"[{dt.now()}] ---------------------------------------------------------")
    print("Total reads out:", len(matches))
    print_stats(matches)
    print(f"[{dt.now()}] ---------------------------------------------------------")
    print("% of reads retained:", round((len(matches)/len(records))*100, 2))
    print(f"[{dt.now()}] ---------------------------------------------------------")
    
    # write results
    print(f"[{dt.now()}] Writing results to {args.outfile} ...")
    SeqIO.write(matches, args.outfile, "fasta")

    print(f"[{dt.now()}] ---------------------------------------------------------")
    print(f"[{dt.now()}] Total run time: {round((time.time()-t0)/60, 2)} minutes.")
    print(f"[{dt.now()}] ---------------------------------------------------------")

if __name__ == "__main__":
    """ This is executed when run from the command line """
    parser = argparse.ArgumentParser()

    # required positional argument
    parser.add_argument("-f", "--fasta_file", action="store", dest="fasta_file", help="(Required) Fasta file to parse")

    # optional argument 
    parser.add_argument("-e", "--exact", action="store_true", default=False)

    # optional argument 
    parser.add_argument("-n", "--num_mismatch", action="store", dest="num_mismatch", type=int, default=5)

    # optional argument 
    parser.add_argument("-o", "--outfile", action="store", dest="outfile")

    # optional verbosity counter (eg. -v, -vv, -vvv, etc.)
    parser.add_argument(
        "-v",
        "--verbose",
        action="count",
        default=0,
        help="Verbosity (-v, -vv, etc)")

    # specify output of "--version"
    parser.add_argument(
        "--version",
        action="version",
        version="%(prog)s (version {version})".format(version=__version__))

    args = parser.parse_args()
    main(args)