from Bio import AlignIO
from Bio.Align import AlignInfo
import argparse
import re

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Calculates a consensus sequence for a multiple sequence alignment in FASTA format")
    parser.add_argument("--input_msa", action="store", required=True,
                                        help="Filename for multiple sequence alignment")
    parser.add_argument("--cutoff", action="store", required=True,
                                        help="Between 0 and 1; cutoff for % sequences need to be identifical \
                                        at each position to return a consensus")
    inputs = parser.parse_args()
    input_msa = inputs.input_msa
    cutoff = float(inputs.cutoff)

alignment = AlignIO.read(input_msa, "fasta")
summary = AlignInfo.SummaryInfo(alignment)
consensus = summary.dumb_consensus(cutoff)
consensus_gap = summary.gap_consensus(cutoff)

writename = re.match("(.*)(?=\.)",input_msa).group()
cutoff_perc = int(cutoff*100)
with open("{}_consensus_seqs_{}perc.fasta".format(writename,cutoff_perc),"w") as f:
    f.write(">{}_naive_consensus;length={};cutoff={}%\n".format(writename,len(consensus),cutoff_perc))
    f.write(str(consensus+"\n"))
    f.write(">{}_gapped_consensus;length={};cutoff={}%\n".format(writename,len(consensus_gap),cutoff_perc))
    f.write(str(consensus_gap+"\n"))

print(">{}_naive_consensus;length={};cutoff={}".format(writename, len(consensus),cutoff))
print(consensus)
print(">{}_gapped_consensus;length={};cutoff={}".format(writename,len(consensus_gap),cutoff))
print(consensus_gap)

