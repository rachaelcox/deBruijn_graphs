# Quick start
```bash
# clone the repository
git clone https://github.com/rachaelcox/deBruijn_graphs.git
```

```bash
# specify desired k-mer size, input fasta, and outfile name
python3 deBruijn_graphs/scripts/generate_deBruijn_kmers.py \
--input_fasta petase_proteins.fasta \
--outfile petase_10mers.csv \
--kmer_size 10
```

```bash
# compute weighted graph 
python3 deBruijn_graphs/scripts/collapse_kmer_uniques.py \
--input_file petase_10mers.csv \
--outfile petase_10mers_counts.csv
```

```bash
# calculate consensus sequence from the weighted graph
python3 deBruijn_graphs/scripts/calc_consensus_seq.py --input_file petase_10mers_counts.csv
```

# Proof-of-concept
https://sites.google.com/utexas.edu/debruijnassemblyofproteins/introduction
