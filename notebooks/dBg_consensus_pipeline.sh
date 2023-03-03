filename=$1
kmersize=$2


#python3 /project/rmcox/deBruijn_graphs/scripts/generate_deBruijn_kmers.py --input_fasta fastas/${x} --kmer_size $y --output_file dbg_output/${x%.fa*}_dBg_${y}.csv

python3 /project/rmcox/deBruijn_graphs/scripts/generate_deBruijn_kmers.py --input_fasta ${filename}.fasta --kmer_size $kmersize

python3 /project/rmcox/deBruijn_graphs/scripts/collapse_kmer_uniques.py --input_file ${filename}_dBg_7mers.csv 

python3 /project/rmcox/deBruijn_graphs/scripts/calc_consensus_seq.py --input_file ${filename}_dBg_7mers_unique.csv

python3 /project/rmcox/deBruijn_graphs/scripts/generate_deBruijn_kmers.py --input_fasta ${filename}_dBg_7mers_unique.consensus.fasta --kmer_size ${kmersize}
