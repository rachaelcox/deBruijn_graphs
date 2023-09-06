# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Oligo reconstruction with de bruijn graphs
# Rachael M. Cox 6/2023
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# project dir
rmcox@hfogcomp01:/stor/work/Marcotte/project/rmcox/deBruijn_graphs/oligo

# =======================================
# general notes
# =======================================

# 4 baits; 3 around ~40nt, 1 blind sequence
# prey library is ~90nt with N15 variable region flanked by 32nt 5' and 43nt 3' 
# N15 sequence interval: [33,42]

# old data dir
/stor/work/Ellington/scratch/Zac/PLA/data/new_data/*Lib12*
/stor/work/Ellington/scratch/Zac/PLA/data/new_data/*Probe12* # same as Lib12?
/stor/work/Ellington/scratch/Zac/PLA/data/new_data/*Probe34* # what are these?
/stor/work/Ellington/scratch/Zac/PLA/data/new_data/*Oligo3-Undetermined* # what is this?w

# organization
mkdir -p data/raw/{oligo_1,oligo_2,oligo_3,oligo_b}
mkdir -p data/qc/{oligo_1,oligo_2,oligo_3,oligo_b}
mkdir -p data/processed/{oligo_1,oligo_2,oligo_3,oligo_b}

ln -s /stor/work/Ellington/scratch/Zac/PLA/data/new_data/*Oligo1* /stor/work/Marcotte/project/rmcox/deBruijn_graphs/oligo_hyb/data/raw/oligo_1
ln -s /stor/work/Ellington/scratch/Zac/PLA/data/new_data/*Oligo2* /stor/work/Marcotte/project/rmcox/deBruijn_graphs/oligo_hyb/data/raw/oligo_2
ln -s /stor/work/Ellington/scratch/Zac/PLA/data/new_data/*Oligo3* /stor/work/Marcotte/project/rmcox/deBruijn_graphs/oligo_hyb/data/raw/oligo_3
ln -s /stor/work/Ellington/scratch/Zac/PLA/data/new_data/*OligoBlind* /stor/work/Marcotte/project/rmcox/deBruijn_graphs/oligo_hyb/data/raw/oligo_b

# exp set up: oligo|wash|rep
# naming scheme is not consistent T_T

# =======================================
# data QC
# =======================================

fastqc --noextract --nogroup -o ../qc/ */* 

for x in 1 2 3 b; do echo "multiqc data/qc/oligo_${x}/ --filename data/qc/multiqc/oligo_${x} --interactive"; done > cmds/multiqc.cmds
cat cmds/multiqc.cmds | parallel -j4

# oo i'm going to try fastp, seems like a good all in 1 tool
# template command
fastp -i R1.fastq.qz -I R2.fastq.gz -o R1.trimmed.fastq.gz -o R2.trimmed.fastq.gz -p --failed_out data/qc/fastp_failed/

for x in 1 2 3 b; do ls -C data/raw/oligo_${x}/*R1* > data/meta/oligo_${x}_R1.txt; done
for x in 1 2 3 b; do ls -C data/raw/oligo_${x}/*R2* > data/meta/oligo_${x}_R2.txt; done
for x in 1 2 3 b; do paste data/meta/oligo_${x}_R1.txt data/meta/oligo_${x}_R2.txt > data/meta/oligo_${x}_R1R2.txt; done

# testing loop template
for x in 1 2 3 b; do
	while IFS=$'\t' read -r -a y; do
    	echo "R1: ${y[0]} | R2: ${y[1]}"
    	echo "O1: ${y[0]##*/} | O2: ${y[1]##*/}"
    	b1="$(basename ${y[0]} .fastq.gz)"
    	b2="$(basename ${y[1]} .fastq.gz)"
    	#echo "B1: data/qc/oligo_${x}/${b1}.fpqc.fastq.gz | B2: data/qc/oligo_${x}/${b2}.fpqc.fastq.gz"
    	merg="${b1%_L00*}.merged.fpqc.fastq.gz"
    	echo "M: ${merg}"
	done < data/meta/oligo_${x}_R1R2.txt
done

# finalized loop
for x in 1 2 3 b; do
	while IFS=$'\t' read -r -a y; do
    	o1="$(basename ${y[0]} .fastq.gz)"
    	o2="$(basename ${y[1]} .fastq.gz)"
    	merg="${o1%_L00*}.merged.fpqc.fastq.gz"
    	echo "fastp -i ${y[0]} -I ${y[1]} -o data/qc/oligo_${x}/${o1}.fpqc.fastq.gz -O data/qc/oligo_${x}/${o2}.fpqc.fastq.gz --failed_out data/qc/fastp/failed/${o1%_L00*}.failed.txt --merge --merged_out data/processed/oligo_${x}/${merg} --length_limit 100 --json data/qc/fastp/oligo_${x}_${o1%_L00*}.json --html data/qc/fastp/oligo_${x}_${o1%_L00*}.html"
	done < data/meta/oligo_${x}_R1R2.txt
done > cmds/fastp.cmds
cat cmds/fastp.cmds | parallel -j24

# qc trimmed and merged reads
for x in 1 2 3 b; do echo "fastqc --noextract --nogroup --threads 16 -o data/qc/oligo_${x}_merged data/processed/oligo_${x}/*"; done > cmds/fastq_fastp.cmds
cat cmds/fastq_fastp.cmds | parallel -j4

for x in 1 2 3 b; do echo "multiqc data/qc/oligo_${x}_merged/* --filename data/qc/multiqc/oligo_${x}_merged --interactive --force"; done > cmds/multiqc_fastp.cmds
cat cmds/multiqc_fastp.cmds | parallel -j4

# reports:
# oligo1: file:///H:/project/rmcox/deBruijn_graphs/oligo_hyb/data/qc/multiqc/oligo_1_merged.html
# oligo2: file:///H:/project/rmcox/deBruijn_graphs/oligo_hyb/data/qc/multiqc/oligo_2_merged.html
# oligo3: file:///H:/project/rmcox/deBruijn_graphs/oligo_hyb/data/qc/multiqc/oligo_3_merged.html
# oligob: file:///H:/project/rmcox/deBruijn_graphs/oligo_hyb/data/qc/multiqc/oligo_b_merged.html

# =======================================
# data munging
# =======================================

# merge experiments
cat file1.gz file2.gz file3.gz > allfiles.gz

# rename dirs with better clarity
# my future self will thank me
for x in 1 2 3 b; do 
	mv data/processed/oligo_${x} data/processed/oligo_${x}_pe_merged
	mkdir data/processed/oligo_${x}_wash_merged
done

for x in 1 2 3 b; do
	ls data/processed/oligo_${x}_pe_merged/*Probe* > data/meta/oligo_${x}_bg.txt
	ls data/processed/oligo_${x}_pe_merged/*0w* > data/meta/oligo_${x}_0w.txt
	ls data/processed/oligo_${x}_pe_merged/*1w* > data/meta/oligo_${x}_1w.txt
	ls data/processed/oligo_${x}_pe_merged/*3w* > data/meta/oligo_${x}_3w.txt
done

for x in 1 2 3 b; do
	echo "{ xargs cat < data/meta/oligo_${x}_bg.txt ; } > data/processed/oligo_${x}_wash_merged/oligo_${x}_bg.merged.fastq.gz" 
	echo "{ xargs cat < data/meta/oligo_${x}_0w.txt ; } > data/processed/oligo_${x}_wash_merged/oligo_${x}_0w.merged.fastq.gz"
	echo "{ xargs cat < data/meta/oligo_${x}_1w.txt ; } > data/processed/oligo_${x}_wash_merged/oligo_${x}_1w.merged.fastq.gz"
	echo "{ xargs cat < data/meta/oligo_${x}_3w.txt ; } > data/processed/oligo_${x}_wash_merged/oligo_${x}_3w.merged.fastq.gz"
done > cmds/concat_merged_reads.sh
cat cmds/concat_merged_reads.sh | parallel -j12

# decompress
for x in 1 2 3 b; do gunzip -f data/processed/oligo_${x}_wash_merged/*gz; done

# convert to fasta format
for x in 1 2 3 b; do
	for y in bg 0w 1w 3w; do
	sed -n '1~4s/^@/>/p;2~4p' data/processed/oligo_${x}_wash_merged/oligo_${x}_${y}.merged.fastq > data/processed/oligo_${x}_wash_merged/oligo_${x}_${y}.merged.fasta
	done
done

# =======================================
# generate de Bruijn graphs
# =======================================

for k in {5..15}; do
	for x in 1 2 3 b; do
	for y in bg 0w 1w 3w; do
	echo "python3 ../scripts/generate_deBruijn_kmers.py --input_fasta data/processed/oligo_${x}_wash_merged/oligo_${x}_${y}.merged.fasta --kmer_size $k --outfile results/${k}mer/oligo_${x}_${y}_${k}mers.csv"
	done
done > cmds/generate_${k}mer_dbgs.sh
done

rm cmds/generate_5to15mer_dbgs.sh
touch cmds/generate_5to15mer_dbgs.sh
for k in {5..15}; do
	cat cmds/generate_${k}mer_dbgs.sh >> cmds/generate_5to15mer_dbgs.sh
done
cat cmds/generate_5to15mer_dbgs.sh | parallel -j16

for k in {5..15}; do
	for x in 1 2 3 b; do
	for y in bg 0w 1w 3w; do
	echo "python3 ../scripts/collapse_kmer_uniques.py --input_file results/${k}mer/oligo_${x}_${y}_${k}mers.csv --outfile results/${k}mer/oligo_${x}_${y}_${k}mers_unique.csv | tee -a logs/oligo_${x}_${y}_${k}mers_collapse_unique.log"
	done
done > cmds/generate_${k}mer_weights.sh
done 

rm cmds/generate_5to15mer_weights.sh
touch cmds/generate_5to15mer_weights.sh
for k in {5..15}; do
	cat cmds/generate_${k}mer_weights.sh >> cmds/generate_5to15mer_weights.sh
done
cat cmds/generate_5to15mer_weights.sh | parallel -j8

python3 ../scripts/collapse_kmer_uniques.py --input_file results/10mer/oligo_1_0w_10mers.csv --outfile results/10mer/oligo_1_0w_10mers_unique.csv
python3 ../scripts/collapse_kmer_uniques.py --input_file results/10mer/oligo_b_3w_10mers.pkl --outfile results/10mer/oligo_b_3w_10mers_unique.csv

python3 ../scripts/calc_consensus_seq.py --input_file results/15mer/oligo_b_3w_15mers_unique.csv

python3 ../scripts/calc_consensus_seq.py --input_file results/15mer/oligo_1_3w_15mers_unique.csv

# yeah this does not work with k=10, highly redundant sections of the N15 library are making it too difficult for the consensus algorithm to find a sequence

# ---------------------------------------
# make graphs from bait and prey sequences
# ---------------------------------------

for k in {10..15}; do python3 ../scripts/generate_deBruijn_kmers.py --input_fasta data/fastas/prey_comp.fasta --kmer_size ${k} --outfile data/bait_kmers/library_comp_${k}mers.csv; done
for k in {10..15}; do python3 ../scripts/generate_deBruijn_kmers.py --input_fasta data/fastas/prey.fasta --kmer_size ${k} --outfile data/bait_kmers/library_${k}mers.csv; done

for k in {11..15}; do
	for x in 1 2 3; do
	python3 ../scripts/generate_deBruijn_kmers.py --input_fasta data/fastas/oligo_${x}_comp.fasta --kmer_size ${k} --outfile data/bait_kmers/oligo_${x}_${k}mers.csv
	done
done


for x in 1 2 3; do python3 ../scripts/generate_deBruijn_kmers.py --input_fasta data/fastas/oligo_${x}.fasta --kmer_size ${k} --outfile data/bait_kmers/oligo_${x}_${k}mers.csv; done
for x in 1 2 3; do python3 ../scripts/collapse_kmer_uniques.py --input_file data/bait_kmers/oligo_${x}_${k}mers.csv --outfile data/bait_kmers/oligo_${x}_${k}mers_unique.csv; done

# ---------------------------------------
# generate consensus sequences from filtered kmers
# ---------------------------------------

for i in 1 2 3; do
	python3 ../scripts/calc_consensus_seq.py --input_file results/thresholded/oligo_${i}_10mer_graph_fdr1e-10.csv --seed_diversity
done

for i in 1 2 3; do
	python3 ../scripts/calc_consensus_seq.py --input_file results/thresholded/oligo_${i}_10mer_graph_fdr1e-10_3log2fc.csv --seed_diversity
done

for i in 1 2 3; do
	python3 ../scripts/calc_consensus_seq.py --input_file results/thresholded/oligo_${i}_10mer_graph_fdr1e-10_3w.csv
done

for i in 1 2 3; do
	cat results/thresholded/oligo_${i}_10mer_graph_fdr1e-10.consensus.fasta results/thresholded/oligo_${i}_10mer_graph_fdr1e-10_3log2fc.consensus.fasta data/fastas/oligo_${i}_comp.fasta data/fastas/prey.fasta > results/alignment_fastas/oligo_${i}_cons_bait.fasta
done

# --------------------------------------------------

python3 ../scripts/generate_deBruijn_kmers.py --input_fasta data/f.fasta --kmer_size 7 --outfile results/f_7mers.csv
python3 ../scripts/generate_deBruijn_kmers.py --input_fasta data/k.fasta --kmer_size 7 --outfile results/k_7mers.csv
python3 ../scripts/generate_deBruijn_kmers.py --input_fasta data/g.fasta --kmer_size 7 --outfile results/g_7mers.csv
python3 ../scripts/generate_deBruijn_kmers.py --input_fasta data/p.fasta --kmer_size 7 --outfile results/p_7mers.csv

python3 ../scripts/collapse_kmer_uniques.py --input_file results/f_7mers.csv --outfile results/f_7mers_uniq.csv --return_ids
python3 ../scripts/collapse_kmer_uniques.py --input_file results/k_7mers.csv --outfile results/k_7mers_uniq.csv --return_ids
python3 ../scripts/collapse_kmer_uniques.py --input_file results/g_7mers.csv --outfile results/g_7mers_uniq.csv --return_ids
python3 ../scripts/collapse_kmer_uniques.py --input_file results/p_7mers.csv --outfile results/p_7mers_uniq.csv --return_ids


python3 ../scripts/generate_deBruijn_kmers.py --input_fasta data/Polyurethenases.fasta --kmer_size 7 --outfile results/Polyurethenases_7mers.csv
python3 ../scripts/collapse_kmer_uniques.py --input_file results/Polyurethenases_7mers.csv --outfile results/Polyurethenases_7mers_unique.csv --return_ids