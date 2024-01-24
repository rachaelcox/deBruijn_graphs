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

# oligo 1 seq
AAGACCTACCTCACATGGCCAACACTCGGACAAAAAAAAAA

# oligo 2 seq
GACCACCAGCAGCAGCCAGCCGACGCAGGGACAAAAAAAAAA

# oligo 3 seq
AGAACTTACATCAACTAAACAACAAATGAACAAAAAAAAAA

# old data dir
/stor/work/Ellington/scratch/Zac/PLA/data/new_data/*Lib12*
/stor/work/Ellington/scratch/Zac/PLA/data/new_data/*Probe12* # same as Lib12?
/stor/work/Ellington/scratch/Zac/PLA/data/new_data/*Probe34* # what are these?
/stor/work/Ellington/scratch/Zac/PLA/data/new_data/*Oligo3-Undetermined* # what is this?

# organization
mkdir -p data/raw/{oligo_1,oligo_2,oligo_3,oligo_b,streptavidin}
mkdir -p data/qc/{oligo_1,oligo_2,oligo_3,oligo_b,streptavidin}
mkdir -p data/processed/{oligo_1,oligo_2,oligo_3,oligo_b,streptavidin}


ln -s /stor/work/Ellington/scratch/Zac/PLA/data/new_data/*Oligo1* /stor/work/Marcotte/project/rmcox/deBruijn_graphs/oligo/data/raw/oligo_1
ln -s /stor/work/Ellington/scratch/Zac/PLA/data/new_data/*Oligo2* /stor/work/Marcotte/project/rmcox/deBruijn_graphs/oligo/data/raw/oligo_2
ln -s /stor/work/Ellington/scratch/Zac/PLA/data/new_data/*Oligo3* /stor/work/Marcotte/project/rmcox/deBruijn_graphs/oligo/data/raw/oligo_3
ln -s /stor/work/Ellington/scratch/Zac/PLA/data/new_data/*OligoBlind* /stor/work/Marcotte/project/rmcox/deBruijn_graphs/oligo/data/raw/oligo_b

# streptavidin
ln -s /stor/work/Ellington/scratch/Zac/PLA/data/new_data/*Lib12-SA* /stor/work/Marcotte/project/rmcox/deBruijn_graphs/oligo/data/raw/streptavidin
ln -s /stor/work/Ellington/scratch/Zac/PLA/data/new_data/*Probe12-SA* /stor/work/Marcotte/project/rmcox/deBruijn_graphs/oligo/data/raw/streptavidin

# exp set up: oligo|wash|rep
# naming scheme is not consistent T_T

# =======================================
# data QC
# =======================================

# i'm going to try fastp, seems like a good all in 1 tool
# template command for fastp
fastp -i R1.fastq.qz -I R2.fastq.gz -o R1.trimmed.fastq.gz -o R2.trimmed.fastq.gz -p --failed_out data/qc/fastp_failed/

# oligo data
for x in 1 2 3 b; do ls -C data/raw/oligo_${x}/*R1* > data/meta/oligo_${x}_R1.txt; done
for x in 1 2 3 b; do ls -C data/raw/oligo_${x}/*R2* > data/meta/oligo_${x}_R2.txt; done
for x in 1 2 3 b; do paste data/meta/oligo_${x}_R1.txt data/meta/oligo_${x}_R2.txt > data/meta/oligo_${x}_R1R2.txt; done

# streptavidin data
ls -C data/raw/streptavidin/*R1* > data/meta/strep_R1.txt
ls -C data/raw/streptavidin/*R2* > data/meta/strep_R2.txt
paste data/meta/strep_R1.txt data/meta/strep_R2.txt > data/meta/strep_R1R2.txt


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

# finalized loop for oligos
for x in 1 2 3 b; do
	while IFS=$'\t' read -r -a y; do
    	o1="$(basename ${y[0]} .fastq.gz)"
    	o2="$(basename ${y[1]} .fastq.gz)"
    	merg="${o1%_L00*}.merged.fpqc.fastq.gz"
    	echo "fastp -i ${y[0]} -I ${y[1]} -o data/qc/oligo_${x}/${o1}.fpqc.fastq.gz -O data/qc/oligo_${x}/${o2}.fpqc.fastq.gz --failed_out data/qc/fastp/failed/${o1%_L00*}.failed.txt --merge --merged_out data/processed/oligo_${x}/${merg} --length_limit 100 --json data/qc/fastp/oligo_${x}_${o1%_L00*}.json --html data/qc/fastp/oligo_${x}_${o1%_L00*}.html"
	done < data/meta/oligo_${x}_R1R2.txt
done > cmds/fastp.cmds
cat cmds/fastp.cmds | parallel -j24

# finalized loop for streptavidin
while IFS=$'\t' read -r -a y; do
	s1="$(basename ${y[0]} .fastq.gz)"
	s2="$(basename ${y[1]} .fastq.gz)"
	merg="${s1%_L00*}.merged.fpqc.fastq.gz"
	echo "fastp -i ${y[0]} -I ${y[1]} -o data/qc/streptavidin/${s1}.fpqc.fastq.gz -O data/qc/streptavidin/${s2}.fpqc.fastq.gz --failed_out data/qc/fastp/failed/${s1%_L00*}.failed.txt --merge --merged_out data/processed/streptavidin/${merg} --length_limit 100 --json data/qc/fastp/streptavidin_${s1%_L00*}.json --html data/qc/fastp/streptavidin_${s1%_L00*}.html"
done < data/meta/strep_R1R2.txt > cmds/fastp_strep.cmds
cat cmds/fastp_strep.cmds | parallel -j5

# qc trimmed and merged reads for oligos
for x in 1 2 3 b; do echo "fastqc --noextract --nogroup --threads 16 -o data/qc/oligo_${x}_merged data/processed/oligo_${x}/*"; done > cmds/fastq_fastp.cmds
cat cmds/fastq_fastp.cmds | parallel -j4

for x in 1 2 3 b; do echo "multiqc data/qc/oligo_${x}_merged/* --filename data/qc/multiqc/oligo_${x}_merged --interactive --force"; done > cmds/multiqc_fastp.cmds
cat cmds/multiqc_fastp.cmds | parallel -j4

# qc trimmed and merged reads for streptavidin
fastqc --noextract --nogroup --threads 16 -o data/qc/streptavidin_merged data/processed/streptavidin_pe_merged/*
multiqc data/qc/streptavidin_merged/* --filename data/qc/multiqc/streptavidin_merged --interactive --force

# reports:
# oligo1: file:///H:/project/rmcox/deBruijn_graphs/oligo/data/qc/multiqc/oligo_1_merged.html
# oligo2: file:///H:/project/rmcox/deBruijn_graphs/oligo/data/qc/multiqc/oligo_2_merged.html
# oligo3: file:///H:/project/rmcox/deBruijn_graphs/oligo/data/qc/multiqc/oligo_3_merged.html
# oligob: file:///H:/project/rmcox/deBruijn_graphs/oligo/data/qc/multiqc/oligo_b_merged.html
# streptavidin: file:///H:/project/rmcox/deBruijn_graphs/oligo/data/qc/multiqc/streptavidin_merged.html

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
mkdir data/processed/streptavidin_wash_merged

# oligo files
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

# streptavidin files
ls data/processed/streptavidin_pe_merged/*Probe* > data/meta/streptavidin_bg.txt
ls data/processed/streptavidin_pe_merged/*0w* > data/meta/streptavidin_0w.txt
ls data/processed/streptavidin_pe_merged/*3w* > data/meta/streptavidin_3w.txt
ls data/processed/streptavidin_pe_merged/*7w* > data/meta/streptavidin_7w.txt

for exp in bg 0w 3w 7w; do
	echo "{ xargs cat < data/meta/streptavidin_${exp}.txt ; } > data/processed/streptavidin_wash_merged/streptavidin_${exp}.merged.fastq.gz"
done > cmds/concat_merged_strep_reads.sh
cat cmds/concat_merged_strep_reads.sh | parallel -j5

# decompress
for x in 1 2 3 b; do gunzip -f data/processed/oligo_${x}_wash_merged/*gz; done
gunzip -f data/processed/streptavidin_wash_merged/*gz

# convert to fasta format
for x in 1 2 3 b; do
	for y in bg 0w 1w 3w; do
	sed -n '1~4s/^@/>/p;2~4p' data/processed/oligo_${x}_wash_merged/oligo_${x}_${y}.merged.fastq > data/processed/oligo_${x}_wash_merged/oligo_${x}_${y}.merged.fasta
	done
done

for exp in bg 0w 3w 7w; do
	sed -n '1~4s/^@/>/p;2~4p' data/processed/streptavidin_wash_merged/streptavidin_${exp}.merged.fastq > data/processed/streptavidin_wash_merged/streptavidin_${exp}.merged.fasta
done

# =======================================
# prune seqs to match N15 flanking regions
# =======================================

n=2
for x in 1 2 3 b; do
	for y in bg 0w 1w 3w; do
	mkdir -p data/pruned/oligo_${x}_${n}mm/
	echo "python3 ../scripts/prune_seqs.py --fasta_file data/processed/oligo_${x}_wash_merged/oligo_${x}_${y}.merged.fasta --outfile data/pruned/oligo_${x}_${n}mm/oligo_${x}_${y}.merged.${n}mm.fasta --num_mismatch ${n} | tee -a logs/prune_oligo_${x}_${y}_${n}mm.log"
done
done > cmds/prune_reads_${n}mm.sh

# =======================================
# generate de Bruijn graphs for oligos
# =======================================

# generate all kmers
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

# collapse to weighted graph
for k in {5..15}; do
	for x in 1 2 3 b; do
	for y in bg 0w 1w 3w; do
	echo "python3 ../scripts/collapse_kmer_uniques.py --input_file results/${k}mer/oligo_${x}_${y}_${k}mers.csv --outfile results/${k}mer/kmer_counts_all/oligo_${x}_${y}_${k}mers_unique.csv | tee -a logs/oligo_${x}_${y}_${k}mers_collapse_unique.log"
	done
done > cmds/generate_${k}mer_weights.sh
done 

rm cmds/generate_5to15mer_weights.sh
touch cmds/generate_5to15mer_weights.sh
for k in {5..15}; do
	cat cmds/generate_${k}mer_weights.sh >> cmds/generate_5to15mer_weights.sh
done
cat cmds/generate_5to15mer_weights.sh | parallel -j6
tail -n +2 cmds/generate_5to15mer_weights.sh | parallel -j6

# uggghhh clay is still hogging hfog1
# going to tee this up on hfog2
cd /project/rmcox/deBruijn_graphs/oligo
for k in {5..15}; do
	mkdir -p results_hfog2/${k}mer/kmer_counts_all
done

for k in {5..15}; do
	for x in 1 2 3 b; do
	for y in bg 0w 1w 3w; do
	echo "python3 ../scripts/collapse_kmer_uniques.py --input_file results/${k}mer/oligo_${x}_${y}_${k}mers.csv --outfile results_hfog2/${k}mer/kmer_counts_all/oligo_${x}_${y}_${k}mers_unique.csv | tee -a logs/hfog2/oligo_${x}_${y}_${k}mers_collapse_unique.log"
	done
done > cmds/generate_${k}mer_weights_hfog2.sh
done 

rm cmds/generate_5to15mer_weights_hfog2.sh
touch cmds/generate_5to15mer_weights_hfog2.sh
for k in {5..15}; do
	cat cmds/generate_${k}mer_weights_hfog2.sh >> cmds/generate_5to15mer_weights_hfog2.sh
done
cat cmds/generate_5to15mer_weights_hfog2.sh | parallel -j12

# wtffff the script is terminating without an error
# (╯°□°)╯︵ ┻━┻
head -1000 results/5mer/oligo_1_0w_5mers.csv > debug/collapse_debug.csv
python3 ../scripts/collapse_kmer_uniques.py --input_file debug/collapse_debug.csv --outfile debug/collapse_debug_unique.csv

# ok I fixed it
# 5 to 9mers on hfog1, going to start at the 10mers on hfog2
cat cmds/generate_5to15mer_weights.sh | grep -n '10mer' # starts at line 81
head -80 cmds/generate_5to15mer_weights.sh > cmds/generate_5to9mer_weights.sh
cat cmds/generate_5to9mer_weights.sh | parallel -j2
tail -n +81 cmds/generate_5to15mer_weights_hfog2.sh > cmds/generate_10to15mer_weights_hfog2.sh
cat cmds/generate_10to15mer_weights_hfog2.sh | parallel -j8

# old notes from first pass:
python3 ../scripts/collapse_kmer_uniques.py --input_file results/10mer/oligo_1_0w_10mers.csv --outfile results/10mer/oligo_1_0w_10mers_unique.csv
python3 ../scripts/collapse_kmer_uniques.py --input_file results/10mer/oligo_b_3w_10mers.pkl --outfile results/10mer/oligo_b_3w_10mers_unique.csv

python3 ../scripts/calc_consensus_seq.py --input_file results/15mer/oligo_b_3w_15mers_unique.csv

python3 ../scripts/calc_consensus_seq.py --input_file results/15mer/oligo_1_3w_15mers_unique.csv

# yeah this does not work with k=10, highly redundant sections of the N15 library are making it too difficult for the consensus algorithm to find a sequence potentially?

# =======================================
# generate de Bruijn graphs for oligos (pruned)
# =======================================

# generate all kmers
for k in {5..15}; do
	for x in 1 2 3 b; do
	for y in bg 0w 1w 3w; do
	echo "python3 ../scripts/generate_deBruijn_kmers.py --input_fasta data/pruned/oligo_${x}_1mm/oligo_${x}_${y}.merged.1mm.fasta --kmer_size $k --outfile results_pruned/${k}mer/oligo_${x}_${y}_1mm_${k}mers.csv"
	done
done > cmds/generate_${k}mer_dbgs_1mm.sh
done

rm cmds/generate_5to15mer_dbgs_1mm.sh
touch cmds/generate_5to15mer_dbgs_1mm.sh
for k in {5..15}; do
	cat cmds/generate_${k}mer_dbgs_1mm.sh >> cmds/generate_5to15mer_dbgs_1mm.sh
done
cat cmds/generate_5to15mer_dbgs_1mm.sh | parallel -j16

# collapse to weighted graph
for k in {5..15}; do
	for x in 1 2 3 b; do
	for y in bg 0w 1w 3w; do
	echo "python3 ../scripts/collapse_kmer_uniques.py --input_file results_pruned/${k}mer/oligo_${x}_${y}_1mm_${k}mers.csv --outfile results_pruned/${k}mer/oligo_${x}_${y}_1mm_${k}mers_unique.csv | tee -a logs/oligo_${x}_${y}_${k}mers_collapse_unique_1mm.log"
	done
done > cmds/generate_${k}mer_weights_1mm.sh
done 

rm cmds/generate_5to15mer_weights_1mm.sh
touch cmds/generate_5to15mer_weights_1mm.sh
for k in {5..15}; do
	cat cmds/generate_${k}mer_weights_1mm.sh >> cmds/generate_5to15mer_weights_1mm.sh
done
cat cmds/generate_5to15mer_weights_1mm.sh | parallel -j16

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

# =======================================
# generate de Bruijn graphs for streptavidin
# =======================================

# generate all kmers
for k in {5..15}; do
	for exp in bg 0w 3w 7w; do
	echo "python3 ../scripts/generate_deBruijn_kmers.py --input_fasta data/processed/streptavidin_wash_merged/streptavidin_${exp}.merged.fasta --kmer_size $k --outfile results/${k}mer/streptavidin_${exp}_${k}mers.csv"
done > cmds/generate_${k}mer_dbgs_strep.sh
done

rm cmds/generate_5to15mer_dbgs_strep.sh
touch cmds/generate_5to15mer_dbgs_strep.sh
for k in {5..15}; do
	cat cmds/generate_${k}mer_dbgs_strep.sh >> cmds/generate_5to15mer_dbgs_strep.sh
done
cat cmds/generate_5to15mer_dbgs_strep.sh | parallel -j24


# collapse to weighted graph
# going to try this with cuda's new pandas GPU parallelization
# gotta be on hfog2 or hfog3 for this
for k in {5..15}; do
	for exp in bg 0w 3w 7w; do
	echo "python3 -m cudf.pandas ../scripts/collapse_kmer_uniques.py --input_file results/${k}mer/streptavidin_${exp}_${k}mers.csv --outfile results/${k}mer/streptavidin_${exp}_${k}mers_unique.csv | tee -a logs/streptavidin_${exp}_${k}mers_collapse_unique.log"
done > cmds/generate_${k}mer_weights_strep.sh
done 

rm cmds/generate_5to15mer_weights_strep.sh
touch cmds/generate_5to15mer_weights_strep.sh
for k in {5..15}; do
	cat cmds/generate_${k}mer_weights_strep.sh >> cmds/generate_5to15mer_weights_strep.sh
done
cat cmds/generate_5to15mer_weights_strep.sh | parallel -j24

# wooooow i think the cudf pandas is actually slower
# starting a parallel run on hfog1
for k in {5..15}; do
	for exp in bg 0w 3w 7w; do
	echo "python3 ../scripts/collapse_kmer_uniques.py --input_file results/${k}mer/streptavidin_${exp}_${k}mers.csv --outfile results/kmer_counts_all/streptavidin_${exp}_${k}mers_unique.csv | tee -a logs/streptavidin_${exp}_${k}mers_collapse_unique.log"
done > cmds/generate_${k}mer_weights_strep.sh
done 

rm cmds/generate_5to15mer_weights_strep.sh
touch cmds/generate_5to15mer_weights_strep.sh
for k in {5..15}; do
	cat cmds/generate_${k}mer_weights_strep.sh >> cmds/generate_5to15mer_weights_strep.sh
done
cat cmds/generate_5to15mer_weights_strep.sh | parallel -j4

# make "all" oligo bait kmer files (each oligo concatenated) for enrichment script
for k in {5..15}; do
	rm oligo_all_${k}mers.csv
	touch oligo_all_${k}mers.csv
	echo "Node1,Node2,ProteinID" >> oligo_all_${k}mers.csv
	for x in 1 2 3; do
	tail -n +2 -q oligo_${x}_${k}mers.csv >> oligo_all_${k}mers.csv
done
done

# make files for top 100 enriched kmers from each sample for discussions with sanchita
for x in `ls results/summarized/*10mers_norm*`; do 
	f=${x##*/}
	rm results/debug/${f%%.csv}_top100.csv
	touch results/debug/${f%%.csv}_top100.csv
	echo "node1,node2,count,label,exp,bg_count,type,exp_pseudo,ctrl_pseudo,F0ctrl_pseudo,F0exp_pseudo,fc,log2fc,F0_ctrl,F0_exp,F1,zscore" >> results/debug/${f%%.csv}_top100.csv
	sort -t, -nrk13 ${x} | head -100 >> results/debug/${f%%.csv}_top100.csv
done
scp -r rmcox@hfogcomp01.ccbb.utexas.edu:/stor/work/Marcotte/project/rmcox/deBruijn_graphs/oligo/results/debug .

# =======================================
# threshold kmer counts by quantile
# =======================================

# results dirs getting messy, time to clean up
# ∴ some of the previous paths might be anachronous

# reorganize un-thresholded files to their own directory
for k in {5..15}; do
	mkdir oligo/results/${k}mer/kmer_counts_all/
	mv oligo/results/${k}mer/*${k}mers_unique.csv oligo/results/${k}mer/kmer_counts_all/
done

# set up directories for thresholded kmers
for k in {5..15}; do
	mkdir oligo/results/${k}mer/kmer_counts_q03cut/
	mkdir oligo/results/${k}mer/kmer_counts_q04cut/
	mkdir oligo/results/${k}mer/kmer_counts_q05cut/
done

for k in {5..15}; do
	ls -1 results/${k}mer/kmer_counts_all/*unique.csv | xargs -n1 basename > file_paths/${k}mer_files.txt
done

for k in {5..15}; do
	while read infile; do
	q="q03"
	echo "python3 ../scripts/threshold_kmer_counts.py --infile results/${k}mer/kmer_counts_all/${infile} --outfile results/${k}mer/kmer_counts_${q}cut/${infile%%unique.csv}${q}.csv --cutoff 0.03 | tee -a logs/${q}cut/threshold_${infile%%unique.csv}${q}.log"
done < file_paths/${k}mer_files.txt
done > cmds/threshold_kmers_q03.cmds
cat cmds/threshold_kmers_q03.cmds | parallel -j8

for k in {5..15}; do
	while read infile; do
	q="q04"
	echo "python3 ../scripts/threshold_kmer_counts.py --infile results/${k}mer/kmer_counts_all/${infile} --outfile results/${k}mer/kmer_counts_${q}cut/${infile%%unique.csv}${q}.csv --cutoff 0.04 | tee -a logs/${q}cut/threshold_${infile%%unique.csv}${q}.log"
done < file_paths/${k}mer_files.txt
done > cmds/threshold_kmers_q04.cmds
cat cmds/threshold_kmers_q04.cmds | parallel -j8

for k in {5..15}; do
	while read infile; do
	q="q05"
	echo "python3 ../scripts/threshold_kmer_counts.py --infile results/${k}mer/kmer_counts_all/${infile} --outfile results/${k}mer/kmer_counts_${q}cut/${infile%%unique.csv}${q}.csv --cutoff 0.05 | tee -a logs/${q}cut/threshold_${infile%%unique.csv}${q}.log"
done < file_paths/${k}mer_files.txt
done > cmds/threshold_kmers_q05.cmds
cat cmds/threshold_kmers_q05.cmds | parallel -j8

# WOOOOW someone is using all the threads on hfog1
# going to have to run these serially
bash cmds/threshold_kmers_q03.cmds
bash cmds/threshold_kmers_q04.cmds
bash cmds/threshold_kmers_q05.cmds

# =======================================
# calculate kmer enrichment & filter
# =======================================

# [08/2023] accomplished in jupyter/R notebooks:
# notebooks/oligo_kmer_enrichment.ipynb
# scripts/plot_oligo_enrichment.R


# organize enrichment files
for k in {5..15}; do
	mkdir oligo/results/${k}mer/enrichment_all/
	
	mkdir oligo/results/${k}mer/enrichment_q03cut/
	mkdir oligo/results/${k}mer/enrichment_q04cut/
	mkdir oligo/results/${k}mer/enrichment_q05cut/
done

for k in {5..15}; do
	mv results/${k}mer/enrichment_all/*${k}mer_enrichment.csv archive/
done

# [01/2024] rewrite to CLI

# template command
python3 ../scripts/calc_kmer_enrichment.py \
--exp_file results/10mer/kmer_counts_all/oligo_3_3w_10mers_unique.csv \
--bg_file results/10mer/kmer_counts_all/oligo_3_bg_10mers_unique.csv \
--lib_file data/bait_kmers/library_10mers.csv \
--oligo_file data/bait_kmers/oligo_3_comp_10mers.csv \
--outfile results/10mer/enrichment_all/oligo_3_3w_enrich.csv

# generate cmds in loop for oligos
for k in {5..15}; do
	for x in 1 2 3; do
	echo "python3 ../scripts/calc_kmer_enrichment.py \
--exp_file results/${k}mer/kmer_counts_all/oligo_${x}_0w_${k}mers_unique.csv \
--bg_file results/${k}mer/kmer_counts_all/oligo_${x}_bg_${k}mers_unique.csv \
--lib_file data/bait_kmers/library_${k}mers.csv \
--oligo_file data/bait_kmers/oligo_${x}_comp_${k}mers.csv \
--outfile results/${k}mer/enrichment_all/oligo_${x}_0w_enrich.csv"
	echo "python3 ../scripts/calc_kmer_enrichment.py \
--exp_file results/${k}mer/kmer_counts_all/oligo_${x}_1w_${k}mers_unique.csv \
--bg_file results/${k}mer/kmer_counts_all/oligo_${x}_bg_${k}mers_unique.csv \
--lib_file data/bait_kmers/library_${k}mers.csv \
--oligo_file data/bait_kmers/oligo_${x}_comp_${k}mers.csv \
--outfile results/${k}mer/enrichment_all/oligo_${x}_1w_enrich.csv"
	echo "python3 ../scripts/calc_kmer_enrichment.py \
--exp_file results/${k}mer/kmer_counts_all/oligo_${x}_3w_${k}mers_unique.csv \
--bg_file results/${k}mer/kmer_counts_all/oligo_${x}_bg_${k}mers_unique.csv \
--lib_file data/bait_kmers/library_${k}mers.csv \
--oligo_file data/bait_kmers/oligo_${x}_comp_${k}mers.csv \
--outfile results/${k}mer/enrichment_all/oligo_${x}_3w_enrich.csv"
done
done > cmds/calc_kmer_enrichment_oligos.sh
cat cmds/calc_kmer_enrichment_oligos.sh | parallel -j12

# generate cmds in loop for streptavidin
for k in {5..15}; do
	for w in 0w 3w 7w; do
	echo "python3 ../scripts/calc_kmer_enrichment.py \
--exp_file results/${k}mer/kmer_counts_all/streptavidin_${w}_${k}mers_unique.csv \
--bg_file results/${k}mer/kmer_counts_all/streptavidin_bg_${k}mers_unique.csv \
--lib_file data/bait_kmers/library_${k}mers.csv \
--oligo_file data/bait_kmers/oligo_all_comp_${k}mers.csv \
--outfile results/${k}mer/enrichment_all/streptavidin_${w}_enrich.csv"
done
done > cmds/calc_kmer_enrichment_strep.sh
cat cmds/calc_kmer_enrichment_strep.sh | parallel -j8

# plot results
# template command
Rscript ../scripts/plot_kmer_enrichment.R --kmer_size 6
for k in {5..15}; do
	echo "Rscript ../scripts/plot_kmer_enrichment.R --kmer_size ${k}"
done > cmds/plot_kmer_enr.sh
cat cmds/plot_kmer_enr.sh | parallel -j6

# filter kmers based on apex z-score
# (look in jupyter)

# generate PCA plots
# template command
Rscript ../scripts/plot_pca.R --kmer_size 6 --top_n 100
for k in {5..15}; do
	echo "Rscript ../scripts/plot_pca.R --kmer_size ${k} --top_n 100"
done > cmds/plot_kmer_pca.sh
cat cmds/plot_kmer_pca.sh | parallel -j3

# ---------------------------------------
# generate consensus sequences from filtered kmers
# ---------------------------------------

for k in {5..15}; do
	mkdir results/${k}mer/consensus_seqs/
done

for k in {5..15}; do
	python3 ../scripts/calc_consensus_seq.py --input_file results/${k}mer/enrichment_all/streptavidin_3w_enrich_fdr_filter.csv --seed_diversity
	python3 ../scripts/calc_consensus_seq.py --input_file results/${k}mer/enrichment_all/streptavidin_3w_enrich_z_filter.csv --seed_diversity
	python3 ../scripts/calc_consensus_seq.py --input_file results/${k}mer/enrichment_all/streptavidin_7w_enrich_fdr_filter.csv --seed_diversity
	python3 ../scripts/calc_consensus_seq.py --input_file results/${k}mer/enrichment_all/streptavidin_7w_enrich_z_filter.csv --seed_diversity
	for i in 1 2 3; do
	python3 ../scripts/calc_consensus_seq.py --input_file results/${k}mer/enrichment_all/oligo_${i}_3w_enrich_fdr_filter.csv --seed_diversity
	python3 ../scripts/calc_consensus_seq.py --input_file results/${k}mer/enrichment_all/oligo_${i}_3w_enrich_z_filter.csv --seed_diversity
done
done

for x in 1 2 3; do
	rm results/oligo_${x}_cons_seqs.fasta
done
for k in 5 8 10; do
	cat results/${k}mer/enrichment_all/oligo_1*z_filter*fasta >> results/oligo_1_cons_seqs.fasta
	cat results/${k}mer/enrichment_all/oligo_2*z_filter*fasta >> results/oligo_2_cons_seqs.fasta
	cat results/${k}mer/enrichment_all/oligo_3*z_filter*fasta >> results/oligo_3_cons_seqs.fasta
done

# ---------------------------------------
# generate 100 consensus sequences from filtered kmers
# ---------------------------------------
k=10
for i in 1 2 3; do
	echo "python3 ../scripts/calc_consensus_seq.py --input_file results/thresholded/oligo_${i}_${k}mer_graph_fdr1e-10.csv --calc_n 100 --seed_diversity --outfile_prefix results/reconstructions/oligo_${i}_${k}mer_graph_fdr1e-10_sd_consensus | tee -a logs/generate_100_oligo_${i}_${k}mer_fdr1e-10_consensus_seqs_sd.log"
done > cmds/generate_100_${k}mer_consensus_seqs_sd.sh
cat cmds/generate_100_${k}mer_consensus_seqs_sd.sh | parallel -j3

k=10
for i in 1 2 3; do
	echo "python3 ../scripts/calc_consensus_seq.py --input_file results/thresholded/oligo_${i}_${k}mer_graph_fdr1e-10_3log2fc.csv --calc_n 100 --seed_diversity --outfile_prefix results/reconstructions/oligo_${i}_${k}mer_graph_fdr1e-10_3log2fc_sd_consensus | tee -a logs/generate_100_oligo_${i}_${k}mer_fdr1e-10_3log2fc_sd_consensus_seqs.log"
done > cmds/generate_100_${k}mer_consensus_seqs_3log2fc_sd.sh
cat cmds/generate_100_${k}mer_consensus_seqs_3log2fc_sd.sh | parallel -j3


# --------------------------------------------------
# MISCELLANEOUS
# --------------------------------------------------

# polyketides
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

# guaymas nylonases and petases
python3 scripts/generate_deBruijn_kmers.py --input_fasta guaymas/data/NylonHits_Guaymas2020_ALLBINS.fasta --kmer_size 8 --outfile guaymas/results/NylonHits_7mers.csv
python3 scripts/collapse_kmer_uniques.py --input_file guaymas/results/NylonHits_7mers.csv --outfile guaymas/results/NylonHits_7mers_unique.csv
python3 scripts/calc_consensus_seq.py --input_file guaymas/results/NylonHits_7mers_unique.csv --outfile_prefix guaymas/results/NylonHits_consensus
python3 scripts/generate_deBruijn_kmers.py --input_fasta guaymas/results/NylonHits_consensus.fasta --kmer_size 8 --outfile guaymas/results/NylonHits_consensus_7mers.csv

python3 scripts/generate_deBruijn_kmers.py --input_fasta guaymas/data/PETHits_Guaymas2020_ALLBINS.fasta --kmer_size 8 --outfile guaymas/results/PETHits_7mers.csv
python3 scripts/collapse_kmer_uniques.py --input_file guaymas/results/PETHits_7mers.csv --outfile guaymas/results/PETHits_7mers_unique_ids.csv --return_ids
python3 scripts/calc_consensus_seq.py --input_file guaymas/results/PETHits_7mers_unique.csv --outfile_prefix guaymas/results/PETHits_consensus
python3 scripts/generate_deBruijn_kmers.py --input_fasta guaymas/results/PETHits_consensus.fasta --kmer_size 8 --outfile guaymas/results/PETHits_consensus_7mers.csv

# phylo-stratified nylonase set
for i in 0 1 2 3 4; do
	python3 scripts/generate_deBruijn_kmers.py --input_fasta /stor/work/Marcotte/project/drbarth/plastics/data/raw/deBruijn_g${i}.fasta --kmer_size 7 --outfile /stor/work/Marcotte/project/rmcox/deBruijn_graphs/guaymas/results/phylo_stratified/nylonase_g${i}_7mers.csv
done

for i in 0 1 2 3 4; do
	python3 scripts/collapse_kmer_uniques.py --input_file guaymas/results/phylo_stratified/nylonase_g${i}_7mers.csv --outfile guaymas/results/phylo_stratified/nylonase_g${i}_7mers_unique.csv
done

# number of sequences per clade:
for i in 0 1 2 3 4; do
	grep -c '^>' /stor/work/Marcotte/project/drbarth/plastics/data/raw/deBruijn_g${i}.fasta
done

# 0: 31
# 1: 62
# 2: 26
# 3: 181
# 4: 206

# top shared per clade:
for i in 0 1 2 3 4; do
	head guaymas/results/phylo_stratified/nylonase_g${i}_7mers_unique.csv
done

# 0: EQGTSF,QGTSFD  8/31 (25.8%)
# 1: ALIGIA,LIGIAI  13/62 (21.0%)
# 2: PLGMTD,LGMTDT  8/26 (30.8%)
# 3: NALVGI,ALVGIL  58/181 (32.0%)
# 4: ILTRRS,LTRRSG  34/206 (16.5%)

# number edges per weighted clade graph
for i in 0 1 2 3 4; do
	wc -l guaymas/results/phylo_stratified/nylonase_g${i}_7mers_unique.csv
done

# --------------------------------------------------
# TO DO
# --------------------------------------------------

# make 100 reconstructions starting at top 100 kmers

# kmers enriched in all oligos?? (no match distribution)

