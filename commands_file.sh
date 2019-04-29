## extract kmers for desired subset(s) of proteins from large kmer file
while read exp; do echo "grep -G \"${exp}[0-9]\" /stor/work/Marcotte/project/rmcox/deBruijn_protein_maps/node_files/oma-seqs_nodes_10mer_allproteins.csv > ${exp}_10mers.csv"; done < exp_list.txt > grep_extract_kmers_1by1.sh

## format kmer file from csv to white-space separated (markov clustering pipeline)
while read exp; do echo "sed 's/,/ /g' /stor/work/Marcotte/project/rmcox/deBruijn_protein_maps/eukarya_analysis/individual_species/${exp}_10mers.csv > /stor/work/Marcotte/project/rmcox/deBruijn_protein_maps/eukarya_analysis/markov_clustering/${exp}_10mers_formatted.csv"; done < exp_list.txt > format_kmers_for_mcl.sh

## initialize markov clustering on kmer files (markov clustering pipeline)
while read exp; do echo "nohup mcl /stor/work/Marcotte/project/rmcox/deBruijn_protein_maps/eukarya_analysis/markov_clustering/${exp}_10mers_formatted.csv -I 7.0 -o /stor/work/Marcotte/project/rmcox/deBruijn_protein_maps/eukarya_analysis/markov_clustering/${exp}_10mers_mcl.csv --abc &> nohup_${exp}_10mers_mcl.out&"; done < exp_list.txt > run_mcl_1by1.sh

## collapse extracted kmers into uniques, while aggregating associated proteins (3rd column) and counting their number (4th column)
while read exp; do echo "python /stor/work/Marcotte/project/rmcox/scripts/collapse_uniques.py --input_file ${exp}_10mers.csv"; done < exp_list.txt > collapse_uniques_1by1.sh

## alphabetize kmers and then collapse into uniques (LGL pipeline)
while read exp; do echo "python /stor/work/Marcotte/project/rmcox/scripts/alphabetize_collapse_uniques.py --input_file ${exp}_10mers.csv"; done < exp_list.txt > alph_collapse_uniques_1by1.sh

## pull kmer columns from alphabetized/uniqued file where column 1 =/= column 2, then erase the header (LGL pipeline)
while read exp; do echo "awk -F',' '\$1!=\$2{print \$1,\$2}' ${exp}_10mers_alph_uniques.csv | tail -n +2 > ${exp}_10mers_lgl.ncol"; done < exp_list.txt > csv_to_ncol_1by1.sh

## create conf files appropriately edited for each experiment (LGL pipeline)
while read exp; do echo "sed 's/exp/${exp}_10mers/g' conf_file_template > conf_file_${exp}_10mers"; done < exp_list.txt > conf_file_edit_1by1.sh

## initialize LGL script for each experiment (LGL pipeline)
while read exp; do echo "nohup /stor/work/Marcotte/project/rmcox/programs/LGL/bin/lgl.pl -c conf_file_${exp}_10mers &> nohup_${exp}_10mers.out&"; done < exp_list.txt > nohup_generateLGLs_1by1.sh



