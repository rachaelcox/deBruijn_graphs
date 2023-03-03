from Bio import Entrez, SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import pandas as pd
import argparse

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Fetches annotations for list of GenBank IDs from the NCBI Nucleotide Database")
    parser.add_argument("-e", "--email", action="store", required=True,
                                        help="Filename containing delimited list of GenBank IDs")
    parser.add_argument("-i", "--id_file", action="store", required=True,
                                        help="Filename containing delimited list of GenBank IDs")
    parser.add_argument("-o", "--outfile", action="store", required=False,
                                        help="Optional name for outfile (default = gb_annots.txt)")
    args = parser.parse_args()

def fetch_gb_annots(email, id_file):

	Entrez.email = str(email) # given email

	# format ID list
	with open(id_file, 'r') as f:
		id_list = f.read().splitlines()

	# download record for each ID
	handle = Entrez.efetch(db="protein", id=id_list, rettype="gb", retmode="text")
	records = SeqIO.parse(handle, "genbank")

	# initialize list of dictionaries to hold annotation info
	annots = []

	# loop over all features in the record
	for record in records:

		# initialize list(s) to hold info
		record_dict = {}
		note_list = []

		# extract protein information where available
		if 'db_source' in record.annotations:
			if 'pdb' in record.annotations['db_source']:
				prot_annot = 'synthetic construct'

		for feature in record.features:
			if feature.type == 'Protein':
				prot_annot = feature.qualifiers['product'][0]
			if feature.type == 'Region':
				note = feature.qualifiers['note'][0]
				note_list.append(note)

		prot_seq = str(record.seq)
		prot_analyzed = ProteinAnalysis(prot_seq) # this does not handle seqs with 'X'
		# going to have to come back later to include protein seq analysis in annot file

		record_dict['sequence'] = prot_seq  # protein sequence
		#record_dict['mw'] = prot_analyzed.molecular_weight()  # molecular weight of protein
		record_dict['notes'] = note_list
		record_dict['organism'] = record.annotations['source']  # organism
		record_dict['protein'] = prot_annot  # protein product name
		record_dict['accession'] = record.id  # accession number

		
		annots.append(record_dict)

	print(record_dict)
	all_annots = pd.DataFrame(annots)

	print(all_annots.head())
	print(all_annots.tail())


	all_annots.to_csv('test.csv', index = False)			

	handle.close()
	return

def write_outfile(outfile):

	if outfile == None:
		writefile = 'gb_annots.txt'
	else:
		writefile = args.outfile

fetch_gb_annots(args.email, args.id_file)

#write_outfile(args.outfile, annot_df)
