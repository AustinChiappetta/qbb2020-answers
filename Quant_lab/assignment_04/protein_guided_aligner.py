#!/usr/bin/env python 3
"""Usage: protein_guided_aligner.py <DNA.fa> <aligned_aa.fa>. Output is nucleotide alignment guided by protein alignment"""
import sys
from fasta_iterator_class import FASTAReader

#1. use FASTAReader to grab seq names and nucleotide seqs from seqdump file, save to nucleotides_parsed as list of tuples
def nuc_parser(DNA_seqs_fa):
	#read in sequence
	nucleotides = FASTAReader(open(DNA_seqs_fa))  #seqdump.txt
	#print('Parsing DNA sequences')
	#print('...')
	nucleotides_parsed = []
	for seq_id, seq in nucleotides:
		nucleotides_parsed += [(seq_id, seq)]
	#print('DNA parsing complete')
	return(nucleotides_parsed)


#2. use FASTAReader to grab seq names and aa seqs from aligned_aa_seqs file, save to aa_parsed as list of tuples
def aa_parser(aa_aligned_fa):
	#read in sequence
	amino_acids = FASTAReader(open(aa_aligned_fa)) #aligned_aa_seqs.txt
	#print('Parsing aligned protein sequences')
	#print('...')
	aa_parsed = []
	for seq_id, seq in amino_acids:
		aa_parsed += [(seq_id, seq)]
	#print('Protein parsing complete')
	return(aa_parsed)
    


#3. use zip to iterate through nucleotides_parsed list and aa_parsed list, making a new list of tuples: [("seq_id", "nuc_seq_with_gaps")]
def aa_guided_nuc_aligner(parsed_DNA, parsed_aa):
	#print('Performing protein-guided nucleotide alignment')
	#print('...')
	nucleotides_with_gaps = []
	for parsed in zip(parsed_DNA, parsed_aa):
		#initialize variables
		nuc_tuple = parsed[0]
		nuc_seq = nuc_tuple[1]
		nuc_ID = nuc_tuple[0]
		aa_tuple = parsed[1]
		aa_seq = aa_tuple[1]
		nuc_gapped_seq = ""
		
		for residue in aa_seq: #interate through aa seq, char by char.
			if residue != "-":
				nuc_gapped_seq += nuc_seq[0:3] #add codon to seq
				nuc_seq = nuc_seq[3:] #shift sequence to next codon

			elif residue == "-":
				nuc_gapped_seq += "---" #add gap

		nucleotides_with_gaps += [(nuc_ID, nuc_gapped_seq)]

	for seq_tuple in nucleotides_with_gaps:
		ID = seq_tuple[0]
		sequence = seq_tuple[1]
		#print(ID + '\n')
		print(ID)
		while len(sequence) > 60:
			#print(sequence[:60] + '\n')
			print(sequence[:60])
			sequence = sequence[60:]
		#print(sequence + '\n')
		print(sequence)


if __name__ == "__main__":
	nuc = sys.argv[1]
	aa = sys.argv[2]

	DNA_parsed = nuc_parser(nuc)
	prot_parsed = aa_parser(aa)
	aa_guided_nuc_aligner(DNA_parsed, prot_parsed)
