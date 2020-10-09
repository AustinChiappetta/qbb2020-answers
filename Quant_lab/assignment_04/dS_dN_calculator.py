#!/usr/bin/env python3

import sys
from fasta_iterator_class import FASTAReader
import math
import copy

#1. open aligned_DNA seq file with FASTAReader and make a list of tuples (ID, sequence)
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

#2. save first list entry (tuple of ID and sequence) as query
def get_query(parsed_DNA):
	query = parsed_DNA[0]
	return(query)

#3. save other list entries as targets
def get_targets(parsed_DNA):
	targets = parsed_DNA[1:]
	return(targets)

def dS_dN(query, targets):
	#4. initialize a mutations_dict dictionary of dS (synonymous) and dN (non-synonymous) changes
		#keys are codon position (1, 2, 3, etc.) and values are a sub dictionary:
			#for each codon, keys are 'dS' and 'dN'
			#values are number fo those changes for each codon
			#{1: {'dS': 23, 'dN': 42}, 2: {'dS': 11, 'dN': 18}, etc.}
	mutations_dict = {}

	codontable = {
	'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
	'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
	'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
	'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
	'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
	'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
	'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
	'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
	'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
	'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
	'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
	'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
	'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
	'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
	'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
	'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W'
	}

	DNA_nucleotides = 'ATCG'

	for t_ID_and_seq in targets: #5. iterate through target sequences
		t_seq = t_ID_and_seq[1]
		q_seq = query[1]
		for i in range(0, len(t_seq), 3):
			t_codon = t_seq[i:i+3]
			q_codon = q_seq[i:i+3]
			codon_num = int(i+3/3)

			if q_codon == "---":
				continue
			if t_codon == "---":
				continue

			if q_codon not in codontable.keys(): #skip codons with ambiguous nucleotides
				continue

			if t_codon not in codontable.keys(): #skip codons with ambiguous nucleotides
				continue

			if t_codon == q_codon: #skip identical codons
				continue
			else: #compare query codon and target codon by codon
				if codontable[t_codon] == codontable[q_codon]: #dS
					mutations_dict.setdefault(codon_num, {})
					mutations_dict[codon_num].setdefault('dS', 0)
					mutations_dict[codon_num]['dS'] += 1

				else: #dN
					mutations_dict.setdefault(codon_num, {})
					mutations_dict[codon_num].setdefault('dN', 0)
					mutations_dict[codon_num]['dN'] += 1

	#6. iterate through mutations_dict and calculate D = dN - dS for each key (codon), also calculate log2(dN/dS) for each
	codons = list(mutations_dict.keys())
	codons.sort()
	D = []
	log2_ratio_dN_dS = []
	for codon in codons:
		mutations_dict[codon].setdefault('dS', 0)
		mutations_dict[codon].setdefault('dN', 0)
		dN = mutations_dict[codon]['dN']
		dS = mutations_dict[codon]['dS']

		diff = dN - dS #calculate D
		#mutations_dict[key]['D'] = diff
		D += [[codon, diff]]

		if dN == 0:
			continue
		elif dS == 0:
			continue
		else:
			ratio = dN/dS
			log2_ratio = math.log(ratio, 2) 
			log2_ratio_dN_dS += [[codon, log2_ratio]]

	return(D, log2_ratio_dN_dS)


if __name__ == "__main__":
	aligned_DNA = sys.argv[1]

	DNA_parsed = nuc_parser(aligned_DNA)
	parsed_query = get_query(DNA_parsed)
	parsed_targets = get_targets(DNA_parsed)
	D, log2_ratio_dN_dS = dS_dN(parsed_query, parsed_targets)

	#D_out = open('D.txt', 'w')
	#D_out.write(print(D))
	#D_out.close()

	with open('D.txt', 'w') as D_out:
		for item in D:
			D_out.write("%s\n" % item)

	#log2_out = open('log2_dN_dS.txt', 'w')
	#log2_out.write(print(log2_ratio_dN_dS))
	#log2_out.close()

	with open('log2_dN_dS.txt', 'w') as log2_out:
		for item in log2_ratio_dN_dS:
			log2_out.write("%s\n" % item)
