#!/usr/bin/env python3
import sys

def kraken_parser(kraken_file):
	"""parses .kraken files and generates an output file for ktImportText.
	Input: file path to .kraken (can take multiple space-separated file paths)
	Output: .parsed text file(s). First column is quantity. Remaining columns are tab-separated taxonomic wedges
	"""

	#iterate through kraken file
	#make dictionary: keys are tax IDs, values are number of times that tax ID occurs in the file
	#then, make parsed file (write to new file: filename.parsed). Should be tab-separated
	#column 1: quantity/occurnces of this tax ID
	#columns 2 and beyond: taxonomic ID (this field is semicolon-separated in the .kraken files, need to make tab-separated)

	taxIDdict = {}

	kraken = open(kraken_file, 'r')
	for line in kraken:
		line = line.strip('\n')
		fields = line.split('\t')
		taxid = fields[1]
		taxIDdict.setdefault(taxid, 0)
		taxIDdict[taxid] += 1
	kraken.close()

	paths = kraken_file.split('/')
	output_name = paths[-1]
	output_name = output_name.strip('.kraken') + ".parsed" #default behavior is to set output name as filename.parsed

	parsed = open(output_name, 'w')
	for key in taxIDdict.keys():
		quantity = str(taxIDdict[key])
		parsed.write(quantity)
		ID = key.split(';')
		for wedge in ID:
			parsed.write('\t' + wedge)
		parsed.write('\n')
	parsed.close()

if __name__ == "__main__":
	files = sys.argv[1:]
	# print(files)
	for file in files:
		kraken_parser(file)





