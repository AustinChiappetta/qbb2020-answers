#!/usr/bin/env python3
import sys
import gzip


E4bedGraph = '/Users/cmdb/qbb2020-answers/Quant_lab/assignment_06/data/E4.0ICM_rep1_1_STEMseq.chr6_bismark_bt2_pe.bedGraph.gz'
E5bedGraph = '/Users/cmdb/qbb2020-answers/Quant_lab/assignment_06/data/E5.5Epi_rep1_1_STEMseq.ch6_bismark_bt2_pe.bedGraph.gz'
#genesbed = '/Users/cmdb/qbb2020-answers/Quant_lab/assignment_06/data/mm10_refseq_genes_chr6_50M_60M.bed'
genesbed = '/Users/cmdb/qbb2020-answers/Quant_lab/assignment_06/data/mm10_refseq_genes_chr6_50M_60M.sorted.bed'


#For each gene, find the fold change in mean methylation signal from E4 to E5.5 cells. If the mean methylation for a gene in E4 is zero, skip it.

#iterate through coordinates file and make list: [[gene start, gene stop, gene name]]
#13th column of coordinates file has the unique genes names, use the first entry for a given gene name
#5th and 6th columns of the coordinates file contain the coordinates
def gene_coords_parser(bedfile):
	genelist = []
	genedata = []
	bed = open(bedfile, 'r')
	for line in bed: #iterate through binary file
		fields = line.split()
		gene_name = fields[12]
		start_coord = int(fields[4])
		end_coord = int(fields[5])
		if gene_name in genelist:
			continue #skip duplicate genes in bed file
		else:
			genelist += [gene_name]
			genedata += [[start_coord, end_coord, gene_name]]

	bed.close()
	return(genedata)


#bismark bedgraph format: <chromosome> <start position> <end position> <methylation percentage>
#Each methylation site is 2 bp
#Use gzip to open the E4bedgraph file, then iterate through it
#make a new dictionary: E4data = {gene1: {site1: (methylation_percentage),
										#site2: (methylation_percentage)},
								#gene2: {site1: (methylation_percentage),
										#site2: (methylation_percentage)}
										#site3: (methylation_percentage)}


def methylation_bedgraph_parser(bedgraph, parsed_gene_data):
	#returns a dictionary: {gene1: methyl_signal, gene2: methyl_signal, etc.}
	methyl_dict = {}
	genedata = parsed_gene_data
	# print('gene list is this long:', len(genedata))
	# for entry in genedata:
	# 	print(entry)

	with gzip.open(bedgraph, 'r') as fs:
		count = 0
		for line in fs:
			gene_list_position = 0
			count += 1
			if b'track' in line: #skip header
				continue
			line = line.strip(b'\n')
			fields = line.split(b'\t')
			start = int(fields[1])
			end = int(fields[2])
			methylation_percentage = float(fields[3])

			if start < genedata[0][0]: #skip site if it is before the first gene
				continue
			if start > genedata[-1][0]: #skip site if it is after the last gene
				continue
	
			# itercount = 0
			while True:
				# itercount += 1
				#if itercount < 10 and count < 10:
				#print('fields is:', fields)
				# print('start is:', start)
				# print('end is:', end)
				# print('methyl_percentage is:', methylation_percentage)
				# print('genedata entry pointing to:', genedata[gene_list_position])
				#print()

				if start > genedata[gene_list_position][1] and end < genedata[gene_list_position + 1][0]: #if this site is between two genes
					# print('Between genes, skipping this site. Start:', start, 'End:', end)
					# print()
					break

				if start >= genedata[gene_list_position][0] and end <= genedata[gene_list_position][1]: #this site is in gene being pointed to
					gene = genedata[gene_list_position][2]
					methyl_dict.setdefault(gene, {})
					methyl_dict[gene].setdefault(start, 0)
					methyl_dict[gene][start] = methylation_percentage
					# print('site is in this gene, data stored. moving on')
					# print()
					break

				elif start == genedata[gene_list_position][1] or end == genedata[gene_list_position][0]: #site overlaps with endge of gene
					gene = genedata[gene_list_position][2]
					methyl_dict.setdefault(gene, {})
					methyl_dict[gene].setdefault(start, 0)
					methyl_dict[gene][start] = methylation_percentage
					# print('site is in this gene, data stored. moving on')
					# print()
					break

				else:
					gene_list_position += 1 #site is not in this gene, move to next position in gene data list
					# print('site not in gene, trying again...')

	# for entry in methyl_dict.keys():
	# 	print(entry, methyl_dict[entry])

	fs.close()


	#calculate methylation signal across each gene, store in a dict, i.e., {gene1: methyl_signal, gene2: methyl_signal}
	#methylation signal = num of bases * percent methylation. Add up methylation signal across each gene (each methylation site is a CpG, so 2 bp long)
	genes_and_summed_signals = {}
	for gene in methyl_dict.keys():
		methyl_signal = 0
		for site in methyl_dict[gene].keys():
			signal = methyl_dict[gene][site]
			methyl_signal += (2*signal)/100

		genes_and_summed_signals.setdefault(gene, 0)
		genes_and_summed_signals[gene] = methyl_signal

	# for entry in genes_and_summed_signals.keys():
	# 	print(entry, genes_and_summed_signals[entry])

	return(genes_and_summed_signals)


#do the same for E5

#Calculate the fold change in signal from E4 to E5.5 (if mean methylation for a gene in E4 is zero, skip it)
#Fold change = (Y - X)/X
#Add the gene name and relative enrichment to a list
#Output the list of gene names and relative methylation enrichments to a new file

def signal_fold_change(data1, data2):
	#returns a list of tuples: [(entry1, fold_change), (entry2, fold_change), etc.]
	dict1 = data1
	dict2 = data2
	enrichments = []
	for entry in dict1.keys():
		signal1 = dict1[entry]
		if entry in dict2:
			signal2 = dict2[entry]
		else:
			signal2 = 0 #entry not in d2, call as zero

		if signal1 == 0:
			continue

		fold_change = (signal2 - signal1)/signal1
		enrichments += [(entry, fold_change)]

	# for tup in enrichments:
	# 	print(tup)

	return(enrichments)


def write_fold_change_to_file(data_list):
	fs = open('E4toE5.5_Methylation_Enrichment.txt', 'w+')
	fs.write("Gene_name" + "\t" + "Relative_Methylation_Enrichment" + "\n")
	for tup in data_list:
		gene = tup[0]
		enrichment = tup[1]
		fs.write(str(gene) + "\t" + str(enrichment) + "\n")

	fs.close()



#with gzip.open(E4bedGraph, 'r') as E4:
#	for line in E4:
#		if b'track' in line:
#			continue
#		#if count > 10:
#			#break
#		line = line.strip(b'\n')
#		fields = line.split(b'\t')
#		if float(fields[2]) - float(fields[1]) == 0:
			#print(fields)
#			zero_count += 1
#		elif float(fields[2]) - float(fields[1]) == 1:
#			one_count += 1
#		elif float(fields[2]) - float(fields[1]) > 1:
#			two_or_greater_count += 1
#		line_count += 1

#print('zero count is:', zero_count)
#print('one count is:', one_count)
#print('two or greater count is:', two_or_greater_count)
#print('line count is:', line_count)


if __name__ == '__main__':
	genedata = gene_coords_parser(genesbed)
	E4methyldict = methylation_bedgraph_parser(E4bedGraph, genedata)
	E5methyldict = methylation_bedgraph_parser(E5bedGraph, genedata)
	fold_change_data = signal_fold_change(E4methyldict, E5methyldict)
	write_fold_change_to_file(fold_change_data)