#!/usr/bin/env python 3

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

#see this link for code: https://matplotlib.org/3.3.2/api/_as_gen/matplotlib.pyplot.bar.html

def CTCF_plots():
	#extract the data
	gained_file = open('/Users/cmdb/qbb2020-answers/Quant_lab/assignment_05/g1e/data/CTCF_sites_gained_count.txt')
	for data in gained_file:
		#print(data)
		field = data.strip()
		#print(field)
		CTCF_gained = int(field)
	gained_file.close()
	#print('gained data is', CTCF_gained)

	lost_file = open('/Users/cmdb/qbb2020-answers/Quant_lab/assignment_05/g1e/data/CTCF_sites_lost_count.txt')
	for data in lost_file:
		field = data.strip()
		CTCF_lost = int(field)
	lost_file.close()
	#print('lost data is', CTCF_lost)

	ER4_features_file = open('/Users/cmdb/qbb2020-answers/Quant_lab/assignment_05/g1e/data/ER4_feature_counts.txt')
	for line in ER4_features_file:
		#print('line is', line)
		fields= line.strip()
		fields = fields.split(' ')
		#print('fields is', fields)
		if 'exon' in fields:
			ER4_exons = int(fields[0])
		elif 'intron' in fields:
			ER4_introns = int(fields[0])
		elif 'promoter' in fields:
			ER4_promoters = int(fields[0])
	ER4_features_file.close()
	#print('ER4 exons data is', ER4_exons)
	#print('ER4 introns data is', ER4_introns)
	#print('ER4 promoters data is', ER4_promoters)

	G1E_features_file = open('/Users/cmdb/qbb2020-answers/Quant_lab/assignment_05/g1e/data/G1E_feature_counts.txt')
	for line in G1E_features_file:
		fields = line.strip()
		fields = fields.split(' ')
		if 'exon' in fields:
			G1E_exons = int(fields[0])
		elif 'intron' in fields:
			G1E_introns = int(fields[0])
		elif 'promoter' in fields:
			G1E_promoters = int(fields[0])
	G1E_features_file.close()
	#print('G1E exons data is', G1E_exons)
	#print('G1E introns data is', G1E_introns)
	#print('G1E promoters data is', G1E_promoters)



	#two subplots (1 row, 0 columns)
	#left plot is number of CTCF sites grouped by cell type and feature type 
		#see this link for code: https://matplotlib.org/3.3.1/gallery/lines_bars_and_markers/barchart.html
	#right plot is CTCF sites lost/gained upon differentiation of G1E to ER4
	
	feature_labels = ['Promoter', 'Exon', 'Intron']
	G1E_features = [G1E_promoters, G1E_exons, G1E_introns]
	ER4_features = [ER4_promoters, ER4_exons, ER4_introns]
	x = np.arange(len(feature_labels))  # the label locations
	width = 0.35
	fig, ax = plt.subplots(ncols = 2, figsize=(10, 5))
	rects1 = ax[0].bar(x - width/2, G1E_features, width, label='G1E')
	rects2 = ax[0].bar(x + width/2, ER4_features, width, label='ER4')

	# Add some text for labels, title and custom x-axis tick labels, etc.
	ax[0].set_ylabel('# of CTCF binding sites')
	ax[0].set_title('CTCF binding sites in each region for G1E and ER4')
	ax[0].set_xticks(x)
	ax[0].set_xticklabels(feature_labels)
	ax[0].legend()

	#lost = ax[1].bar(x - width, CTCF_lost, width, label='CTCF sites lost')
	#gained = ax[1].bar(x + width, CTCF_gained, width, label='CTCF sites gained')
	#ax[1].bar(['lost', 'gained'], [CTCF_lost, CTCF_gained])
	lost = ax[1].bar(['lost'], CTCF_lost, 0.5, label='CTCF sites lost')
	gained = ax[1].bar(['gained'], CTCF_gained, 0.5, label='CTCF sites gained')
	ax[1].set_ylabel('# of CTCF binding sites')
	ax[1].set_title('CTCF binding sites lost/gained upon differentiation of G1E to ER4')
	#ax[1].set_xticks(x)
	#ax[1].set_xticklabels(feature_labels)
	ax[1].legend()

	fig.tight_layout()
	plt.show()
	plt.savefig('G1E_ER4.png')


if __name__ == "__main__":
	CTCF_plots()


