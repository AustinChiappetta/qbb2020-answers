#!/usr/bin/env python2

import hifive
import h5py
import numpy as np
import matplotlib.pyplot as plt

#fetch the expression data from the bed files
pos_compartment_fpkms = []
neg_compartment_fpkms = []
fs = open('/Users/cmdb/qbb2020-answers/Quant_lab/assignment_09/pos_genes_fpkm.bed', 'r')
count = 0
for line in fs:
	line = line.strip('\n')
	# if count < 10:
	# 	print('line is:', line)
	fields = line.split('\t') #file is tab-separated
	# if count < 10:
	# 	print('fields is:', fields)
	fpkm = float(fields[4])
	pos_compartment_fpkms += [fpkm]
	count += 1
fs.close()

fs = open('/Users/cmdb/qbb2020-answers/Quant_lab/assignment_09/neg_genes_fpkm.bed', 'r')
for line in fs:
	line = line.strip('\n')
	fields = line.split('\t') #file is tab-separated
	fpkm = float(fields[4])
	neg_compartment_fpkms += [fpkm]
fs.close()


#now make violin plots
fig, ax = plt.subplots()
ax.set_title("Expression by compartment")
ax.set_ylabel("FPKM")
ax.set_yscale('log')
labels = ['Pos compartment scores\nB compartment (lower expression)', 'Neg compartment scores\nA compartment (higher expression)']
ax.set_xticks(np.arange(1, len(labels) + 1))
ax.set_xticklabels(labels)
ax.violinplot([pos_compartment_fpkms, neg_compartment_fpkms])
fig.savefig("expression_violin_plots.png")
plt.close(fig)




