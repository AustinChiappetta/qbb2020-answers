#!/usr/bin/env python2

import hifive
import h5py
import numpy as np
import matplotlib.pyplot as plt
import pyBigWig

bw = pyBigWig.open('data/WT_H3K27me3.bw')
#bw.stats('chr1', 1100000, 1200000, type='sum')

A_expression = []
A_H3K27me3 = []
B_expression = []
B_H3K27me3 = []
#read through pos_genes_fpkm.bed. 
#1) Grab chromosome label (field 0), start pos (field 1), end pos (field 2), 
#gene name (field 3), and expression (field 4)
#2) enter values into bw.stats and store output (H3K27me3 level). If value is None, replace with 0
#make data list: [[gene name, expression, H3K27me3 level]]
#do the same for neg_genes_fpkm.bed
fs = open('/Users/cmdb/qbb2020-answers/Quant_lab/assignment_09/pos_genes_fpkm.bed', 'r')
for line in fs:
	line = line.strip('\n')
	fields = line.split('\t')
	chr = fields[0]
	start = int(fields[1])
	end = int(fields[2])
	gene_name = fields[3]
	fpkm = float(fields[4])
	H3K27me3 = bw.stats(chr, start, end, type = 'sum')
	if H3K27me3 == None:
		H3K27me3 = 0
	B_expression += [fpkm]
	B_H3K27me3 += [H3K27me3]
fs.close()

fs = open('/Users/cmdb/qbb2020-answers/Quant_lab/assignment_09/neg_genes_fpkm.bed', 'r')
for line in fs:
	line = line.strip('\n')
	fields = line.split('\t')
	chr = fields[0]
	start = int(fields[1])
	end = int(fields[2])
	gene_name = fields[3]
	fpkm = float(fields[4])
	H3K27me3 = bw.stats(chr, start, end, type = 'sum')
	if H3K27me3 == None:
		H3K27me3 = 0
	A_expression += [fpkm]
	A_H3K27me3 += [H3K27me3]
fs.close()



#Make 2 scatter plots
fig, axes = plt.subplots(ncols = 2, sharey = True, sharex = True)
# fig.suptitle('Relationship between expresssion and H3K27me3')
axes[0].set_yscale('log')
# axes[0].set_xscale('log')
axes[0].scatter(A_expression, A_H3K27me3)
axes[0].set_xlabel('Expression (FPKM)')
axes[0].set_ylabel('H3K27me3 signal')
axes[0].set_title('A compartment expression \nvs repression')

axes[1].set_yscale('log')
# axes[1].set_xscale('log')
axes[1].scatter(B_expression, B_H3K27me3)
axes[1].set_xlabel('Expression (FPKM)')
axes[1].set_ylabel('H3K27me3 signal')
axes[1].set_title('B compartment expression \nvs repression')

fig.tight_layout()
fig.savefig("expression_H3K27me3.png")
plt.close(fig)





