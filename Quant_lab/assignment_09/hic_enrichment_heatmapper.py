#!/usr/bin/env python2

import hifive
import h5py
import numpy as np
import matplotlib.pyplot as plt

hic = hifive.HiC('/Users/cmdb/qbb2020-answers/Quant_lab/assignment_09/data/hic_project', 'r')
data = hic.cis_heatmap('chr13', 1000000, datatype = 'fend', arraytype = 'full', diagonalincluded = True)

#add a pseudocount to the raw and expected values to resolve the divide by zero issue
enrichment = data
enrichment[:, :, 0] += 1
enrichment[:, :, 1] += 1 
enrichment = enrichment[:, :, 0] / enrichment[:, :, 1] #enrichment is now a 2D array
enrichment = np.log10(enrichment) #take a log of the values for plotting

#print(enrichment.shape)
fig, ax = plt.subplots()
ax.set_title("Chr13 Enrichment Heatmap (log10-transformed)")
ax.set_xlabel("Chr13 nucleotide pos")
ax.set_ylabel("Chr13 nucleotide pos")
c = ax.pcolor(enrichment)
fig.colorbar(c, ax=ax)
fig.savefig("Chr13_Enrichment_Heatmap.png")
plt.close(fig)



