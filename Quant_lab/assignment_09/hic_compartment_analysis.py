#!/usr/bin/env python2

import hifive
import h5py
import numpy as np
import matplotlib.pyplot as plt

hic = hifive.HiC('/Users/cmdb/qbb2020-answers/Quant_lab/assignment_09/data/hic_project', 'r')
data = hic.cis_heatmap('chr13', 1000000, datatype = 'fend', arraytype = 'full', diagonalincluded = True)

Comp = hifive.hic_domains.Compartment(hic, 100000, chroms=['chr13'], out_fname='tmp.hdf5')
Comp.write_eigen_scores('hic_comp.bed')

X = Comp.positions['chr13']
Y = Comp.eigenv['chr13']

sub = X[:, 0]


# print(X)
# print(sub)
#print(Y)
# print(Y.shape)
# print(sub.shape)

fig, ax = plt.subplots()
ax.set_title("Compartment Analysis of Chr13")
ax.set_xlabel("Genomic positions")
ax.set_ylabel("Compartment scores")
ax.scatter(sub, Y)
#fig.show()
fig.savefig("compartment_scores.png")
plt.close(fig)



