#!/bin/bash

for NODEFILE in nodes_bin.1	nodes_bin.3	nodes_bin.5	nodes_bin.7 nodes_bin.2	nodes_bin.4	nodes_bin.6	nodes_bin.8
do
	while IFS= read -r line; do
		grep ${line:1} /Users/cmdb/qbb2020-answers/Quant_lab/assignment_10/week13_data/KRAKEN/assembly.kraken >> ${NODEFILE}.taxa
	done < $NODEFILE
done		