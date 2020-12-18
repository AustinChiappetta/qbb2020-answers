#!/bin/bash

for FILE in nodes_bin.1.taxa.parsed	nodes_bin.3.taxa.parsed	nodes_bin.5.taxa.parsed	nodes_bin.7.taxa.parsed nodes_bin.2.taxa.parsed	nodes_bin.4.taxa.parsed	nodes_bin.6.taxa.parsed	nodes_bin.8.taxa.parsed
do
	sort -k 1 -h -r ${FILE} > ${FILE}.sorted
done