#!/bin/bash

for BIN in bin.1 bin.3 bin.5 bin.7 bin.2 bin.4 bin.6 bin.8
do
	grep "NODE" ${BIN}.fa > nodes_${BIN}
done