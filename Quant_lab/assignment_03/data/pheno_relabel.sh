#!/bin/bash

#awk '{sub(/\./," ",$1)};1' <file

awk 'NR == 1 {print $0}; NR>1 {gsub(/_/, "\t"); print $0}' /Users/cmdb/qbb2020-answers/Quant_lab/assignment_03/data/BYxRM_PhenoData.txt > PhenoData.relabeled.txt