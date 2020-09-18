#!/bin/bash

for SAMPLE in A01_09.fastq A01_23.fastq A01_27.fastq A01_35.fastq A01_62.fastq A01_11.fastq A01_24.fastq A01_31.fastq A01_39.fastq A01_63.fastq
do 
	bwa mem -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}" -o ${SAMPLE}.sam /Users/cmdb/qbb2020-answers/Quant_lab/assignment_02/ref_yeast/sacCer3.fa /Users/cmdb/qbb2020-answers/Quant_lab/assignment_02/reads/${SAMPLE}
done
