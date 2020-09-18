#!/bin/bash

for SAM in A01_09.fastq.sam A01_24.fastq.sam A01_35.fastq.sam A01_63.fastq.sam A01_11.fastq.sam A01_27.fastq.sam A01_39.fastq.sam A01_23.fastq.sam A01_31.fastq.sam A01_62.fastq.sam
do
	samtools sort -O bam -o ${SAM}.bam ${SAM}
done