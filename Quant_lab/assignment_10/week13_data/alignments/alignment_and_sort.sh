#!/bin/bash

for SAMPLE in SRR492183 SRR492188 SRR492190 SRR492194 SRR492186 SRR492189 SRR492193 SRR492197
do 
	#bwa mem -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}" -o ${SAMPLE}.sam /Users/cmdb/qbb2020-answers/Quant_lab/assignment_10/week13_data/assembly.fasta /Users/cmdb/qbb2020-answers/Quant_lab/assignment_10/week13_data/READS/${SAMPLE}_1.fastq /Users/cmdb/qbb2020-answers/Quant_lab/assignment_10/week13_data/READS/${SAMPLE}_2.fastq
	bwa mem -o ${SAMPLE}.sam /Users/cmdb/qbb2020-answers/Quant_lab/assignment_10/week13_data/assembly.fasta /Users/cmdb/qbb2020-answers/Quant_lab/assignment_10/week13_data/READS/${SAMPLE}_1.fastq /Users/cmdb/qbb2020-answers/Quant_lab/assignment_10/week13_data/READS/${SAMPLE}_2.fastq
	samtools sort -O bam -o ${SAMPLE}.sorted.bam /Users/cmdb/qbb2020-answers/Quant_lab/assignment_10/week13_data/alignments/${SAMPLE}.sam
done
