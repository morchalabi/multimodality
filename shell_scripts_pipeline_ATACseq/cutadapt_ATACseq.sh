#!/bin/bash

dir=../data/fastq

for f in A01 A02 A03 A04
do
	echo processing ${f}

	cutadapt \
	-j 32 \
	-m 35 \
	-a "CTGTCTCTTATACACATCT;max_error_rate=0.2" \
	-A "CTGTCTCTTATACACATCT;max_error_rate=0.2" \
	-O 5 \
	-o ${dir}/${f}_R1_001_cutadapt.fastq.gz \
	-p ${dir}/${f}_R2_001_cutadapt.fastq.gz \
	${dir}/${f}_R1_001.fastq.gz \
	${dir}/${f}_R2_001.fastq.gz > ${dir}/report_cutadapt.txt
done

