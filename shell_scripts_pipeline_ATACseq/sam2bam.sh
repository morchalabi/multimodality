#!/bin/bash

source ~/.profile

# directory of SAM files
sam_dir=$HOME/data/projects/TAD_Enhnacers/data/HiC/GSM2790405

# converting SAM file to BAM file

for sam_fl in alignment_reads_1 alignment_reads_2 
do
	samtools view -Sub -@ 22 ${sam_dir}/${sam_fl}/alg.sam > ${sam_dir}/${sam_fl}/alg.bam
	echo Done for ${sam_fl}
done
