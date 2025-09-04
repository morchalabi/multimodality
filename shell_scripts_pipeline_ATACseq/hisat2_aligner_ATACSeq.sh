#!/bin/bash

HISAT2_INDEXES=/nfsdata/Morteza/projects/genomes/hg38/hisat2_index/

dir=../data/fastq/

for f in A01 A02 A03 A04
do
	echo processing ${f}

	hisat2 \
	-x ${HISAT2_INDEXES}/genome_snp_tran \
	-1 ${dir}/${f}_R1_001_cutadapt.fastq.gz \
	-2 ${dir}/${f}_R2_001_cutadapt.fastq.gz \
	-S ${dir}/${f}_001_cutadapt.sam \
	--no-spliced-alignment \
	-k 12 \
	--secondary \
	--summary-file ${dir}/${f}_001_cutadapt_summary.txt \
	--met-file ${dir}/${f}_001_cutadapt_metrics.txt \
	--reorder \
	-p 20

	# sam to bam then bam is coordinate-sorted: -h: include header; -b: output is bam
	samtools view -h -b ${dir}/${f}_001_cutadapt.sam > ${dir}/${f}_001_cutadapt.bam
	samtools sort ${dir}/${f}_001_cutadapt.bam -o ${dir}/${f}_001_cutadapt_sorted.bam

	# removing MT (make sure it is not chrMT in bam) reads, multimappers, unaligned reads, and duplicates
	sambamba view -h -t 10 -f bam \
	-F "ref_name != 'MT' and mapping_quality > 10 and not unmapped and not duplicate" \
	${dir}/${f}_001_cutadapt_sorted.bam > ${dir}/${f}_001_cutadapt_sorted_filtered.bam

	# indexing bam to generte bai (for UCSC GB)
	samtools index ${dir}/${f}_001_cutadapt_sorted_filtered.bam
done
