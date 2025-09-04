#!/bin/bash

for cond in siNT si673_4h si673_12h si673_24h
do
	for rep in 1 2 3
	do
		printf "file ../data/%s/MACS2/%s\n" "${cond}" "$rep"

		# sorting

		fl=$(basename -a ../data/${cond}/MACS2/${cond}_${rep}*treat*.bdg)
		#sort -k1,1 -k2,2n ../data/${cond}/MACS2/${fl} > ../data/${cond}/MACS2/sorted_${fl}
		
		fl=$(basename -a ../data/${cond}/MACS2/${cond}_${rep}*control*.bdg)
		sort -k1,1 -k2,2n ../data/${cond}/MACS2/${fl} > ../data/${cond}/MACS2/sorted_${fl}
	done

	# conversion

	fls=$(basename -a ../data/${cond}/MACS2/sorted_*.bdg)
	for fl in ${fls}
	do
        	bedGraphToBigWig ../data/${cond}/MACS2/${fl} \
        	/nfsdata/Morteza/projects/genomes/hg38/hisat2_index/ref/genome.chrom.sizes \
        	../data/${cond}/MACS2/${fl}.bw
	done
done

