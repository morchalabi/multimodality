#!/bin/bash

fls_=$(basename -a ../out/*summits*.bed)

for fl_ in ${fls_}
do
	printf "processing %s" "${fl_}"

	findMotifsGenome.pl ../out/${fl_} \
	/nfsdata/Morteza/projects/CNOT/data/homo_sapien/index_hg38/ref_genome/hg38.fa \
	../out/motif/${fl_} \
	-size 200 \
	-mask
done

