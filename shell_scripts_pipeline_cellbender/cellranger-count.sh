#!/bin/bash

source aws_credentials.sh

for i_ in W001 W002 W003 W004 W005 W006 W007 W008 W009 W010 W011 W012
do
	echo ${i_}

	mkdir out/${i_}
	aws s3 cp s3://wilms-tumor/raw_data/GEX/${i_}/ out/${i_} --recursive

	cellranger count --id=${i_} \
	--fastqs=/data/out/${i_} \
	--transcriptome=/data/refdata-gex-GRCh38-2020-A \
	--include-introns \
	--expect-cells=1000 \
	--localcores=14

	aws s3 cp ${i_}/outs/ s3://wilms-tumor/cellranger/v6.1.2/GEX/${i_}/ --recursive

	rm -rf ${i_} out/${i_}

done
