#!/bin/bash

# N.B. source from /data/runs/
# put multi_config.csv in /data/
# put fasq files in /data/fsqs/
# put refrence geneome and VDJ reference directories in /data/
# put aws_credentials.sh in /data/
# set parameters accordingly in multi_config.csv

source /data/aws_credentials.sh

for i_ in SC_TCR_10X SC_5PGEX_10X
do
        echo ${i_}

        mkdir ../fasqs/${i_}
        aws s3 cp s3://cpoi-pal/fastq/${i_}/ fasqs/${i_}/ --recursive
done

cellranger multi --id=Human_TCell_GEX_VDJ --csv=../multi_config.csv

cd Human_TCell_GEX_VDJ
aws s3 cp ./ s3://cpoi-pal/cellranger/v6.1.2/ --recursive
