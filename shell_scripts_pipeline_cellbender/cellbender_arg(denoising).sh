#!/bin/sh

pat=$1
cells=$2
drop=$3

aws s3 sync s3://snrna-seq/cellranger/${pat}/ data/${pat}/ --exclude '*' --include 'raw_feature_bc_matrix.h5' --quiet

sleep 15

cellbender remove-background \
                 --input data/${pat}/raw_feature_bc_matrix.h5 \
                 --output data/${pat}/${pat}.h5 \
                 --expected-cells ${cells} \
                 --total-droplets-included ${drop} \
                 --cuda \
                 --epochs 300

aws s3 sync data/${pat}/ s3://snrna-seq/cellbender/${pat}/ --exclude "*.out" --exclude '*.err' --exclude 'raw_feature*' --exclude '.*' --quiet
