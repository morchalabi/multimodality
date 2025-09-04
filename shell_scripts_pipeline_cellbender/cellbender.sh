#!/bin/bash

source aws_credentials.sh
source ~/tools/miniconda3/etc/profile.d/conda.sh
conda activate CellBender

s_=(PN1 PN2 PN3 PN4 PN5 PN6 PN7 PN8 PN9 PN10 PN11 PN12 PN13 PN14 PN15 PN16 PN17 PN18 PN19 PN20 PN21 PN22 PN23 PN24 PN25 PN26 PN27 PN28)

for (( i_=0; i_<=28; i_++ ))
do
        f_=${s_[$i_]}
        aws s3 cp s3://pancreas-seq/new_samples/cellranger/v7.0/${f_}/ . --recursive --exclude "*" --include "metrics_summary.csv" --include "raw_feature_bc_matrix.h5"
        # metrics_summary is a csv file delimited by ', but it also uses "" and , for quoting and separating numbers!
        # this block reads 2nd line and extracts first field which is the estimated number of cells
        # it (1) replaces ", with _ then (2) replaces " with nothing, (3) replaces , with nothing, (4) replaces _ with space, finally 1st element is first field
        c_=$(sed -n '2p' "metrics_summary.csv")
        c_=(${c_//\",/_})
        c_=(${c_//\"/})
        c_=(${c_//,/})
        c_=(${c_//_/ })
        c_=${c_[0]}

        echo processing for sample ${f_} with ${c_} cells

        cellbender remove-background \
                 --input raw_feature_bc_matrix.h5 \
                 --output filtered_feature_bc_matrix_Cellbender.h5 \
                 --expected-cells ${c_} \
                 --epochs 150

        aws s3 cp . s3://pancreas-seq/new_samples/cellranger/v7.0/${f_}/ --recursive --exclude "*" --include "*Cellbender*"

        rm raw_feature_bc_matrix.h5 *Cellbender* metrics_summary.csv
done

conda deactivate