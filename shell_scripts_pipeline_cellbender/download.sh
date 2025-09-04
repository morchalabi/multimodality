#!/bin/bash

source ./aws_credentials.sh

for i_ in W040 W132 W241 W340 W495 W612 W050 W180 W184 W370 W441 W944
do
	echo ${i_}
	aws s3 cp s3://wilms-tumor/cellranger/v6.1.2/GEX/${i_}/ /Users/morchalabi/downloads --recursive --exclude "*" --include "filtered_feature_bc_matrix_Cellbender_filtered.h5" --include "filtered_feature_bc_matrix.h5"
	mkdir ~/downloads/${i_}
	mv ~/downloads/*.h5 ~/downloads/${i_}
done
