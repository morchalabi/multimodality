# multimodality
Showcasing coding for variant calling, python, ATACseq, cellbender and spatial transcriptomics, 

1) **FGA_TMB**: it contains an R script for calculating FGA and TMB in Wilms tumor project. It also contains a bash script pipeline of GATK-based somatic short mutations (SSM) calling that is required for TMB (nonsynonymous mutations). The result of SSM calling can be accessed on GB: https://genome.ucsc.edu/s/mor.chalabi%40gmail.com/TMB_WILMS.
2) **python_scripts**: it contains some python scripts from different projects
	j.ipynb and j.py are from PDAC tropism repository (scripts/Fig_4)
	API.py from EnrichR repository
	mtx2loom.py from veloPipe repository

4) **shell_scripts_pipeline_ATACseq**: it contains bash scripts of ATAC-seq pipeline. These are from the ATACseq repository.

5) **shell_scripts_pipeline_cellbender**: it contains bash scripts of cellbender (denoising snRNAseq data) pipeline

6) **spatial_transcriptomics**: it contains scripts for RCTD-based deconvolution of spatial transcriptome sequencing data in the PDAC tropism project (scripts/Fig_4)
