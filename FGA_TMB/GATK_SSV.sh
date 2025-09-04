#!/bin/bash

#To generate common biallelic SNPs (ref allele + 1 alternative allele):
#bcftools view \
#  -m2 -M2 \                     # Keep only biallelic sites
#  -v snps \                     # Keep only SNPs
#  -i 'AF[0] >= 0.01' \          # Filter by allele frequency (AF ≥ 0.01)
#  af-only-gnomad.hg38.vcf.gz -Oz -o common_biallelic.hg38.vcf.gz

# input files
REFERENCE=/data/tools/reference_genome/USCS/hg38.fa
TUMOR_BAM=/data/Wilms/data/W612/W612.bam
GERMLINE_RESOURCE=/data/tools/reference_genome/broad_gnomAD/af-only-gnomad.hg38.vcf.gz
COMMON_SNPS=/data/tools/reference_genome/broad_gnomAD/common_biallelic.hg38.vcf.gz
FUNCOTATOR_DS=/data/tools/reference_genome/Funcotator_dataSources/

# output prefix and directory
OUTDIR=/data/Wilms/out/W612/SSV
OUTPUT_PREFIX=${OUTDIR}/W612

# Step 1: Mutect2 in tumor-only mode
gatk Mutect2 \
  -R $REFERENCE \
  -I $TUMOR_BAM \
  --germline-resource $GERMLINE_RESOURCE \
  -O ${OUTPUT_PREFIX}.unfiltered.vcf.gz

# Step 2: GetPileupSummaries (using common SNPs resource)
gatk GetPileupSummaries \
  -R $REFERENCE \
  -I $TUMOR_BAM \
  -V $COMMON_SNPS \
  -L $COMMON_SNPS \
  -O ${OUTPUT_PREFIX}.pileups.table

# Step 3: CalculateContamination
gatk CalculateContamination \
  -I ${OUTPUT_PREFIX}.pileups.table \
  -O ${OUTPUT_PREFIX}.contamination.table

# Step 4: FilterMutectCalls
gatk FilterMutectCalls \
  -R $REFERENCE \
  -V ${OUTPUT_PREFIX}.unfiltered.vcf.gz \
  --contamination-table ${OUTPUT_PREFIX}.contamination.table \
  --min-allele-fraction 0.01 \
  --unique-alt-read-count 3 \
  --min-median-mapping-quality 30 \
  --max-alt-allele-count 1 \
  -O ${OUTPUT_PREFIX}.filtered.all.vcf.gz

# Step 5: keep only PASS and weak alleles
bcftools view \
  -i 'FILTER="PASS" || FILTER="weak_evidence"' \
  ${OUTPUT_PREFIX}.filtered.all.vcf.gz -Oz -o ${OUTPUT_PREFIX}.filtered.vcf.gz

# Step 6: Funcotator annotation
gatk Funcotator \
  -R $REFERENCE \
  -V ${OUTPUT_PREFIX}.filtered.vcf.gz \
  --ref-version hg38 \
  --data-sources-path $FUNCOTATOR_DS \
  --output ${OUTPUT_PREFIX}_filtered_PASS-weak-evidence_annotated.vcf.gz \
  --output-file-format VCF
