#!/bin/bash

## Download vcf files from the 1000 Genomes Project
## Download vcf index files from the 1000 Genomes Project

# Directory to save vcf_files
raw_vcf="/path/to/raw_vcf"
mkdir -p "$raw_vcf"

# Chromosomes of interest
chromosome=("2" "3" "6" "8" "9" "10" "11" "20")


# URL to FTP server
ftp_url="ftp://ftp.ensembl.org/pub/data_vcf_files/homo_sapiens/GRCh38/variation_genotype"

# Loop through the chromosomes
for chr in "${chromosome[@]}"; do
    vcf_file="ALL.chr${chr}_GRCh38.genotypes.20170504.vcf.gz"
    vcf_index="${vcf_file}.tbi"

    # Download VCF file
    echo "Downloading ${vcf_file}..."
    curl -o "${raw_vcf}/${vcf_file}" "${ftp_url}/${vcf_file}" || echo "Download failed for ${vcf_file}" 

    # Download VCF index file
    echo "Downloading ${vcf_index}..."
    curl  -o "${raw_vcf}/${vcf_index}" "${ftp_url}/${vcf_index}" || echo "Download failed for ${vcf_index}" 

done
