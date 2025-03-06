#!/bin/bash

# Processes the intergrated variant call sets from the 1000 Genomes Project 
# to be separated into the populations of interest based on the sample IDs
# Sample IDs for each population is downlaoded and stored in a {subpop}_id txt file (e.g. BEB_id.txt)
# Outputs a vcf.gz file that contains only the subpopulations for each chromosome 


# List of populations
populations=("BEB" "GIH" "ITU" "PJL" "STU")

# Loop over each population
for pop in "${populations[@]}"; do
    # Construct input and output file names based on population
    input_file="all_chr6.vcf.gz"
    output_file="${pop}_chr6.vcf.gz"
    subpop_file="../subpop_id/${pop}_id.txt"

    # Run bcftools command
    bcftools view -S $subpop_file --force-samples -Oz -o $output_file --threads 18 $input_file
done
