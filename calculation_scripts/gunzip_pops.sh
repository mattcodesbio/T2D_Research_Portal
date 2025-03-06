#!/bin/bash

# Uncompress all the split subpopulation vcf.gz files into vcf format
# SweeD only takes in vcf format 

# List of populations
populations=("BEB" "GIH" "ITU" "PJL" "STU")

# Loop over each population
for pop in "${populations[@]}"; do
    # Construct input and output file names based on population
    input_file="${pop}_chr2.vcf.gz"
    output_file="../subpop_chr2_vcf/${pop}_chr2.vcf"

    # Run gunzip command
    gunzip -c $input_file > $output_file
done
