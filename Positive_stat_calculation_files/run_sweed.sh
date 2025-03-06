#!/bin/bash

# Script to run SweeD-P on the 1000 Genomes Project vcf files for the 5 subpopulations and 8 chromosomes of interest
# Grid size of 100000 (i.e. SweeD scans the chromosome 100000 times) -> Fine resolution across the chromosome 
# Adjust threads according to the system specs. 

# List of populations
populations=("BEB" "GIH" "ITU" "PJL" "STU")

# List of chromosomes to process
chromosomes=("2" "3" "6" "8" "9" "10" "11" "20")

# Loop over each population
for pop in "${populations[@]}"; do
    # Loop over each chromosome
    for chr in "${chromosomes[@]}"; do
        # Construct the input file path based on population and chromosome
        input_file="/path/to/vcf_data/subpop_chr${chr}_vcf/${pop}_chr${chr}.vcf"
        
        # Construct the name for the output file
        output_name="${pop}_chr${chr}_100000"
        
        # Run the SweeD-P command for the current population and chromosome
        ./SweeD-P -name $output_name -input $input_file -grid 100000 -threads 18
    done
done
