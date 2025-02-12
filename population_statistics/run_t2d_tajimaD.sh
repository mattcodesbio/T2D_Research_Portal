#!/bin/bash

# Define the populations and gene regions
populations=("BEB" "GIH" "ITU" "PJL" "STU")
gene_regions=("ADCY5" "C3orf65" "CDC123" "CDKAL1" "CDKN1C" "CDKN2B" "COBLL1" "HHEX" "HNF4A" "IGF2BP2" "SLC30A8" "TCF7L2")

# Window size for Tajima’s D
WINDOW_SIZE=10000

# Loop through each population
for pop in "${populations[@]}"; do
    echo "🔍 Processing population: $pop"
    
    # Define population directory
    pop_dir="./$pop"
    
    # Create an output directory for Tajima's D results
    mkdir -p "$pop_dir/${pop}_t2d_tajimaD"

    # Loop through each gene region
    for gene in "${gene_regions[@]}"; do
        echo "Processing gene region: $gene"
        vcf_file="$pop_dir/${pop}_t2dvcfs/${pop}_${gene}_t2dsnp.vcf.gz"

        if [[ -f "$vcf_file" ]]; then
            echo "Calculating Tajima's D for $pop in $gene..."
            vcftools --gzvcf "$vcf_file" --TajimaD $WINDOW_SIZE --out "$pop_dir/${pop}_t2d_tajimaD/${pop}_${gene}_t2d"
        else
            echo "Skipping $pop in $gene (VCF file not found: $vcf_file)"
        fi
    done
done

echo "Tajima's D calculations completed for all populations!"


