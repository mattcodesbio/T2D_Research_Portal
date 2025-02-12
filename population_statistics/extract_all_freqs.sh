#!/bin/bash

# Define populations and gene regions
populations=("BEB" "GIH" "ITU" "PJL" "STU")
gene_regions=("ADCY5" "C3orf65" "CDC123" "CDKAL1" "CDKN1C" "CDKN2B" "COBLL1" "HHEX" "HNF4A" "IGF2BP2" "SLC30A8" "TCF7L2")

# Loop through each population and gene
for pop in "${populations[@]}"; do
    for gene in "${gene_regions[@]}"; do
        # Define VCF file path
        vcf_file="./$pop/${pop}_${gene}.vcf.gz"
        output_dir="./$pop/${pop}_freqs"
        output_file="$output_dir/${pop}_${gene}_freq.txt"

        # Check if VCF file exists
        if [[ -f "$vcf_file" ]]; then
            echo "Extracting allele frequency for $pop in $gene..."
            
            # Create output directory if it doesn't exist
            mkdir -p "$output_dir"

            # Run bcftools query to extract AF
            bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%AF\n' "$vcf_file" > "$output_file"
            
            echo "Saved: $output_file"
        else
            echo "Skipping $pop in $gene (VCF file not found: $vcf_file)"
        fi
    done
done

