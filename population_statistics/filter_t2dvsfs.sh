#!/bin/bash

# Define populations and gene regions
populations=("BEB" "GIH" "ITU" "PJL" "STU")
gene_regions=("ADCY5" "C3orf65" "CDC123" "CDKAL1" "CDKN1C" "CDKN2B" "COBLL1" "HHEX" "HNF4A" "IGF2BP2" "SLC30A8" "TCF7L2")

# Path to the T2D-associated SNPs positions file
t2d_snps_file="t2dsnps_positions.txt"

# Loop through each population and gene
for pop in "${populations[@]}"; do
    for gene in "${gene_regions[@]}"; do
        # Define input and output VCF file paths
        input_vcf="./$pop/${pop}_${gene}.vcf.gz"
        output_dir="./$pop/${pop}_t2dvcfs"
        output_vcf="$output_dir/${pop}_${gene}_t2dsnp.vcf.gz"

        # Check if the input VCF file exists
        if [[ -f "$input_vcf" ]]; then
            echo "Filtering T2D SNPs for $pop in $gene..."

            mkdir -p "$output_dir"

            # Use bcftools to filter the VCF
            bcftools view -T "$t2d_snps_file" -Oz -o "$output_vcf" "$input_vcf"

            echo "Filtered VCF saved: $output_vcf"
        else
            echo "Skipping $pop in $gene (VCF file not found: $input_vcf)"
        fi
    done
done

echo "T2D SNP filtering complete!"