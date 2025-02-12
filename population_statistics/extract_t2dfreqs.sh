#!/bin/bash

# Define populations and gene regions
populations=("BEB" "GIH" "ITU" "PJL" "STU")
gene_regions=("ADCY5" "C3orf65" "CDC123" "CDKAL1" "CDKN1C" "CDKN2B" "COBLL1" "HHEX" "HNF4A" "IGF2BP2" "SLC30A8" "TCF7L2")

# Path to the T2D-associated SNPs file
t2d_snps_file="t2dsnps_positions.txt"

# Loop through each population and gene
for pop in "${populations[@]}"; do
    for gene in "${gene_regions[@]}"; do
        # Define input and output file paths
        input_freq_file="./$pop/${pop}_freqs/${pop}_${gene}_freq.txt"
        output_dir="./$pop/${pop}_t2dfreqs"
        output_freq_file="$output_dir/${pop}_${gene}_t2dsnp_freq.txt"

        # Check if the input frequency file exists
        if [[ -f "$input_freq_file" ]]; then
            echo "Extracting T2D SNP frequencies for $pop in $gene..."

            mkdir -p "$output_dir"

            # Use awk to extract matching positions
            awk 'NR==FNR {snp[$1"\t"$2]; next} ($1"\t"$2) in snp' "$t2d_snps_file" "$input_freq_file" > "$output_freq_file"

            echo "Extracted SNP frequencies saved: $output_freq_file"
        else
            echo "Skipping $pop in $gene (Frequency file not found: $input_freq_file)"
        fi
    done
done

echo "T2D SNP frequency extraction complete!"
