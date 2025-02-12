#!/bin/bash

# Define populations and genes
populations=("BEB" "GIH" "ITU" "PJL" "STU")
gene_regions=("ADCY5" "C3orf65" "CDC123" "CDKAL1" "CDKN1C" "CDKN2B" "COBLL1" "HHEX" "HNF4A" "IGF2BP2" "SLC30A8" "TCF7L2")

# Loop through population pairs
for ((i=0; i<${#populations[@]}; i++)); do
    for ((j=i+1; j<${#populations[@]}; j++)); do
        pop1="${populations[$i]}"
        pop2="${populations[$j]}"
        echo "Calculating Fst for $pop1 vs $pop2"

        # Create output directory
        output_dir="Fst_results/${pop1}_vs_${pop2}"
        mkdir -p "$output_dir"

        # Loop through gene regions
        for gene in "${gene_regions[@]}"; do
            echo "Processing gene region: $gene"

            # Define VCF files
            vcf1="./$pop1/${pop1}_t2dvcfs/${pop1}_${gene}_t2dsnp.vcf.gz"
            vcf2="./$pop2/${pop2}_t2dvcfs/${pop2}_${gene}_t2dsnp.vcf.gz"
            merged_vcf="$output_dir/${pop1}_${pop2}_${gene}_t2dsnp.vcf.gz"

            # Ensure VCF files exist
            if [[ -f "$vcf1" && -f "$vcf2" ]]; then
                echo "Indexing VCF files..."
                bcftools index "$vcf1"
                bcftools index "$vcf2"

                echo "Merging VCFs for $pop1 vs $pop2 in $gene..."
                bcftools merge "$vcf1" "$vcf2" -Oz -o "$merged_vcf"

                # Check if merge was successful
                if [[ -f "$merged_vcf" ]]; then
                    echo "Running Fst for $pop1 vs $pop2 in $gene..."
                    vcftools --gzvcf "$merged_vcf" \
                             --weir-fst-pop "./population_samples/${pop1}_samples.txt" \
                             --weir-fst-pop "./population_samples/${pop2}_samples.txt" \
                             --out "$output_dir/${pop1}_${pop2}_${gene}_Fst"
                else
                    echo "Skipping Fst for $pop1 vs $pop2 in $gene (Merge failed)"
                fi

            else
                echo "Skipping $gene (VCF file missing for one or both populations)"
            fi
        done
    done
done

echo "Fst calculations completed!"

