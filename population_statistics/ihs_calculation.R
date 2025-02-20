# iHS Calculation Script in R

# Load necessary libraries
if (!requireNamespace("rehh", quietly = TRUE)) install.packages("rehh")
if (!requireNamespace("vcfR", quietly = TRUE)) install.packages("vcfR")
library(rehh)
library(vcfR)

# Base directory for VCF files
base_dir <- "~/Documents/BIO727P-GROUP_PROJECT/T2D_Web_development_project/population_statistics"

# List of populations and genes
populations <- c("BEB", "GIH", "ITU", "PJL", "STU")
genes <- c("ADCY5", "C3orf65", "CDC123", "CDKAL1", "CDKN1C", "CDKN2B", "COBLL1", "HHEX", "HNF4A", "IGF2BP2", "SLC30A8", "TCF7L2")

# Loop through each population and gene
for (pop in populations) {
  for (gene in genes) {
    vcf_file <- file.path(base_dir, pop, paste0(pop, "_", gene, ".vcf.gz"))
    output_dir <- file.path(base_dir, paste0(pop, "_", gene, "_ihs"))
    if (!dir.exists(output_dir)) dir.create(output_dir)
    
    cat("Processing", vcf_file, "\n")
    
    hap_data <- data2haplohh(hap_file = vcf_file, vcf_reader = "vcfR")
    
    if (is.null(hap_data)) {
      cat("Skipping", vcf_file, "due to errors.\n")
      next
    }
    
    scan_data <- scan_hh(hap_data)
    ihs_result <- ihh2ihs(scan_data)
    
    output_file <- file.path(output_dir, paste0(pop, "_", gene, "_ihs_results.txt"))
    write.table(ihs_result$ihs, file = output_file, sep = "\t", row.names = FALSE, quote = FALSE)
    
    cat("iHS results saved to", output_file, "\n")
  }
}

cat("iHS calculation completed for all populations and genes.\n")