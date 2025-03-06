rm(list = ls())

# Load libraries
library(tidyverse)
library(data.table)

# Define parameters
base_dir <- "/path/to/SweeD_output_file/"
block_size <- 50  # Num. of positions per block for spatial adjustment
conf_level <- 1.96  # z-score for 95 CI

# Initialise vectors to store CLR values 
clr_values_all <- c()
block_means <- c()

# Process SweeD report files 
file_paths <- list.files(base_dir, 
                         pattern = "SweeD_Report.*$", 
                         full.names = TRUE, 
                         recursive = TRUE)

for (file_path in file_paths) {
  # Read and sort genomic data
  df <- fread(file_path) %>%
    arrange(Position) 
  
  clr_scores <- df$Likelihood
  clr_values_all <- c(clr_values_all, clr_scores)
  
  # Calculate block means per chromosome
  n_blocks <- length(clr_scores) %/% block_size
  if(n_blocks > 0) {
    blocks <- split(clr_scores, 
                    cut(seq_along(clr_scores), 
                        breaks = n_blocks, 
                        labels = FALSE))
    chr_block_means <- sapply(blocks, mean, na.rm = TRUE)
    block_means <- c(block_means, chr_block_means)
  }
}

# Calculate significance thresholds
threshold <- quantile(clr_values_all, 0.95, na.rm = TRUE)
se <- sd(block_means)/sqrt(length(block_means))
adjusted_threshold <- threshold + conf_level * se

# Output results
cat(str_glue("
  Unadjusted 95th percentile: {round(threshold, 2)}
  Block SE: {round(se, 2)}
  Spatially adjusted threshold: {round(adjusted_threshold, 2)}
  Effective blocks: {length(block_means)}
"))

# Visualise distribution of CLR values 
ggplot(data.frame(CLR = clr_values_all), aes(x = CLR)) +
  geom_histogram(bins = 30, fill = "steelblue", alpha = 0.8) +
  geom_vline(xintercept = adjusted_threshold, 
             color = "red", 
             linetype = "dashed") +
  labs(title = "CLR Distribution with Spatially Adjusted Threshold",
       x = "Composite Likelihood Ratio",
       y = "Frequency") +
  theme_minimal()
