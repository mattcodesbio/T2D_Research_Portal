# All data used here was extracted from: https://zenodo.org/records/3234689
# Filter the data for the SAS subsets
BEB_data <- subset(pop_size, population == "BEB")
STU_data <- subset(pop_size, population == "STU")
PJL_data <- subset(pop_size, population == "PJL")
ITU_data <- subset(pop_size, population == "ITU")
GIH_data <- subset(pop_size, population == "GIH")

# Population data list
pop_subsets_list <- list(BEB_data, STU_data, PJL_data, ITU_data, GIH_data)

# Convert gens_ago and population_size to log base 10 for each population
pop_subsets_list <- lapply(pop_subsets_list, function(data) 
  {data$population_size <- log10(data$population_size)
  data$gens_ago <- log10(data$gens_ago)
  return(data) 
  })

# Assign the converted data to the respective pop 
BEB_data <- pop_subsets_list[[1]]
STU_data <- pop_subsets_list[[2]]
PJL_data <- pop_subsets_list[[3]]
ITU_data <- pop_subsets_list[[4]]
GIH_data <- pop_subsets_list[[5]]

# Plot the staircase graph for the first population 
plot(BEB_data$gens_ago, BEB_data$population_size, type = "s", 
     main = "SAS Population Size Over Generations",
     xlab = "Timeline (Log10)", 
     ylab = "Population Size (Log10)", 
     col = "orange", 
     lwd = 2,
     xlim = c(2, max(c(BEB_data$gens_ago, STU_data$gens_ago, 
                           PJL_data$gens_ago, ITU_data$gens_ago, 
                           GIH_data$gens_ago))),  # Set x-axis starting from log10(100) 
     ylim = c(2, max(c(BEB_data$population_size, STU_data$population_size, 
                       PJL_data$population_size, ITU_data$population_size,
                       GIH_data$population_size)))) # Set y-axis limit starting at log10(100)

# Add grid
grid()

# Add the other populations 
# STU
lines(STU_data$gens_ago, STU_data$population_size, type = "s", 
      col = "coral1", lwd = 2)

# PJL
lines(PJL_data$gens_ago, PJL_data$population_size, type = "s", 
      col = "violetred", lwd = 2)

# ITU
lines(ITU_data$gens_ago, ITU_data$population_size, type = "s", 
      col = "turquoise2", lwd = 2)

# ('GIH') 
lines(GIH_data$gens_ago, GIH_data$population_size, type = "s", 
      col = "yellowgreen", lwd = 2)

# Legend 
legend("topright", legend = c("BEB", "STU", "PJL", "ITU", "GIH"), 
       col = c("orange", "coral1", "violetred", "turquoise2", "yellowgreen"), 
       lty = 1, 
       lwd = 2,
       cex = 0.7,  # Smaller text
       x.intersp = 0.3,  # Reduce horizontal spacing
       y.intersp = 0.4,  # Reduce vertical spacing
       inset = 0.02,  # Move inside the plot
       box.lwd = 0.3,  # Box outline thinner
       text.width = 0.1)  # Reduce extra space to the right
