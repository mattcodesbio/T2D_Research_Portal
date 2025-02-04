import pandas as pd


# Cleaning data csv:

 # Load your data
ps_data = pd.read_csv('snps_with_positions.csv')

# Check for missing or invalid values  # Drop rows with missing values
ps_data['case_control'] = ps_data['case_control'].apply(lambda x: x.split(' ; '))  # Split case/control counts


# calculate Minor Allele Frequency:

# if rare alleles are under selection a lower MAF would indicate a positive selection signal 

