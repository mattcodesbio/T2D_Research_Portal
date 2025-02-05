from docx import Document
import pandas as pd
import requests
import time
import re

# Load the .docx document
doc_path = "/Users/Coding/Downloads/mmc1.docx"  # Update with your file path
doc = Document(doc_path)

# Extract tables from the document
tables = doc.tables

# Function to extract table data
def extract_table_data(table):
    data = []
    for row in table.rows:
        data.append([cell.text.strip() for cell in row.cells])
    return data

# Iterate through all tables and extract data
table_data = [extract_table_data(table) for table in tables]

# Identify the correct table (manually check)
for idx, table in enumerate(table_data):
    print(f"Table {idx + 1}: First row -> {table[0] if table else 'Empty'}")

# Select the correct table index (Update this based on manual inspection)
snp_table_index = 3  # Change this to the actual index found
snp_data = table_data[snp_table_index]

# Define column headers based on the expected format
column_names = ["gene", "chromosome", "snp_id", "effect", "risk_allele", "detection_method",
                "population", "case_control", "control_hwe_p_value", "odds_ratio",
                "p_value", "reference_id"]

# Ensure we have enough columns before applying headers
if len(snp_data) > 1:  # Only apply headers if there's actual data
    df = pd.DataFrame(snp_data[1:], columns=column_names)  # Skip first row (header)
else:
    df = pd.DataFrame(snp_data, columns=column_names)  # Use full data if only one row

# Function to extract the odds ratio before the brackets
def extract_odds_ratio(odds_ratio_str):
    match = re.match(r'([0-9.]+)\(', odds_ratio_str)  # Match the number before the opening parenthesis
    if match:
        return float(match.group(1))  # Convert it to float
    return None  # Return None if the format doesn't match

# Apply the function to the 'odds_ratio' column
df['odds_ratio'] = df['odds_ratio'].apply(extract_odds_ratio)

# Convert numeric columns safely
df["odds_ratio"] = pd.to_numeric(df["odds_ratio"], errors="coerce")
# Handle '<' values in the 'p_value' column (e.g., <0.0001) and convert others to numeric
df['p_value'] = df['p_value'].apply(lambda x: x if isinstance(x, str) and '<' in x else pd.to_numeric(x, errors='coerce'))

# Replace NaN values with "N/A"
df['p_value'].fillna("N/A", inplace=True)
df["reference_id"] = pd.to_numeric(df["reference_id"], errors="coerce")

# **Replace missing values with "N/A"**
df.fillna("N/A", inplace=True)
# **Drop rows with missing SNP IDs**
df.dropna(subset=['snp_id'], inplace=True)

# **Convert empty SNP IDs to "N/A" for debugging**
df['snp_id'].replace('', 'N/A', inplace=True)

# Display first few rows
print("\nExtracted SNPs Preview:")
print(df.head())

# Handle missing values
print("\nMissing values per column AFTER replacement:")
print(df.isnull().sum())  # Should now be 0 for all fields

# Drop duplicates, keeping only the first occurrence
df.drop_duplicates(subset=["snp_id", "population"], keep="first", inplace=True)

# Check for exact duplicates across all columns
exact_duplicates = df[df.duplicated()]
print("\n🔍 **Exact Duplicates Across All Columns:**")
print(exact_duplicates)

# Check for duplicates based on snp_id and population (which causes SQL error)
dup_snp_population = df[df.duplicated(subset=["snp_id", "population"], keep=False)]
print("\n⚠️ **Duplicate SNPs for the Same Population:**")
print(dup_snp_population)

# Function to fetch SNP position (start and end) from Ensembl using SNP ID
def fetch_snp_position(snp_id):
    url = f"https://rest.ensembl.org/variation/human/{snp_id}?content-type=application/json"
    response = requests.get(url)

    if response.status_code == 200:
        snp_data = response.json()
        if 'mappings' in snp_data:
            for mapping in snp_data['mappings']:
                chromosome = mapping['seq_region_name']
                start_position = mapping['start']
                end_position = mapping['end']
                return chromosome, start_position, end_position
        else:
            print(f"No mapping data found for SNP ID {snp_id}")
            return None
    else:
        print(f"Error fetching SNP data for {snp_id}")
        return None

# Adding start and end position to the DataFrame
df['chromosome'] = df['chromosome'].fillna('')  # Handle any empty chromosomes (just in case)
df['start_position'] = None
df['end_position'] = None

# Loop through each SNP and fetch its position data
for index, row in df.iterrows():
    snp_id = row['snp_id']
    print(f"Fetching position for SNP ID {snp_id}")
    position_data = fetch_snp_position(snp_id)
    
    if position_data:
        chromosome, start, end = position_data
        df.at[index, 'chromosome'] = chromosome
        df.at[index, 'start_position'] = start
        df.at[index, 'end_position'] = end
    else: 
        # Handle cases where no data is returned
        df.at[index, 'chromosome'] = 'Unknown'
        df.at[index, 'start_position'] = 'Unknown'
        df.at[index, 'end_position'] = 'Unknown'
    
    # Sleep to avoid hitting API rate limit
    time.sleep(0.5)

unknown_chromosomes_count = df[df['chromosome'] == 'Unknown'].shape[0]
print(f"Number of SNPs with 'Unknown' chromosome: {unknown_chromosomes_count}")

# Check the updated DataFrame
print("\nUpdated SNPs Preview with Positions:")
print(df.head())

# Save the updated dataframe as CSV for database import
df.to_csv("snps_with_positions.csv", index=False)

# Apply the cleaning to the 'population' column to keep only the second-to-last and last parts
df['population'] = df['population'].apply(lambda x: '/'.join(x.split('/')[-2:]) if len(x.split('/')) > 1 else x)

# Preview the updated DataFrame
print(df['population'].head())

# Save the updated dataframe as CSV for database import
df.to_csv("snps_with_positions.csv", index=False)
