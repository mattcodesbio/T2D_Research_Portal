import pandas as pd
from models import db, SNP
import re

# Function to retrieve SNP information from the database

from models import db, SNP

def get_snp_info(snp_id=None, chromosome=None, start=None, end=None, gene_name=None):
    query = SNP.query

    # Filtering conditions based on user input
    if snp_id:
        query = query.filter(SNP.snp_id.ilike(f"%{snp_id}%"))  # Case-insensitive search
    if chromosome:
        query = query.filter(SNP.chromosome.like(f"{chromosome}%"))
    if start:
        query = query.filter(SNP.start_position == start)
    if end:
        query = query.filter(SNP.end_position == end)
    if gene_name:
        query = query.filter(SNP.gene.ilike(f"%{gene_name}%"))  # Case-insensitive search

    snps = query.all()

    # Return a list of SNP dictionaries without odds_ratio and reference_id
    return [{
        'snp_id': snp.snp_id,
        'gene': snp.gene if snp.gene else 'N/A',
        'chromosome': snp.chromosome,
        'start_position': snp.start_position if snp.start_position else 'N/A',
        'end_position': snp.end_position if snp.end_position else 'N/A',
        'risk_allele': snp.risk_allele,
        'p_value': snp.p_value if snp.p_value else 'N/A',
        'population': snp.population
    } for snp in snps] if snps else None


# Function to load SNPs into the database from CSV

def clean_p_value(value):
    match = re.search(r"[\d\.]+", str(value))  # Extract numeric part
    return float(match.group()) if match else None  # Convert to float or None

from sqlalchemy.exc import IntegrityError

def load_snps_from_csv(csv_file):

    from main import app

    df = pd.read_csv(csv_file)

    df['snp_id'] = df['snp_id'].str.strip().str.lower()
    df['population'] = df['population'].str.strip().str.lower()

    # Remove duplicates based on snp_id and population before inserting
    df.drop_duplicates(subset=["snp_id", "population"], keep="first", inplace=True)

    # Handle missing risk_allele values
    df['risk_allele'].fillna('N/A', inplace=True)  # Replace NaN with 'N/A'
    df = df[df['risk_allele'] != '']  # Remove rows with empty risk_allele

    # Ensure chromosome column is filled
    df = df[df['chromosome'] != 'Unknown']  # Or you can skip rows with no chromosome info using df.dropna(subset=["chromosome"])
    df = df[df['chromosome'] != 'N/A'] 

    with app.app_context():
        for _, row in df.iterrows():
            try:
                # Check if SNP with same snp_id and population already exists
                existing_snp = SNP.query.filter_by(snp_id=row["snp_id"], population=row["population"]).first()
                if existing_snp:
                    print(f"SNP {row['snp_id']} for population {row['population']} already exists.")
                    continue  # Skip the existing SNP

                # Create a new SNP entry
                snp_entry = SNP(
                    snp_id=row["snp_id"],
                    gene=row["gene"],
                    chromosome=row["chromosome"],
                    start_position=row["start_position"],
                    end_position=row["end_position"],
                    risk_allele=row["risk_allele"],
                    p_value=clean_p_value(row["p_value"]),
                    population=row["population"],
                )
                db.session.add(snp_entry)

            except Exception as e:
                print(f"Error processing SNP {row['snp_id']}: {e}")
        
        db.session.commit()  # Commit changes after all rows are processed
        print("SNP data loaded successfully.")



# Run this function once to populate the database
if __name__ == "__main__":
    load_snps_from_csv("snps_with_positions.csv")



