from models import db, SNP
import pandas as pd
import re

# Function to retrieve SNP information from the database
def get_snp_info(dbSNP=None, chromosome=None, position=None, gene=None, position_range=None):
    query = SNP.query

    # Filtering conditions based on user input
    if dbSNP:
        query = query.filter(SNP.dbSNP.ilike(f"%{dbSNP}%"))  # Case-insensitive search
    if chromosome:
        query = query.filter(SNP.Chromosome == str(chromosome))
    if gene:
        query = query.filter(SNP.Gene.ilike(f"%{gene}%"))  # Case-insensitive gene search
    if position:
        query = query.filter(SNP.Position == position)
    if position_range:  # Tuple (start, end)
        query = query.filter(SNP.Position.between(*position_range))

    snps = query.all()

    # Return list of SNP dictionaries with relevant fields
    return [{
        'dbSNP': snp.dbSNP,
        'Gene': snp.Gene,
        'Chromosome': snp.Chromosome,
        'Position': snp.Position,
        'P_Value': f"{snp.P_Value:.2e}",  # Format in scientific notation
        'Minor_Allele': snp.Minor_Allele,
        'zScore': round(snp.zScore, 2) if snp.zScore else None,
        'varId': snp.varId
    } for snp in snps] if snps else []

# load SNPs from CSV
def clean_p_value(value):
    """Extract numeric value from potential string formatting"""
    try:
        return float(re.search(r"[\d\.eE+-]+", str(value)).group())
    except (AttributeError, ValueError):
        return None


from models import db, SNP
import pandas as pd
import re
from flask import current_app  # Use current_app instead of importing app

def load_snps_from_csv(csv_file):
    df = pd.read_csv(csv_file)
    
    # Clean data
    df = df.drop_duplicates(subset=["dbSNP"])
    df['P_Value'] = df['P_Value'].apply(clean_p_value)
    df['zScore'] = pd.to_numeric(df['zScore'], errors='coerce')
    df['Position'] = pd.to_numeric(df['Position'], errors='coerce').astype('Int64')
    df = df.dropna(subset=["dbSNP", "Chromosome", "Position", "P_Value"])

    # Use Flask's application context
    with current_app.app_context():
        try:
            db.session.bulk_insert_mappings(SNP, df.to_dict(orient='records'))
            db.session.commit()
            print(f"Loaded {len(df)} SNPs from {csv_file}")
        except Exception as e:
            db.session.rollback()
            print(f"Error: {str(e)}")


# populate the database
if __name__ == "__main__":
    load_snps_from_csv("DATA/SEC_DATA.csv")
