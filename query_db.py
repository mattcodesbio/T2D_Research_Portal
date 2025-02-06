from app import create_app  # Import create_app from your app.py
from models import db, SNP, Gene, Population  # Import necessary models

app = create_app()  # Initialize the app using create_app

with app.app_context():  # Set up the app context to interact with the database
    # Query SNP data from the database
    snp = SNP.query.filter_by(rsID='rs12345').first()
    
    if snp:
        print(f"SNP: {snp.rsID}, Chromosome: {snp.chromosome}")
    else:
        print("SNP not found!")

