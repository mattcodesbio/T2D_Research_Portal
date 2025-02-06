from app import create_app, db
from models import SNP, Gene, Population

app = create_app()

with app.app_context():
    db.create_all()  # Create tables
    # Add sample data to tables, if needed
    snp = SNP(rsID='rs12345', chromosome='11')
    db.session.add(snp)

    gene = Gene(gene_name='INS')
    db.session.add(gene)

    population = Population(population_name='British Bangladeshi')
    db.session.add(population)

    db.session.commit()  # Commit to save changes

