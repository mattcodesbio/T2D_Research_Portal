from flask_sqlalchemy import SQLAlchemy

db = SQLAlchemy()

class SNP(db.Model):
    __tablename__ = 'snps'
    
    snp_id = db.Column(db.String(20), primary_key=True, nullable=False)  # SNP ID (e.g., rs123456)
    chromosome = db.Column(db.String(5), nullable=False)  # Chromosome
    grch38_start = db.Column(db.Integer, nullable=False)  # Genomic position (GRCh38)
    gene_name = db.Column(db.String(100), nullable=True)  # Gene name(s)
    p_value = db.Column(db.Float, nullable=True)  # P-value from GWAS
    reference_allele = db.Column(db.String(10), nullable=True)  # Reference Allele
    alternative_allele = db.Column(db.String(10), nullable=True)  # Alternative Allele
    consequence = db.Column(db.String(50), nullable=True)  # Variant consequence

    def __repr__(self):
        return f"<SNP {self.snp_id}, Chromosome {self.chromosome}, Position {self.grch38_start}>"

# Tajima's D Table
class TajimaD(db.Model):
    __tablename__ = 'tajima_d_results'
    
    id = db.Column(db.Integer, primary_key=True)
    population = db.Column(db.String(50), nullable=False)  # Population name
    chromosome = db.Column(db.String(5), nullable=False)  # Chromosome
    bin_start = db.Column(db.Integer, nullable=False)  # Start position of the bin
    bin_end = db.Column(db.Integer, nullable=False)  # End position of the bin
    n_snps = db.Column(db.Integer, nullable=False)  # Number of SNPs in the window
    tajima_d = db.Column(db.Float, nullable=False)  # Tajima's D value

    def __repr__(self):
        return f"<TajimaD {self.population} Chromosome {self.chromosome} Bin {self.bin_start}-{self.bin_end}>"

