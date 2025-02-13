# from flask_sqlalchemy import SQLAlchemy

# db = SQLAlchemy()  # Initialize the database

# class SNP(db.Model):
#     __tablename__ = 'snps'
    
#     snp_id = db.Column(db.String(20), nullable=False)  # SNP ID (e.g., rs123456)
#     gene = db.Column(db.String(50), nullable=True)       # Associated gene
#     chromosome = db.Column(db.String(20), nullable=False)  # Chromosome e.g., '9q31.1'
#     start_position = db.Column(db.Integer, nullable=True)  # Genomic start position
#     end_position = db.Column(db.Integer, nullable=True)    # Genomic end position
#     risk_allele = db.Column(db.String(5), nullable=False)  # Risk allele
#     p_value = db.Column(db.Float, nullable=True)  # P-value
#     population = db.Column(db.String(100), nullable=False)  # Study population

#     # Use a composite primary key (snp_id + population)
#     __table_args__ = (db.PrimaryKeyConstraint('snp_id', 'population'),)

#     def __repr__(self):
#         return f"<SNP {self.snp_id}, Gene: {self.gene}, Chromosome: {self.chromosome}>"
from flask_sqlalchemy import SQLAlchemy

db = SQLAlchemy()

class SNP(db.Model):
    __tablename__ = 'snps'

    
    Gene = db.Column(db.String(50), nullable=True)           # Gene
    Chromosome = db.Column(db.String(20), nullable=False)    # Chromosome
    dbSNP = db.Column(db.String(20), primary_key=True)       # dbSNP (primary key)
    Reference = db.Column(db.String(1), nullable=False)      # Reference
    Alternate = db.Column(db.String(1), nullable=False)      # Alternate
    Position = db.Column(db.Integer, nullable=False)         # Position
    P_Value = db.Column(db.Float, nullable=False)            # P_Value
    Minor_Allele = db.Column(db.String(1), nullable=False)   # Minor_Allele
    varId = db.Column(db.String(50), unique=True)            # varId
    zScore = db.Column(db.Float)                             # zScore

    def __repr__(self):
        return f"<SNP {self.dbSNP} ({self.Gene})>"
