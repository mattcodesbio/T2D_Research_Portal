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
