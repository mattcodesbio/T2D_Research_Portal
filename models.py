from flask_sqlalchemy import SQLAlchemy

# Initialize SQLAlchemy instance
# This is shared across different files in the Flask app
# The instance is linked to the Flask app in 'main.py'
db = SQLAlchemy()




# Define SNP Database table schema
class SNP(db.Model):    
    """
    This table stores information about Type 2 Diabetes (T2D)-associated SNPs (Single Nucleotide Polymorphisms)
    obtained from a genome-wide association study (GWAS).
    """
    __tablename__ = 'snps'
    
    snp_id = db.Column(db.String(20), primary_key=True, nullable=False)  # SNP ID (e.g., rs123456)
    chromosome = db.Column(db.String(5), nullable=False, index = True)  # Chromosome
    grch38_start = db.Column(db.Integer, nullable=False, index = True)  # Genomic position based on GRCh38 reference genome
    gene_name = db.Column(db.String(100), nullable=True, index = True)  # Mapped gene(s) for the SNP (if available)
    p_value = db.Column(db.Float, nullable=True)  # P-value from GWAS study
    reference_allele = db.Column(db.String(10), nullable=True)  # Reference Allele
    alternative_allele = db.Column(db.String(10), nullable=True)  # Alternative Allele
    consequence = db.Column(db.String(50), nullable=True)  # Variant consequence

    fst_results = db.relationship('FstSNP', backref='snp', lazy=True)

    def __repr__(self):
        return f"<SNP {self.snp_id}, Chromosome {self.chromosome}, Position {self.grch38_start}>"


# Define Tajima's D database Table schema
class TajimaD(db.Model):
    """
    Stores Tajima's D statistics for different populations and genomic regions.
    """
    __tablename__ = 'tajima_d_results'
    
    id = db.Column(db.Integer, primary_key=True) # Unique ID for each record
    population = db.Column(db.String(50), nullable=False, index = True)  # Name of the population studied
    chromosome = db.Column(db.String(5), nullable=False, index = True)  # Chromosome
    bin_start = db.Column(db.Integer, nullable=False, index = True)  # Start position of the genomic region/bin
    bin_end = db.Column(db.Integer, nullable=False, index = True)  # End position of the bin (10kb window size)
    n_snps = db.Column(db.Integer, nullable=False) # Number of SNPs found within this genomic bin
    tajima_d = db.Column(db.Float, nullable=False)  # Tajima's D statistic value

    def __repr__(self):
        return f"<TajimaD {self.population} Chromosome {self.chromosome} Bin {self.bin_start}-{self.bin_end}>"


class FstSNP(db.Model):
    __tablename__ = 'fst_snp'
    id = db.Column(db.Integer, primary_key=True)
    snp_id = db.Column(db.String(20),db.ForeignKey('snps.snp_id'), nullable=False, index = True)
    gene = db.Column(db.String(50), nullable=False, index = True)
    chromosome = db.Column(db.String(5), nullable=True, index = True)
    position = db.Column(db.Integer, nullable=True, index = True)  # Fixed line
    fst_beb = db.Column(db.Float, nullable=False)
    fst_gih = db.Column(db.Float, nullable=False)
    fst_itu = db.Column(db.Float, nullable=False)
    fst_pjl = db.Column(db.Float, nullable=False)
    fst_stu = db.Column(db.Float, nullable=False)

    def __repr__(self):
        return f"<FstSNP {self.snp_id} Chr{self.chromosome}:{self.position}>"
    

# CLR Table 
class CLRTest(db.Model):
    """
    Stores Composite Likelihood Ratio (CLR) test results for selective sweeps.
    CLR is used to detect recent positive selection.
    """
    __tablename__ = 'clr_results'
    
    id = db.Column(db.Integer, primary_key=True)
    population = db.Column(db.String(50), nullable=False, index = True)
    chromosome = db.Column(db.String(5), nullable=False, index = True)
    position = db.Column(db.Integer, nullable=False, index = True)
    clr = db.Column(db.Float, nullable=False) # CLR value
    alpha = db.Column(db.Float, nullable=False) # Alpha statistic (Measures the strength of the likelihood/selection signal)

    def __repr__(self):
        return f"<CLRTest {self.population} Chr{self.chromosome}:{self.position} CLR={self.clr}>"

#### Important SQLAlchemy concepts ####

# 'db.Model' : This is the base class for definning a database table
# 'db.column' : Defines the a field in the databse table
# Primary keys ('primary_key=True') : Uniquely identifies each record in the table 
# Data types:
    # 'db.string(20)' : This field takes string datatypes with a max lengh of 20
    # 'db.Integer': This field takes integer values
    # 'db.Float': This field takes floating point values



