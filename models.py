from flask_sqlalchemy import SQLAlchemy

db = SQLAlchemy()

class SNP(db.Model):
    __tablename__ = 'snp'
    id = db.Column(db.Integer, primary_key=True)
    snp_id = db.Column(db.String(50), unique=True, nullable=False)
    chromosome = db.Column(db.String(10), nullable=False)
    position = db.Column(db.Integer, nullable=False)
    gene = db.Column(db.String(50), nullable=False)
    ontology = db.Column(db.String(50), nullable=False)
    population = db.Column(db.String(100))
    region = db.Column(db.String(100))

    def __repr__(self):
        return f"<SNP {self.snp_id}>"

class Gene(db.Model):
