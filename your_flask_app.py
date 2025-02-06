from flask import Flask
from flask_sqlalchemy import SQLAlchemy

app = Flask(__name__)

# Database configuration
app.config['SQLALCHEMY_DATABASE_URI'] = 'sqlite:///t2d_project.db'  # Path to  database
app.config['SQLALCHEMY_TRACK_MODIFICATIONS'] = False  # Disable modification tracking

# Initialize the database object
db = SQLAlchemy(app)

# Define models here e.g.:
class SNP(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    rsID = db.Column(db.String(20), unique=True, nullable=False)
    chromosome = db.Column(db.String(5))
    position_start = db.Column(db.Integer)
    position_end = db.Column(db.Integer)
    p_value = db.Column(db.Float)

# Create tables in the database (only run once or after making changes to models)
with app.app_context():
    db.create_all()

