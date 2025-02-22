from flask import Flask, render_template, request, jsonify
from flask_sqlalchemy import SQLAlchemy
from functions import load_snps_from_csv, load_tajima_d_results

# Initialize the Flask app
app = Flask(__name__)

# Configuration
app.config['SECRET_KEY'] = 'your_secret_key'
app.config['SQLALCHEMY_DATABASE_URI'] = 'sqlite:///database.db'
app.config['SQLALCHEMY_TRACK_MODIFICATIONS'] = False

# Initialize the database with Flask
from models import db, SNP, TajimaD, Fst  # Import db after initializing app

db.init_app(app)  # This ensures the SQLAlchemy instance is tied to the app

# Create database tables if they don’t exist
with app.app_context():
    db.create_all()

    # Load SNP and positive selection data into the database
    # load_snps_from_csv("updated_snp_data_with_mapped_genes.csv")
    # load_tajima_d_results("population_statistics/tajima_d_10kb_results")
    # load_fst_results("population_statistics/fst_results")

# Import routes AFTER initializing db to avoid circular imports
from routes import *

if __name__ == '__main__':
    app.run(debug=True)