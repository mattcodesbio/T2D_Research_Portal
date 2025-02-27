# main.py - The entry point for the Flask web application
# It initializes the Flask app, configures the database, and imports routes

from flask import Flask, render_template, request, jsonify
from flask_sqlalchemy import SQLAlchemy
from functions import load_snps_from_csv, load_tajima_d_results

# Initialize the Flask app
app = Flask(__name__)

# Configuration settings for Flask and the database
app.config['SECRET_KEY'] = 'your_secret_key' #Secret key used for security of the session
app.config['SQLALCHEMY_DATABASE_URI'] = 'sqlite:///database.db' #Database connection using SQLite
app.config['SQLALCHEMY_TRACK_MODIFICATIONS'] = False # Disables modification tracking to improve performance

# Initialize the database with Flask
from models import db, SNP, TajimaD  # Import database models after initializing Flask app

db.init_app(app)  # This ensures the SQLAlchemy database is linked to the Flask app

# Create database tables if they don’t exist
with app.app_context():
    db.create_all()
            
    # Uncomment the following lines if you want to preload SNP and Tajima's D data into the database
    # load_snps_from_csv("/path/to/file/e.g/csv")
    load_tajima_d_results("/mnt/c/Users/hleeh/Documents/BIO727P/tajima_d_10kb_chr2/tajima_d_10kb_chr2")
    

# Import routes AFTER initializing db to avoid circular imports
from routes import *

# Start the Flask app 
if __name__ == '__main__':
    app.run(debug=True) # Runs the Flask application in debug mode

