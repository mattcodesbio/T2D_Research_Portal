from flask import Flask
from models import db, SNP
from functions import load_snps_from_csv
import os

app = Flask(__name__)
app.config['SQLALCHEMY_DATABASE_URI'] = 'sqlite:///database.db'
app.config['SQLALCHEMY_TRACK_MODIFICATIONS'] = False
app.config['SECRET_KEY'] = 'your_secret_key'

# Initialize database
db.init_app(app)

# Import routes AFTER initializing db to avoid circular imports
with app.app_context():
    from routes import *

# Create tables and load data
with app.app_context():
    db.create_all()
    if not SNP.query.first():
        csv_path = os.path.join("DATA", "SEC_DATA.csv") 
        print(f"Loading CSV from: {csv_path}") # debugging purposes 
        load_snps_from_csv(csv_path)

if __name__ == '__main__':
    app.run(debug=True)