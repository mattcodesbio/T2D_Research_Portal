from flask import Flask, jsonify, request
from flask_sqlalchemy import SQLAlchemy
import pandas as pd
import os

# Load Excel file
excel_file = 'FILE.xlsx'  # Change this to your file path
df = pd.read_excel(excel_file)

# Convert the DataFrame to CSV
df.to_csv('file.csv', index=False)  # This will save the file as CSV without the index


app = Flask(__name__)

# Setup the database URI
app.config['SQLALCHEMY_DATABASE_URI'] = 'sqlite:///database.db'  # SQLite file-based DB
app.config['SQLALCHEMY_TRACK_MODIFICATIONS'] = False  # Disable track modifications

# Initialize the database
db = SQLAlchemy(app)

# Define the model based on the structure of your CSV
class DataEntry(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    snp_id = db.Column(db.String(20), primary_key=True)
    chromosome = db.Column(db.Integer, nullable=False)
    population = db.Column(db.String(60), nullable=False)
    position = db.Column(db.Integer, nullable=False)
    gene_id = db.Column(db.String(20), db.ForeignKey('gene.gene_id'))
    p_value = db.Column(db.Float)
    positive_selection_stat = db.Column(db.Float, nullable=True)
    external_link = db.Column(db.String(200), nullable=True)

    def __init__(self, column1, column2):
        self.snp_id = snp_id
        self.chromosome = chromosome
        self.population = population
        self.position = position
        self.gene_id = gene_id
        self.p_value = p_value
        self.positive_selection_stat = positive_selection_stat
        self.external_link = external_link 

# Route to fetch all entries
@app.route('/data', methods=['GET'])
def get_data():
    entries = DataEntry.query.all()
    return jsonify([entry.to_dict() for entry in entries])

# Convert DataFrame to database entries
def load_csv_to_db(csv_file):
    df = pd.read_csv(csv_file)
    for _, row in df.iterrows():
        # Assuming the DataFrame columns match the model's columns
        entry = DataEntry(column1=row['column1'], column2=row['column2'])
        db.session.add(entry)
    db.session.commit()

@app.before_first_request
def create_tables():
    db.create_all()

# Add a simple route to trigger CSV upload
@app.route('/upload_csv', methods=['POST'])
def upload_csv():
    csv_file = request.files['file']
    csv_file.save('uploaded_file.csv')  # Save file locally
    load_csv_to_db('uploaded_file.csv')  # Load CSV data into DB
    return 'CSV data loaded into database!', 200

if __name__ == '__main__':
    app.run(debug=True)

