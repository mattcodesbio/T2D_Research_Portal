from flask import Flask, render_template
from flask_sqlalchemy import SQLAlchemy

# Initialize the app
app = Flask(__name__)

# Configuration
app.config['SECRET_KEY'] = 'your_secret_key'
app.config['SQLALCHEMY_DATABASE_URI'] = 'sqlite:///database.db'

# Initialize the database
db = SQLAlchemy(app)

# Import routes
from routes import *

if __name__ == '__main__':
    app.run(debug=True) 