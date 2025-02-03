from main import app
from flask import render_template, request
import requests, sys
from functions import get_snp_info
import re
from urllib.parse import quote_plus  # For URL encoding

@app.route('/', methods=['GET', 'POST'])
def home():
    if request.method == 'POST':
        snp_id = request.form.get('snp_name')
        chromosome = request.form.get('chromosome')
        start = request.form.get('start')
        end = request.form.get('end')   
        gene_name = request.form.get('gene_name')

        # Extract only the numeric part of the chromosome for searching
        if chromosome:
            chromosome = ''.join(filter(str.isdigit, chromosome))  # Extract number before 'q'

        print(f"Search Query: SNP={snp_id}, Chromosome={chromosome}, Start={start}, End={end}, Gene={gene_name}")

        # Get results
        snp_info = get_snp_info(snp_id, chromosome, start, end, gene_name)
        return render_template('snp_results.html', snp_info=snp_info)

    return render_template('index.html')


@app.route('/about')
def about():
    return render_template('about.html')


# Define your route for displaying gene terms
@app.route('/gene_terms/<gene_name>')
def gene_terms(gene_name):
    # Extract gene name inside parentheses if any
    gene_name = extract_gene_name(gene_name)

    # Get GO terms for the gene
    go_terms = get_gene_ontology_terms(gene_name)

    if go_terms:
        return render_template('ontology.html', gene_name=gene_name, terms=go_terms)
    else:
        return render_template('ontology.html', gene_name=gene_name, terms=None)
    

def extract_gene_name(gene_name):
    # Regular expression to capture the portion inside parentheses
    match = re.search(r'\((.*?)\)', gene_name)
    if match:
        return match.group(1)  # Return the gene name inside the parentheses
    return gene_name  # Return the original gene name if no parentheses are found


# Function to fetch gene ontology terms from Ensembl API
def get_gene_ontology_terms(gene_name):
    server = "https://rest.ensembl.org"
    ext = f"/lookup/symbol/human/{gene_name}?content-type=application/json"  # Modify to include the gene name

    # Sending GET request to Ensembl API
    r = requests.get(server + ext, headers={"Content-Type": "application/json"})

    # Handling error if the request fails
    if not r.ok:
        r.raise_for_status()

    # Parsing the response from the API
    decoded = r.json()

    # Print out the entire decoded response for inspection
    print(f"API Response for {gene_name}: {decoded}")

    # Check if 'go_terms' is part of the response and return it
    if 'go_terms' in decoded and decoded['go_terms']:
        return decoded['go_terms']
    else:
        return f"No Gene Ontology terms found for gene: {gene_name}"
