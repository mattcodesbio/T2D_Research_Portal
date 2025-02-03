from main import app
from flask import render_template, request
from functions import get_snp_info

@app.route('/', methods=['GET', 'POST'])
def home():
    if request.method == 'POST':
        # Retrieving form data
        snp_name = request.form.get('snp_name')
        chromosome = request.form.get('chromosome')
        start = request.form.get('start')
        end = request.form.get('end')
        gene_name = request.form.get('gene_name')

        print(f"Form Data: SNP={snp_name}, Chromosome={chromosome}, Start={start}, End={end}, Gene={gene_name}")
        
        # Process the form data by calling `get_snp_info`
        snp_info = get_snp_info(snp_name, chromosome, start, end, gene_name)
        print(snp_info) 

        # Render SNP results if information is found
        if snp_info:
            return render_template('snp_results.html', snp_info=snp_info)
        else:
            return render_template('snp_results.html', message="No SNP information found for your query.")

    # Render the home page for GET requests
    return render_template('index.html')

@app.route('/about')
def about():
    return render_template('about.html')


