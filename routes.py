
# from flask import render_template, request
# from main import app  
# from functions import get_snp_info
# import requests, re
# from urllib.parse import quote_plus

# @app.route('/', methods=['GET', 'POST'])
# def home():
#     if request.method == 'POST':
#         # Get form inputs
#         snp_id = request.form.get('snp_name')
#         chromosome = request.form.get('chromosome')
#         start = request.form.get('start')
#         end = request.form.get('end')
#         gene_name = request.form.get('gene_name')

#         # Clean chromosome input (extract numeric part)
#         if chromosome:
#             chromosome = ''.join(filter(str.isdigit, chromosome))

#         # Process position range
#         position_range = None
#         position = None
#         try:
#             if start and end:
#                 position_range = (int(start), int(end))
#             elif start:
#                 position = int(start)
#         except ValueError:
#             pass  # Handle invalid numeric input

#         # Query SNPs
#         snp_info = get_snp_info(
#             dbSNP=snp_id,
#             chromosome=chromosome,
#             position=position,
#             gene=gene_name,
#             position_range=position_range
#         )

#         return render_template('snp_results.html', snp_info=snp_info)
    
#     return render_template('index.html')

# @app.route('/about')
# def about():
#     return render_template('about.html')


# @app.route('/gene_terms/<gene_name>')
# def gene_terms(gene_name):
#     # Extract gene name inside parentheses if any
#     gene_name = extract_gene_name(gene_name)

#     # Get GO terms for the gene using Ensembl API
#     go_terms = get_gene_ontology_terms(gene_name)

#     if go_terms:
#         terms = []
        
#         # Loop through all terms and extract relevant information
#         for term in go_terms:
#             term_data = {
#                 'name': term.get('name', 'N/A'),
#                 'accession': term.get('accession', 'N/A'),
#                 'namespace': term.get('namespace', 'N/A'),
#                 'definition': term.get('definition', 'No definition available'),
#                 'synonyms': term.get('synonyms', ['No synonyms available']),
#                 'children': term.get('children', []),
#                 'parents': term.get('parents', [])  # Collect parent terms for each entry
#             }
#             terms.append(term_data)

#         # Collect parent terms from all the terms
#         all_parents = []
#         for term in terms:
#             all_parents.extend(term['parents'])  # Flatten the list of parent terms

#         return render_template('ontology.html', gene_name=gene_name, terms=terms, all_parents=all_parents)
#     else:
#         return render_template('ontology.html', gene_name=gene_name, terms=None, all_parents=None)


# def extract_gene_name(gene_name):
#     # Extract the gene symbol from parentheses, e.g., "Adiponectin, C1Q and Collagen Domain Containing (ADIPOQ)" -> "ADIPOQ"
#     match = re.search(r'\((.*?)\)', gene_name)
#     if match:
#         return match.group(1).strip()  # Return the gene symbol inside the parentheses
#     return gene_name  # Return the original gene name if no parentheses are found

# def get_gene_ontology_terms(gene_name):
#     url = f"https://rest.ensembl.org/ontology/name/{quote_plus(gene_name)}?content-type=application/json"
#     response = requests.get(url)
    
#     if response.status_code == 200:
#         go_terms = response.json()
#         return go_terms  # Return the entire JSON data of GO terms
#     else:
#         return None  # Handle failed requests

# @app.route('/select_population', methods=['POST'])
# def select_population():
#     selected_populations = request.form.getlist('selected_populations')
    
#     # Assuming we have a dictionary or function that maps population names to geographical information
#     population_info = get_population_info(selected_populations)
    
#     return render_template('population_info.html', population_info=population_info)

# def get_population_info(populations):
#     # Example dictionary mapping populations to their geographical descriptions
#     population_map = {
#         'haryana/northern india': "Haryana is a northern state of India, known for its agricultural economy and historical sites.",
#         'tamilnadu/southern india': "Tamil Nadu, India, known for its industrial growth and proximity to the Western Ghats."
#         # Add more population descriptions as needed
#     }

#     # Collect the information for selected populations
#     info = {}
#     for population in populations:
#         if population in population_map:
#             info[population] = population_map[population]
    
#     return info
from flask import render_template, request, abort
from main import app  
from functions import get_snp_info
import requests, re
from urllib.parse import quote_plus

@app.route('/', methods=['GET', 'POST'])
def home():
    if request.method == 'POST':
        # Get form inputs
        snp_id = request.form.get('snp_name')
        chromosome = request.form.get('chromosome')
        start = request.form.get('start')
        end = request.form.get('end')
        gene_name = request.form.get('gene_name')

        # Clean chromosome input
        if chromosome:
            chromosome = ''.join(filter(str.isdigit, chromosome))

        # Validate and process position inputs
        position_range = None
        position = None
        try:
            if start and end:
                # Convert to integers and validate range
                start_int = int(start)
                end_int = int(end)
                if end_int < start_int:
                    abort(400, "End position cannot be smaller than start position")
                position_range = (start_int, end_int)
            elif start:
                position = int(start)
        except ValueError:
            abort(400, "Invalid position values - must be integers")

        # Query SNPs
        snp_info = get_snp_info(
            dbSNP=snp_id,
            chromosome=chromosome,
            position=position,
            gene=gene_name,
            position_range=position_range
        )

        return render_template('snp_results.html', snp_info=snp_info)
    
    # GET request - show empty form
    return render_template('index.html')

@app.route('/about')
def about():
    return render_template('about.html')

@app.route('/gene_terms/<gene_name>')
def gene_terms(gene_name):
    gene_name = extract_gene_name(gene_name)
    go_terms = get_gene_ontology_terms(gene_name)

    if go_terms:
        terms = []
        all_parents = []
        
        for term in go_terms:
            term_data = {
                'name': term.get('name', 'N/A'),
                'accession': term.get('accession', 'N/A'),
                'namespace': term.get('namespace', 'N/A'),
                'definition': term.get('definition', 'No definition available'),
                'synonyms': term.get('synonyms', ['No synonyms available']),
                'parents': term.get('parents', [])
            }
            terms.append(term_data)
            all_parents.extend(term['parents'])

        return render_template('ontology.html', 
                             gene_name=gene_name,
                             terms=terms,
                             all_parents=all_parents)
    
    return render_template('ontology.html', 
                         gene_name=gene_name,
                         terms=None,
                         all_parents=None)

def extract_gene_name(gene_name):
    match = re.search(r'\((.*?)\)', gene_name)
    return match.group(1).strip() if match else gene_name

def get_gene_ontology_terms(gene_name):
    url = f"https://rest.ensembl.org/ontology/name/{quote_plus(gene_name)}?content-type=application/json"
    response = requests.get(url)
    return response.json() if response.status_code == 200 else None

@app.route('/select_population', methods=['POST'])
def select_population():
    selected_populations = request.form.getlist('selected_populations')
    population_info = get_population_info(selected_populations)
    return render_template('population_info.html', population_info=population_info)

def get_population_info(populations):
    population_map = {
        'haryana/northern india': "Haryana is a northern state of India...",
        'tamilnadu/southern india': "Tamil Nadu, India..."
    }
    return {pop: population_map[pop] for pop in populations if pop in population_map}

@app.errorhandler(400)
def bad_request(error):
    return render_template('error.html',
                         error_message=error.description,
                         error_code=400), 400