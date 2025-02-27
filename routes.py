from flask import render_template, request, jsonify, Response, session
from models import SNP, TajimaD, CLRTest
from functions import (get_snp_info, get_gene_ontology_terms, get_t2d_snps, 
                       get_gene_coordinates_ensembl, get_tajima_d_data,
                       get_clr_data, download_clr_data, download_tajima_d_data)
from main import app
import json
import numpy as np
import logging

# Configure the logger
logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)


@app.route('/', methods=['GET', 'POST'])
def home():
    if request.method == 'POST':
        snp_id = request.form.get('snp_name')
        chromosome = request.form.get('chromosome')
        start = request.form.get('start')
        end = request.form.get('end')
        gene_name = request.form.get('gene_name')

        logger.debug(f"Received POST request with parameters: snp_id={snp_id}, chromosome={chromosome}, start={start}, end={end}, gene_name={gene_name}")

        snps = get_snp_info(snp_id, chromosome, start, end, gene_name)
        logger.debug(f"SNPs fetched: {snps}")

        snp_info = []
        if snps:
            for snp in snps:
                logger.debug(f"Processing SNP: {snp}")

                tajima_d_results = (
                    TajimaD.query
                    .filter(TajimaD.chromosome == snp['chromosome'])
                    .filter(TajimaD.bin_start <= snp['grch38_start'], TajimaD.bin_end >= snp['grch38_start'])
                    .all()
                )
                logger.debug(f"Tajima's D results: {tajima_d_results}")

                clr_results = (
                    CLRTest.query
                    .filter(CLRTest.chromosome == snp['chromosome'])
                    .all()
                )
                logger.debug(f"CLR results: {clr_results}")

                positive_selection = {}

                for tajima in tajima_d_results:
                    if tajima.population not in positive_selection:
                        positive_selection[tajima.population] = {}
                    positive_selection[tajima.population]["tajima_d"] = tajima.tajima_d

                closest_clr_dict = {}
                for clr in clr_results:
                    if clr.population not in closest_clr_dict or abs(clr.position - snp['grch38_start']) < abs(closest_clr_dict[clr.population].position - snp['grch38_start']):
                        closest_clr_dict[clr.population] = clr

                for population, closest_clr in closest_clr_dict.items():
                    if population not in positive_selection:
                        positive_selection[population] = {}
                    positive_selection[population]["clr"] = closest_clr.clr
                    positive_selection[population]["alpha"] = closest_clr.alpha

                positive_selection_list = [
                    {
                        "population": pop,
                        "tajima_d": data.get("tajima_d", "N/A"),
                        "clr": data.get("clr", "N/A"),
                        "alpha": data.get("alpha", "N/A"),
                    }
                    for pop, data in positive_selection.items()
                ]

                snp_info.append({
                    "snp_id": snp['snp_id'],
                    "chromosome": snp['chromosome'],
                    "grch38_start": snp['grch38_start'],
                    "gene_name": snp['gene_name'],
                    "p_value": snp['p_value'],
                    "reference_allele": snp['reference_allele'],
                    "alternative_allele": snp['alternative_allele'],
                    "positive_selection": positive_selection_list
                })

        logger.debug(f"SNP info to be rendered: {snp_info}")
        return render_template('snp_results.html', snp_info=snp_info)

    return render_template('index.html')


# @app.route('/', methods=['GET', 'POST'])
# def home():
#     """
#     Creates the home page, allowing users to search for SNPs by:
#      - SNP ID
#      - Chromosome
#      - genomic position range
#      - Gene name

#     Handles SNP search and fetches associated Tajima's D and CLR data.

#     Returns:
#         - index.html for GET requests (search form).
#         - snp_results.html for POST requests (search results).
#     """
#     if request.method == 'POST':
#         # Extract search parameters
#         snp_id = request.form.get('snp_name')
#         chromosome = request.form.get('chromosome')
#         start = request.form.get('start')
#         end = request.form.get('end')   
#         gene_name = request.form.get('gene_name')
      
#         # Fetch SNPs
#         snps = get_snp_info(snp_id, chromosome, start, end, gene_name)

#         snp_info = []
#         if snps:
#             for snp in snps:
#                 # Fetch the closest Tajima's D for each population where SNP falls within the bin
#                 tajima_d_results = (
#                     TajimaD.query
#                     .filter(TajimaD.chromosome == snp['chromosome'])
#                     .filter(TajimaD.bin_start <= snp['grch38_start'], TajimaD.bin_end >= snp['grch38_start'])
#                     .all()
#                 )

#             # Fetch all CLR results for the chromosome
#             clr_results = (
#                 CLRTest.query
#                 .filter(CLRTest.chromosome == snp['chromosome'])
#                 .all()
#             )

#             positive_selection = {}

#             # Store Tajima’s D values
#             for tajima in tajima_d_results:
#                 if tajima.population not in positive_selection:
#                     positive_selection[tajima.population] = {}
#                 positive_selection[tajima.population]["tajima_d"] = tajima.tajima_d

#             # Find the closest CLR position to query SNP
#             closest_clr_dict = {}
#             for clr in clr_results:
#                 if clr.population not in closest_clr_dict or abs(clr.position - snp['grch38_start']) < abs(closest_clr_dict[clr.population].position - snp['grch38_start']):
#                     closest_clr_dict[clr.population] = clr

#             # Store CLR and Alpha values from the closest CLR entry
#             for population, closest_clr in closest_clr_dict.items():
#                 if population not in positive_selection:
#                     positive_selection[population] = {}
#                 positive_selection[population]["clr"] = closest_clr.clr
#                 positive_selection[population]["alpha"] = closest_clr.alpha

#             # Convert dictionary to list for rendering
#             positive_selection_list = [
#                 {
#                     "population": pop,
#                     "tajima_d": data.get("tajima_d", "N/A"),
#                     "clr": data.get("clr", "N/A"),
#                     "alpha": data.get("alpha", "N/A"),
#                 }
#                 for pop, data in positive_selection.items()
#             ]

#             snp_info.append({
#                 "snp_id": snp['snp_id'],
#                 "chromosome": snp['chromosome'],
#                 "grch38_start": snp['grch38_start'],
#                 "gene_name": snp['gene_name'],
#                 "p_value": snp['p_value'],
#                 "reference_allele": snp['reference_allele'],
#                 "alternative_allele": snp['alternative_allele'],
#                 "positive_selection": positive_selection_list
#             })

#         return render_template('snp_results.html', snp_info=snp_info)

#     return render_template('index.html')




@app.route('/population_analysis', methods=['GET', 'POST'])
def population_analysis():
    """
    Allows analysis of Tajima's D and CLR statistics for selected populations.

    GET: Loads the population analysis form.
    POST: Fetches Tajima's D, CLR, and SNP data for selected populations & chromosomes.

    Returns:
        - population_analysis.html with the requested population statistics.
    """
    if request.method == 'GET':
        return render_template(
            'population_analysis.html', 
            selected_population_info={}, # No population selected yet
            tajima_d_data=json.dumps({}), # Empty JSON data for initial rendering 
            t2d_snp_data=json.dumps({}),
            clr_data=json.dumps({}),
            selected_chromosome="10" # Default chromosome if not selected
        )

    selected_populations = request.form.getlist('selected_population')
    selected_chromosome = request.form.get('selected_chromosome', '10')

    # Dictionary of population descriptions
    population_descriptions = {
        "BEB": "Bengali population from Bangladesh",
        "GIH": "Gujarati Indian population from Houston, USA",
        "ITU": "Indian Telugu population from the UK",
        "PJL": "Punjabi population from Lahore, Pakistan",
        "STU": "Sri Lankan Tamil population from the UK"
    }

    # Ensure selected_population_info is a dictionary
    selected_population_info = {
        pop: population_descriptions.get(pop, "Unknown population") for pop in selected_populations
    }

    # Fetch SNPs associated with Type 2 Diabetes (T2D) for the selected chromosome
    t2d_snp_data = get_t2d_snps(selected_chromosome)

    #Fetch Tajima's D data for the selected chromosome and populations
    tajima_d_data, _ = get_tajima_d_data(selected_chromosome, None, selected_populations)

    # Fetch CLR data for the selected chromosome and populations
    clr_data, _ = get_clr_data(selected_chromosome, None, selected_populations)

    return render_template(
            'population_analysis.html',
            tajima_d_data=json.dumps(tajima_d_data, indent=2),  # Ensure JSON is correctly formatted
            t2d_snp_data=json.dumps(t2d_snp_data, indent=2),
            clr_data=json.dumps(clr_data, indent=2),
            selected_population_info=selected_population_info,
            selected_chromosome=selected_chromosome
        )



@app.route('/population_analysis_region', methods=['GET'])
def population_analysis_region():
    """
    Fetches Tajima's D and CLR statistics for a user-defined genomic region and returns SNP data.
    Query Parameters:
        - chromosome (str): Chromosome number (default: "10").
        - start (int): Start position of the region.
        - end (int): End position of the region.
        - gene_name (str): Gene name to fetch coordinates from Ensembl.
        - selected_population (list): List of selected populations.
    Returns:
        - JSON containing Tajima's D, CLR, and T2D SNP data for the requested region.
    """
    # Retrieve query parameters
    start = request.args.get("start", type=int)
    end = request.args.get("end", type=int)
    gene_name = request.args.get("gene_name", type=str)
    chromosome = request.args.get("chromosome", "10")
    selected_populations = request.args.getlist("selected_population")

     # If a gene name is provided, fetch its genomic coordinates from Ensembl
    if gene_name:
        gene_info = get_gene_coordinates_ensembl(gene_name)
        if gene_info and gene_info.get("start") and gene_info.get("end"):  # Ensure values exist
            chromosome, start, end = gene_info["chromosome"], gene_info["start"], gene_info["end"]
        else:
            return jsonify({"error": f"Gene '{gene_name}' not found in Ensembl"}), 400

    if start is None or end is None:
        return jsonify({"error": "Invalid region parameters"}), 400

    # Fetch data
    tajima_d_data, summary_stats = get_tajima_d_data(chromosome, (start, end), selected_populations)
    t2d_snp_data = get_t2d_snps(chromosome, start, end)
    clr_data, clr_summary_stats = get_clr_data(chromosome, (start, end), selected_populations)

    return jsonify({
        "tajima_d_data": tajima_d_data,
        "t2d_snp_data": t2d_snp_data,
        "clr_data": clr_data,
        "summary_stats": summary_stats,
        "clr_summary_stats": clr_summary_stats
    })



@app.route('/download_tajima_d', methods=['GET'])
def download_tajima_d():
    """
    Generates a downloadable text file containing Tajima's D values.

    Query Parameters:
        - chromosome: Chromosome number.
        - start: Start position (optional).
        - end: End position (optional).
        - selected_population: Population(s) to filter by.

    Returns:
        - Generates .txt file with Tajima's D values, average mean and standard deviation for a requested genomic region.
    """
    return download_tajima_d_data()



@app.route('/download_clr', methods=['GET'])
def download_clr():
    """
    Generates a downloadable text file containing CLR (Composite Likelihood Ratio) values.
    
    Query Parameters:
        - chromosome (str): Chromosome number.
        - start (int, optional): Start position of the region.
        - end (int, optional): End position of the region.
        - selected_population (list): List of selected populations.

    Returns:
        - A text file containing CLR values for the specified genomic region.
    """

    return download_clr_data()



@app.route('/about')
def about():
    """
    Renders the "About" page.

    Returns:
        - about.html
    """
    return render_template('about.html')



@app.route('/gene_terms/<gene_name>')
def gene_terms(gene_name):
  """
    Fetches Gene Ontology (GO) terms for a given gene.

    Parameters:
        - gene_name: Gene name.

    Returns:
        - ontology.html with GO terms.
    """
  go_terms = get_gene_ontology_terms(gene_name)
  return render_template("ontology.html", gene_name=gene_name, go_terms=go_terms if go_terms else {})


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

# tajima_bp = Blueprint("tajima", __name__)

# @tajima_bp.route("/tajima_d", methods=["GET"])
# def get_tajima_d():
#     """
#     Retrieve Tajima's D results based on query parameters.
#     Example: /api/tajima_d?population=BEB&chromosome=10&start=50000&end=1500000
#     """
#     population = request.args.get("population", type=str)
#     chromosome = request.args.get("chromosome", type=str)
#     start = request.args.get("start", type=int)
#     end = request.args.get("end", type=int)

#     query = TajimaD.query

#     if population:
#         query = query.filter(TajimaD.population == population)
#     if chromosome:
#         query = query.filter(TajimaD.chromosome == chromosome)
#     if start and end:
#         query = query.filter(TajimaD.bin_start >= start, TajimaD.bin_end <= end)

#     results = query.all()

#     return jsonify([
#         {
#             "population": row.population,
#             "chromosome": row.chromosome,
#             "bin_start": row.bin_start,
#             "bin_end": row.bin_end,
#             "n_snps": row.n_snps,
#             "tajima_d": row.tajima_d
#         }
#         for row in results
#     ])

