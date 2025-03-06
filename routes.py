from flask import render_template, request, jsonify, Response, session
from models import SNP, TajimaD, CLRTest, FstSNP
from functions import (get_snp_info, get_gene_ontology_terms, get_t2d_snps, 
                       get_gene_coordinates_ensembl, get_tajima_d_data,
                       get_clr_data, download_clr_data, download_tajima_d_data,
                       download_fst_data, get_tajima_d_snp_results, get_clr_snp_results,
                       get_fst_data, process_snp_results)
from main import app
import json
import numpy as np
import logging
import uuid
from store import analysis_store
# Configure the logger for debugging
logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)


@app.route('/', methods=['GET', 'POST'])
def home():
    """
    Creates the home page, allowing users to search for SNPs by:
     - SNP ID
     - Chromosome
     - Genomic position range
     - Gene name

    For POST requests, it retrieves SNPs, queries associated Tajima's D and CLR data via helper functions,
    and then adds FST data for each SNP as in the original functionality.
    Returns snp_results.html for POST requests and index.html for GET requests.
    """
    if request.method == 'POST':
        # Extract query parameters from the form
        snp_id = request.form.get('snp_name')
        chromosome = request.form.get('chromosome')
        start = request.form.get('start')
        end = request.form.get('end')
        gene_name = request.form.get('gene_name')

        # Retrieve SNPs based on user input
        snps = get_snp_info(snp_id, chromosome, start, end, gene_name)
        if not snps:
            return render_template('snp_results.html', snp_info=[])

        # Extract positions and unique chromosomes from SNPs for downstream queries
        snp_positions = [snp['grch38_start'] for snp in snps]
        snp_chromosomes = {snp['chromosome'] for snp in snps}

        # Retrieve Tajima's D data (structured by population) using a helper function
        tajima_d_dict = get_tajima_d_snp_results(snp_chromosomes, snp_positions)
        # Retrieve CLR data using a helper function
        clr_dict = get_clr_snp_results(snp_chromosomes)

        # Process the SNP results for rendering using the new functions
        # This function returns a list of SNP dictionaries with 'positive_selection' entries
        # that contain Tajima's D, CLR, and Alpha values.
        snp_info = process_snp_results(snps, tajima_d_dict, clr_dict)

        # Preserve the original FST functionality by looping through each SNP result.
        # For each SNP, query the FstSNP table by snp_id and add the FST value for each population.
        for snp in snp_info:
            fst_data = FstSNP.query.filter_by(snp_id=snp['snp_id']).first()
            positive_selection_list = []
            for selection in snp.get("positive_selection", []):
                pop = selection.get("population")
                fst_value = None
                if fst_data:
                    # Use the population code (lowercase) to retrieve the corresponding FST value.
                    fst_value = getattr(fst_data, f'fst_{pop.lower()}', None)
                    fst_value = round(fst_value, 4) if fst_value is not None else None
                selection["fst"] = fst_value if fst_value is not None else "N/A"
                positive_selection_list.append(selection)
            # Update the SNP's positive_selection with the FST data
            snp["positive_selection"] = positive_selection_list

        return render_template('snp_results.html', snp_info=snp_info)

    # GET request returns the search page
    return render_template('index.html')



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

    # Fetch FST data for the selected chromosome and populations
    fst_data, fst_summary_stats = get_fst_data(selected_chromosome, None, selected_populations)

    return render_template(
        'population_analysis.html',
        selected_population_info=selected_population_info,
        tajima_d_data=json.dumps(tajima_d_data, indent=2),
        t2d_snp_data=json.dumps(t2d_snp_data, indent=2),
        clr_data=json.dumps(clr_data, indent=2),
        fst_data=json.dumps(fst_data, indent=2),
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
    user_start = request.args.get("start", type=int)
    user_end = request.args.get("end", type=int)
    gene_name = request.args.get("gene_name", type=str)
    user_selected_chromosome = request.args.get("chromosome", type=str)
    selected_populations = request.args.getlist("selected_population")

    selected_chromosome = user_selected_chromosome

    if gene_name:
        gene_info = get_gene_coordinates_ensembl(gene_name, chromosome=user_selected_chromosome)
        
        if gene_info and gene_info.get("start") and gene_info.get("end"):
            gene_start = gene_info["start"]
            gene_end = gene_info["end"]

            # Validate chromosome match
            if gene_info["chromosome"] != user_selected_chromosome:
                return jsonify({
                    "error": f"Gene '{gene_name}' is not located on chromosome {user_selected_chromosome}."
                }), 400

            # If the user provides a valid sub-region, use it; otherwise, use full gene coordinates
            if user_start and user_end:
                if gene_start <= user_start <= gene_end and gene_start <= user_end <= gene_end:
                    start = user_start
                    end = user_end
                else:
                    return jsonify({
                        "error": f"Selected region ({user_start}-{user_end}) is outside the gene boundaries ({gene_start}-{gene_end})."
                    }), 400
            else:
                start = gene_start
                end = gene_end

        else:
            return jsonify({"error": f"Gene '{gene_name}' not found in Ensembl."}), 400

    else:
        # If no gene is provided, use the user-selected chromosome as the region's chromosome
        start = user_start
        end = user_end
        

    if start is None or end is None:
        return jsonify({"error": "Invalid region parameters. Provide a gene name or valid start/end coordinates."}), 400

    # Fetch data
    tajima_d_data, summary_stats = get_tajima_d_data(selected_chromosome, (start, end), selected_populations)
    t2d_snp_data = get_t2d_snps(selected_chromosome, start, end)
    clr_data, clr_summary_stats = get_clr_data(selected_chromosome, (start, end), selected_populations)
    fst_data, fst_summary_stats = get_fst_data(selected_chromosome, (start, end), selected_populations)

    analysis_id = str(uuid.uuid4())
    analysis_store[analysis_id] = {
        "tajima_d_data": tajima_d_data,
        "summary_stats": summary_stats,
        "clr_data": clr_data,
        "clr_summary_stats": clr_summary_stats,
        "fst_data": fst_data,
        "fst_summary_stats": fst_summary_stats,
        "selected_chromosome": selected_chromosome,
        "start": start,
        "end": end,
        "selected_populations": selected_populations
    }

    print(f"DEBUG - Stored Analysis ID: {analysis_id}")

    return jsonify({
        "analysis_id": analysis_id,
        "tajima_d_data": tajima_d_data,
        "t2d_snp_data": t2d_snp_data,
        "clr_data": clr_data,
        "summary_stats": summary_stats,
        "clr_summary_stats": clr_summary_stats,
        "fst_data": fst_data,
        "fst_summary_stats": fst_summary_stats,
        "selected_chromosome": selected_chromosome,
        "start": start,
        "end": end,
        "selected_populations" : selected_populations
    })
    
@app.route('/download_fst', methods=['GET'])
def download_fst():
    """
    Endpoint to download FST data as a text file.
    """
    return download_fst_data()


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
    Fetches Gene Ontology (GO) terms and Ensembl gene information for a given gene.
    """

    chromosome = request.args.get("chromosome")  # Get chromosome from URL parameters
    position = request.args.get("position", type=int)
    rsID = request.args.get("rsID")

    ensembl_info = get_gene_coordinates_ensembl(gene_name, rsID=rsID, chromosome=chromosome, position=position)

    if not ensembl_info or "ensembl_id" not in ensembl_info:
        return jsonify({"error": f"Gene '{gene_name}' not found in Ensembl"}), 400


    # Fetch GO Terms
    go_terms = get_gene_ontology_terms(gene_name)

    return render_template(
        "ontology.html",
        gene_name=gene_name,
        ensembl_info=ensembl_info,  # Use ensembl_info directly
        go_terms=go_terms if go_terms else {}
    )


@app.route('/population')
def population():
    # Query your Population table and pass data to a template
    
    return render_template('population_info.html')



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

