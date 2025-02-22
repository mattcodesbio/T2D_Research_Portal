from flask import render_template, request, jsonify, session
from models import SNP, TajimaD
from functions import get_snp_info, get_gene_ontology_terms, get_gene_coordinates_ensembl
from main import app
import json
import pandas as pd

@app.route('/', methods=['GET', 'POST'])
def home():
    if request.method == 'POST':
        snp_id = request.form.get('snp_name')
        chromosome = request.form.get('chromosome')
        start = request.form.get('start')
        end = request.form.get('end')   
        gene_name = request.form.get('gene_name')

        query = SNP.query

        if snp_id:
            query = query.filter(SNP.snp_id.ilike(f"%{snp_id}%"))
        if chromosome:
            query = query.filter(SNP.chromosome == str(chromosome))
        if start:
            query = query.filter(SNP.grch38_start >= int(start))
        if end:
            query = query.filter(SNP.grch38_start <= int(end))
        if gene_name:
            query = query.filter(SNP.gene_name.ilike(f"%{gene_name}%"))

        snps = query.all()

        snp_info = []
        for snp in snps:
            # Fetch the closest Tajima's D for each population where SNP falls within the bin
            tajima_d_results = (
                TajimaD.query
                .filter(TajimaD.chromosome == snp.chromosome)
                .filter(TajimaD.bin_start <= snp.grch38_start, TajimaD.bin_end >= snp.grch38_start)
                .all()
            )

            positive_selection = []
            for tajima in tajima_d_results:
                positive_selection.append({
                    "population": tajima.population,
                    "tajima_d": tajima.tajima_d
                })

            snp_info.append({
                "snp_id": snp.snp_id,
                "chromosome": snp.chromosome,
                "grch38_start": snp.grch38_start,
                "gene_name": snp.gene_name,
                "p_value": snp.p_value,
                "reference_allele": snp.reference_allele,
                "alternative_allele": snp.alternative_allele,
                "positive_selection": positive_selection
            })

        return render_template('snp_results.html', snp_info=snp_info)

    return render_template('index.html')



@app.route('/population_analysis', methods=['GET', 'POST'])
def population_analysis():
    """
    Handles whole chromosome Tajima's D analysis.
    - GET: Loads the page with default values.
    - POST: Handles user selection for specific populations and chromosome.
    """
    if request.method == 'GET':  
        # Default: Load the page with no filters
        return render_template(
            'population_analysis.html',
            selected_population_info={},
            tajima_d_data=json.dumps({}),
            t2d_snp_data=json.dumps({}),
            selected_chromosome="10"  # Default chromosome
        )

    # Handle POST request for user selection
    selected_populations = request.form.getlist('selected_population')
    selected_chromosome = request.form.get('selected_chromosome', '10')

    # Population descriptions
    population_descriptions = {
        "BEB": "Bengali population from Bangladesh",
        "GIH": "Gujarati Indian population from Houston, USA",
        "ITU": "Indian Telugu population from the UK",
        "PJL": "Punjabi population from Lahore, Pakistan",
        "STU": "Sri Lankan Tamil population from the UK"
    }

    # Get selected population descriptions
    selected_population_info = {pop: population_descriptions.get(pop, "Unknown population") for pop in selected_populations}

    # Fetch Tajima’s D values for the selected chromosome and populations
    tajima_d_results = (
        TajimaD.query
        .filter(TajimaD.population.in_(selected_populations), TajimaD.chromosome == selected_chromosome)
        .all()
    )


    # Organize Tajima's D data per population
    tajima_d_data = {}
    for pop in selected_populations:
        tajima_d_data[pop] = [
            {"bin_start": row.bin_start, "bin_end": row.bin_end, "tajima_d": row.tajima_d}
            for row in tajima_d_results if row.population == pop
        ]

    # Fetch T2D-associated SNPs for the chromosome
    t2d_snps = (
        SNP.query
        .filter(SNP.chromosome == selected_chromosome)
        .filter(SNP.p_value < 5e-8)  # Significant T2D SNPs
        .all()
    )

    # Store SNPs in a dictionary
    t2d_snp_data = [{"snp_id": snp.snp_id, "position": snp.grch38_start} for snp in t2d_snps]

    return render_template(
        'population_analysis.html',
        selected_population_info=selected_population_info,
        tajima_d_data=json.dumps(tajima_d_data),
        t2d_snp_data=json.dumps(t2d_snp_data),
        selected_chromosome=selected_chromosome
    )


@app.route('/population_analysis_region', methods=['GET'])
def population_analysis_region():
    """
    Fetches Tajima's D for a user-defined region and T2D SNPs.
    """
    start = request.args.get("start", type=int)
    end = request.args.get("end", type=int)
    gene_name = request.args.get("gene_name", type=str)
    chromosome = request.args.get("chromosome", "10")
    selected_populations = request.args.getlist("selected_population")

    # If a gene name is provided, get the region from Ensembl
    if gene_name:
        gene_info = get_gene_coordinates_ensembl(gene_name)
        if gene_info:
            chromosome, start, end = gene_info["chromosome"], gene_info["start"], gene_info["end"]
        else:
            return jsonify({"error": f"Gene '{gene_name}' not found in Ensembl"}), 400

    # Ensure valid region
    if start is None or end is None:
        return jsonify({"error": "Invalid region parameters"}), 400

    # Fetch Tajima's D values for the selected populations in the region
    tajima_d_results = (
        TajimaD.query
        .filter(TajimaD.chromosome == chromosome)
        .filter(TajimaD.bin_start >= start, TajimaD.bin_end <= end)
        .filter(TajimaD.population.in_(selected_populations))  # Only return selected populations
        .all()
    )

    region_tajima_d = {pop: [] for pop in selected_populations}
    for tajima in tajima_d_results:
        region_tajima_d[tajima.population].append({
            "bin_start": tajima.bin_start,
            "bin_end": tajima.bin_end,
            "tajima_d": tajima.tajima_d
        })

    # Fetch T2D SNPs in the region
    t2d_snps = (
        SNP.query
        .filter(SNP.chromosome == chromosome)
        .filter(SNP.grch38_start >= start, SNP.grch38_start <= end)
        .filter(SNP.p_value < 5e-8)  # Ensure it's a significant T2D SNP
        .all()
    )

    # Map SNPs to their nearest Tajima's D bin
    t2d_snp_data = []
    for snp in t2d_snps:
        closest_tajima_d = (
            TajimaD.query
            .filter(TajimaD.chromosome == chromosome)
            .filter(TajimaD.bin_start <= snp.grch38_start, TajimaD.bin_end >= snp.grch38_start)
            .first()  # Get the closest match
        )

        t2d_snp_data.append({
            "snp_id": snp.snp_id,
            "position": snp.grch38_start,
            "tajima_d": closest_tajima_d.tajima_d if closest_tajima_d else None
        })

    return jsonify({
        "tajima_d_data": region_tajima_d,
        "t2d_snp_data": t2d_snp_data
    })


@app.route('/about')
def about():
    return render_template('about.html')

@app.route('/gene_terms/<gene_name>')
def gene_terms(gene_name):
    go_terms = get_gene_ontology_terms(gene_name)
    
    if not go_terms:
        return render_template("ontology.html", gene_name=gene_name, go_terms={})

    return render_template("ontology.html", gene_name=gene_name, go_terms=go_terms)


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

