import pandas as pd
from models import db, SNP, TajimaD, CLRTest, FstSNP
import requests
import numpy as np
import os
from flask import jsonify, Response, request
from sqlalchemy import func, text
import bisect
def get_snp_info(snp_id=None, chromosome=None, start=None, end=None, gene_name=None):
    """
    Retrieves SNP information from the SNP database based on the user's query.

    Args:
        snp_id (str): SNP ID (e.g., 'rs12345').
        chromosome (str): Chromosome (e.g., '1').
        start (int): Start position on the chromosome.
        end (int): End position on the chromosome.
        gene_name (str): Mapped Gene name.

    Returns:
        list: A list of dictionaries, each containing retrieved information from the SNP database about a SNP.
              Returns None if no SNPs are found.
    """
    query = SNP.query # start a query on the SNP table

        
    # Apply filters based on the user's input    
    if snp_id:
        query = query.filter(SNP.snp_id.ilike(f"%{snp_id}%")) # Case-insensitive search for SNP ID
    if chromosome:
        query = query.filter(SNP.chromosome == str(chromosome)) # query chromosome with matching SNP chromosome
    if start:
        query = query.filter(SNP.grch38_start >= int(start)) # Filter SNPs starting from this position
    if end:
        query = query.filter(SNP.grch38_start <= int(end)) # Filter SNPs up to this position
    if gene_name:
        query = query.filter(SNP.gene_name.ilike(f"%{gene_name}%")) # Search for SNPs by gene name

    snps = query.all() # carry out the query

    # Convert the retrieved database objects into a dictionary format
    results = []
    for snp in snps:
        results.append({
            'snp_id': snp.snp_id,
            'chromosome': snp.chromosome,
            'grch38_start': snp.grch38_start,
            'gene_name': snp.gene_name,
            'p_value': snp.p_value,
            'reference_allele': snp.reference_allele,
            'alternative_allele': snp.alternative_allele,
            'consequence': snp.consequence
        })
    
    return results if results else None # Return results if found, else return None



def get_tajima_d_snp_results(snp_chromosomes, snp_positions):
    """
    Retrieve Tajima's D results based on the user's SNP query in the home route.

    This function queries the database for Tajima's D statistics within the specified 
    SNP chromosome(s) and position range.

    Parameters:
    snp_chromosomes: A list of str containing the chromosome names where the SNPs are located.
    snp_positions : A list of str containing the SNP positions on the chromosomes.

    Returns:
        A nested dictionary where:
        - Keys are tuples of (chromosome, bin_start, bin_end).
        - Values are dictionaries mapping population names to Tajima's D values.
    """

    # Query the database for Tajima's D statistics based on the SNP chromosome and position range
    results = (
        TajimaD.query
        .filter(TajimaD.chromosome.in_(snp_chromosomes))
        .filter(TajimaD.bin_start <= max(snp_positions))
        .filter(TajimaD.bin_end >= min(snp_positions))
        .all()
    )

    tajima_d_dict = {}
    for tajima in results:
        key = (tajima.chromosome, tajima.bin_start, tajima.bin_end)
        if key not in tajima_d_dict:
                tajima_d_dict[key] = {}
        tajima_d_dict.setdefault(key, {})[tajima.population] = tajima.tajima_d

    return tajima_d_dict


def get_clr_snp_results(snp_chromosomes):
    """
    Retrieve CLR results based on the user's SNP query in the home route.

    This function queries the database for CLR statistics based on SNP chromosome input.
    It efficiently retrieves results by filtering for a range of chromosome values.

    Parameters:
    snp_chromosomes: A list of chromosome names or numbers for which CLR results are needed.

    SQL Query:
    - Optimised to retrieve CLR values from the indexed 'clr_results' table 
        (CREATE INDEX idx_clr_chrom_pop_pos ON clr_results (chromosome, population, position))
    - Uses SQL parameter binding to prevent injections. 

    Returns:
        A nested dictionary where:
        - Keys are tuples of (chromosome, population).
        - Values are dictionaries mapping population names to CLR and alpha values.
    """
    
    # SQL query to retrieve relevant CLR data from the database
    query = text("""
        SELECT chromosome, population, position, clr, alpha
        FROM clr_results
        WHERE chromosome BETWEEN :chrom_min AND :chrom_max
    """)

    # Query parameters to filter chromosome range dynamically
    params = {"chrom_min": min(snp_chromosomes), "chrom_max": max(snp_chromosomes)}
    clr_results = db.session.execute(query, params).fetchall()

    clr_dict = {}
    for row in clr_results:
            chrom, pop, pos, clr_value, alpha = row
            if (chrom, pop) not in clr_dict:
                clr_dict[(chrom, pop)] = {}
            clr_dict[(chrom, pop)][pos] = {"clr": clr_value, "alpha": alpha}
    return clr_dict


def get_clr_closest_position(positions, snp_position):
    """
    Helper function for process_snp_results to efficiently search for the closest CLR position using binary search. 

    Args:
        positions (list[int]): A sorted list of CLR positions.
        snp_position (int): The SNP position to find the closest CLR position for.

    Returns:
        int: The closest CLR position to `snp_position`.
    """
    # Find the insertion point for snp_position
    idx = bisect.bisect_left(positions, snp_position)

    # If snp_position is smaller than all positions, return the first position
    if idx == 0:
        return positions[0]
    # If snp_position is larger than all positions, return the last position
    if idx == len(positions):
        return positions[-1]

    # Get positions around the insertion point
    left, right = positions[idx - 1], positions[idx]
    return left if snp_position - left <= right - snp_position else right #faster for chr10 search ** pick one and delete
    # return left if abs(left - snp_position) <= abs(right - snp_position) else right #faster for tcf7l2 search ** pick one and delete


def process_snp_results(snps, tajima_d_dict, clr_dict):
    """
    Process SNPs and associate them with Tajima's D and CLR results in the home route.
    
    Args:
        snps: A list of SNP dictionaries containing SNP information.
        tajima_d_dict: A dictionary containing Tajima's D results.
        clr_dict: A dictionary containing CLR results.

    Returns:
        list: A list of dictionaries aggregating all the information for each SNP.
    """
    
    populations = ['BEB', 'GIH', 'ITU', 'PJL', 'STU']
    snp_info = []

    for snp in snps:
        snp_position, snp_chrom = snp['grch38_start'], snp['chromosome']

        # Fetch Tajima's D for this SNP
        positive_selection = {
            pop: {"tajima_d": tajima_d}
            for (chrom, bin_start, bin_end), pop_data in tajima_d_dict.items()
            if chrom == snp_chrom and bin_start <= snp_position <= bin_end
            for pop, tajima_d in pop_data.items()
        }

        # Fetch CLR for this SNP
        for pop in populations:
            if (snp_chrom, pop) in clr_dict:
                positions = sorted(clr_dict[(snp_chrom, pop)].keys())
                if positions:
                    closest_pos = get_clr_closest_position(positions, snp_position)
                    if closest_pos is not None:
                        positive_selection.setdefault(pop, {}).update(clr_dict[(snp_chrom, pop)][closest_pos])
        
        # Convert to list format for rendering
        positive_selection_list = [
            {"population": pop, **{k: data.get(k, "N/A") for k in ["tajima_d", "clr", "alpha"]}}
            for pop, data in positive_selection.items()
        ]

        # Store SNP results
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

    return snp_info




def get_tajima_d_data(chromosome, region=None, populations=None):
    """
    Fetch Tajima's D statistics for a whole chromosome or a region.

    Args:
        chromosome (str): Chromosome to filter by.
        region (tuple, optional argument): Start and end positions of the region. Defaults to None.
        populations (list, optional argument): List of populations to filter by. Defaults to None.

    Returns:
        tuple: A tuple containing two dictionaries:
            - tajima_d_data: Tajima's D values for each bin and respective population.
            - summary_stats: Mean and standard deviation of Tajima's D across the genomic region requested for each population.
    """
    query = TajimaD.query.filter(TajimaD.chromosome == chromosome) # Query Tajima's D table, filters for the chromosome specified

    if populations:
        query = query.filter(TajimaD.population.in_(populations)) # Filter by selected populations

    if region:
        start, end = region
        query = query.filter(TajimaD.bin_start >= start, TajimaD.bin_end <= end) # Filter by genomic region

    # Initialise a dictionary with all selected populations
    tajima_d_data = {pop: [] for pop in populations} if populations else {}

    for tajima in query.all():
        if tajima.population not in tajima_d_data:
            tajima_d_data[tajima.population] = []
        tajima_d_data[tajima.population].append({
            "bin_start": tajima.bin_start,
            "bin_end": tajima.bin_end,
            "tajima_d": tajima.tajima_d
        })
            
    # Calculate mean and standard deviation for each population    
    summary_stats = {
        pop: {
            "mean": round(np.mean([d["tajima_d"] for d in data]), 4) if data else None,
            "std_dev": round(np.std([d["tajima_d"] for d in data]), 4) if data else None
        }
        for pop, data in tajima_d_data.items()
    }

    return tajima_d_data, summary_stats

def get_clr_data(chromosome, region=None, populations=None):
    """
    Fetch CLR statistics for a chromosome or a genomic region.

    Args:
        chromosome (str): Chromosome to filter by.
        region (tuple, optional): Start and end positions of the region. Defaults to None.
        populations (list, optional): List of populations to filter by. Defaults to None.

    Returns:
        tuple: A dictionary of CLR values and summary statistics.
    """
    query = CLRTest.query.filter(CLRTest.chromosome == chromosome) # Query CLRTest table, filters for the chromosome specified

    if populations:
        query = query.filter(CLRTest.population.in_(populations)) # Filter by population

    if region:
        start, end = region
        query = query.filter(CLRTest.position >= start, CLRTest.position <= end) # Filter by genomic region

    # Initialise the dictionary with all selected populations
    clr_data = {pop: [] for pop in populations} if populations else {}

    for clr in query.all():
        if clr.population not in clr_data:
            clr_data[clr.population] = []
        clr_data[clr.population].append({
            "position": clr.position,
            "clr": clr.clr,
            "alpha": clr.alpha
        })
            
    # Calculate summary statistics    
    clr_summary_stats = {
        pop: {
            "mean_clr": round(np.mean([d["clr"] for d in data]), 4) if data else None,
            "std_dev_clr": round(np.std([d["clr"] for d in data]), 4) if data else None,
            "mean_alpha": round(np.mean([d["alpha"] for d in data]), 4) if data else None,
            "std_dev_alpha": round(np.std([d["alpha"] for d in data]), 4) if data else None
        }
        for pop, data in clr_data.items()
    }

    return clr_data, clr_summary_stats


def get_t2d_snps(chromosome, start=None, end=None):
    """
    Fetch significant T2D SNPs within a whole chromosome or requested region.

    Args:
        chromosome (str): Chromosome to filter by.
        start (int, optional argument): Start position of the region. Defaults to None.
        end (int, optional argument): End position of the region. Defaults to None.

    Returns:
        list: A list of dictionaries, each containing the SNP ID and position.
    """
    query = SNP.query.filter(SNP.chromosome == chromosome) # Query SNP table by chromosome

    if start and end:
        query = query.filter(SNP.grch38_start >= start, SNP.grch38_start <= end) # Filter by genomic range

    return [{"snp_id": snp.snp_id, "position": snp.grch38_start} for snp in query.all()] # Convert database results into a list of dictionaries
        


def download_tajima_d_data():
    """
    Generates a text file containing Tajima's D statistics for a region.

    This function retrieves Tajima's D data for a specified chromosome and genomic region,
    optionally filtered by population. It handles requests with either explicit start/end 
    coordinates or a gene name (in which case it fetches gene coordinates from Ensembl). 
    The function generates a plain text file with Tajima's D values and summary statistics 
    (mean and standard deviation) for each selected population.

    Request Parameters:
        chromosome (str): Chromosome (e.g., "10").
        gene_name (str, optional): Gene name to define the region.
        start (int, optional): Start position of the genomic region.
        end (int, optional): End position of the genomic region.
        selected_population (list, optional): List of populations to include.

    Returns:
        flask.Response: A Flask Response object containing the text file with Tajima's D data.
                       If no data is found, returns a response with an error message.
    """
        
    # Extract parameters
    user_selected_chromosome = request.args.get("chromosome")
    gene_name = request.args.get("gene_name")  # Optional gene name
    start = request.args.get("start", type=int)
    end = request.args.get("end", type=int)
    populations = request.args.getlist("selected_population")  # List of selected populations

    # If gene name is provided and start/end not provided, fetch its genomic coordinates
    if gene_name and (start is None or end is None):
        gene_info = get_gene_coordinates_ensembl(gene_name)
        if gene_info:
            # Validate that the gene is on the user-selected chromosome
            if gene_info["chromosome"] != user_selected_chromosome:
                return jsonify({
                    "error": f"Gene '{gene_name}' is not located on the selected chromosome {user_selected_chromosome}"
                }), 400
            # Use gene coordinates
            chromosome = gene_info["chromosome"]
            start = gene_info["start"]
            end = gene_info["end"]
        else:
            return jsonify({"error": f"Gene '{gene_name}' not found in Ensembl"}), 400
    else:
        # No gene provided; use the user-selected chromosome
        chromosome = user_selected_chromosome

    # Ensure start and end are valid numbers
    if start is None or end is None:
        return jsonify({"error": "Invalid region parameters. Must provide a valid start and end position or a gene name."}), 400

    # Fetch Tajima's D data safely
    try:
        tajima_d_data, summary_stats = get_tajima_d_data(chromosome, (start, end), populations)
    except Exception as e:
        return jsonify({"error": f"Error retrieving Tajima's D data: {str(e)}"}), 500

    # Handle case where no data exists
    if not tajima_d_data or all(len(entries) == 0 for entries in tajima_d_data.values()):
        response_text = f"No Tajima's D data found for Chromosome {chromosome}, Region {start}-{end}.\n"
        response = Response(response_text, mimetype="text/plain")
        response.headers["Content-Disposition"] = f"attachment; filename=TajimaD_chr{chromosome}_{start}_{end}.txt"
        return response

    # Generate file content
    file_content = ["Population\tBin Start\tBin End\tTajima's D"]
    for pop, entries in tajima_d_data.items():
        for entry in entries:
            file_content.append(f"{pop}\t{entry['bin_start']}\t{entry['bin_end']}\t{entry['tajima_d']:.4f}")

    # Add summary statistics
    file_content.append("\nSummary Statistics")
    file_content.append("Population\tMean Tajima's D\tStd Dev Tajima's D")
    for pop, stats in summary_stats.items():
        mean_tajima_d = stats["mean"] if stats["mean"] is not None else "N/A"
        std_dev_tajima_d = stats["std_dev"] if stats["std_dev"] is not None else "N/A"
        file_content.append(f"{pop}\t{mean_tajima_d}\t{std_dev_tajima_d}")

    # Return the file response
    response = Response("\n".join(file_content), mimetype="text/plain")
    response.headers["Content-Disposition"] = f"attachment; filename=TajimaD_chr{chromosome}_{start}_{end}.txt"
    return response

def download_clr_data():
    """
    Generates a text file containing CLR (Composite Likelihood Ratio) statistics for a specified genomic region.

    This function retrieves CLR statistics for a given chromosome and region, optionally filtered by population. 
    The user can specify the region using explicit start and end coordinates or by providing a gene name, 
    in which case the function fetches the corresponding coordinates from Ensembl. The function generates 
    a plain text file containing CLR values along with summary statistics (mean and standard deviation) for 
    each selected population.

    Request Parameters:
        chromosome (str): Chromosome (e.g., "10").
        gene_name (str, optional): Gene name to define the region.
        start (int, optional): Start position of the genomic region.
        end (int, optional): End position of the genomic region.
        selected_population (list, optional): List of populations to include.

    Returns:
        flask.Response: A Flask Response object containing the plain text file with CLR data.
                        If an error occurs (e.g., missing parameters, gene not found, data retrieval failure), 
                        returns a JSON error message with an appropriate HTTP status code.
    """
    # Extract parameters
    user_selected_chromosome = request.args.get("chromosome")
    gene_name = request.args.get("gene_name")  # Optional gene name
    start = request.args.get("start", type=int)
    end = request.args.get("end", type=int)
    populations = request.args.getlist("selected_population")

    # If gene name is provided, fetch its genomic coordinates

    # If gene name is provided and start/end not provided, fetch its genomic coordinates
    if gene_name and (start is None or end is None):
        gene_info = get_gene_coordinates_ensembl(gene_name)
        if gene_info:
            # Validate that the gene is on the user-selected chromosome
            if gene_info["chromosome"] != user_selected_chromosome:
                return jsonify({
                    "error": f"Gene '{gene_name}' is not located on the selected chromosome {user_selected_chromosome}"
                }), 400
            # Use gene coordinates
            chromosome = gene_info["chromosome"]
            start = gene_info["start"]
            end = gene_info["end"]
        else:
            return jsonify({"error": f"Gene '{gene_name}' not found in Ensembl"}), 400
    else:
        # No gene provided; use the user-selected chromosome
        chromosome = user_selected_chromosome

    # Fetch CLR data safely
    try:
        clr_data, clr_summary_stats = get_clr_data(chromosome, (start, end), populations)
    except Exception as e:
        return jsonify({"error": f"Error retrieving CLR data: {str(e)}"}), 500

    # Handle case where no data exists
    if not clr_data or all(len(entries) == 0 for entries in clr_data.values()):
        response_text = f"No CLR data found for Chromosome {chromosome}, Region {start}-{end}.\n"
        response = Response(response_text, mimetype="text/plain")
        response.headers["Content-Disposition"] = f"attachment; filename=CLR_chr{chromosome}_{start}_{end}.txt"
        return response

    # Generate file content
    file_content = ["Population\tPosition\tCLR\tAlpha"]
    for pop, entries in clr_data.items():
        for entry in entries:
            file_content.append(f"{pop}\t{entry['position']}\t{entry['clr']:.4f}\t{entry['alpha']:.4f}")

    # Add summary statistics
    file_content.append("\nSummary Statistics")
    file_content.append("Population\tMean CLR\tStd Dev CLR\tMean Alpha\tStd Dev Alpha")
    for pop, stats in clr_summary_stats.items():
        mean_clr = stats["mean_clr"] if stats["mean_clr"] is not None else "N/A"
        std_dev_clr = stats["std_dev_clr"] if stats["std_dev_clr"] is not None else "N/A"
        mean_alpha = stats["mean_alpha"] if stats["mean_alpha"] is not None else "N/A"
        std_dev_alpha = stats["std_dev_alpha"] if stats["std_dev_alpha"] is not None else "N/A"
        file_content.append(f"{pop}\t{mean_clr}\t{std_dev_clr}\t{mean_alpha}\t{std_dev_alpha}")

    # Return the file response
    response = Response("\n".join(file_content), mimetype="text/plain")
    response.headers["Content-Disposition"] = f"attachment; filename=CLR_chr{chromosome}_{start}_{end}.txt"
    return response

def load_snps_from_csv(csv_file):
    """
    Loads SNP data from a CSV file into the database.

    Expected CSV Structure:
    The function assumes that the input CSV file contains the following columns:
    - dbSNP: SNP identifier (e.g., rs123456)
    - chromosome: Chromosome number (e.g., 10)
    - GRCh38_start: SNP genomic start position (GRCh38 assembly)
    - Mapped_Genes: Gene(s) mapped to the SNP (comma-separated list)
    - pValue: P-value from the GWAS study
    - reference: Reference allele (e.g., A)
    - alt: Alternative allele (e.g., T)
    - consequence: Functional consequence of the variant

    Example Usage:
    load_snps_from_csv("path/to/snp_data.csv")

    The function will:
    - Reads the CSV file using Pandas.
    - If an SNP exists, it **updates** the existing record.
    - If an SNP doesn’t exist, it **adds** a new record.
    """
        
    df = pd.read_csv(csv_file) # Load CSV data into a Pandas DataFrame

    with db.session.begin():
        for _, row in df.iterrows():
            # Process the "Mapped_Genes" column for consistency
            mapped_genes = None
            if pd.notna(row["Mapped_Genes"]):
                mapped_genes = row["Mapped_Genes"].strip()  # Remove surrounding spaces
                if mapped_genes.startswith('"') and mapped_genes.endswith('"'):  
                    mapped_genes = mapped_genes[1:-1]  # Remove surrounding quotes
                mapped_genes = mapped_genes.replace(", ", ",")  # fix spacing
                mapped_genes = " ".join(mapped_genes.split(","))  # Convert commas to spaces
                
            # Check if the SNP already exists in the database
            existing_snp = SNP.query.filter_by(snp_id=row["dbSNP"]).first()

            if existing_snp:
                # Update existing SNP record
                existing_snp.chromosome = str(row["chromosome"])
                existing_snp.grch38_start = int(row["GRCh38_start"])
                existing_snp.gene_name = mapped_genes
                existing_snp.p_value = row["pValue"]
                existing_snp.reference_allele = row["reference"]
                existing_snp.alternative_allele = row["alt"]
                existing_snp.consequence = row["consequence"]
            else:
                # Insert new SNP data into the database    
                db.session.add(SNP(
                    snp_id=row["dbSNP"],
                    chromosome=str(row["chromosome"]),
                    grch38_start=int(row["GRCh38_start"]),
                    gene_name=mapped_genes,  # Cleaned up gene names
                    p_value=row["pValue"],
                    reference_allele=row["reference"],
                    alternative_allele=row["alt"],
                    consequence=row["consequence"]
                ))

    db.session.commit() # Commit changes to the database



def get_gene_ontology_terms(gene_name):
    """
    Retrieves GO terms for a given gene name using MyGene.info and QuickGO APIs.

    Args:
        gene_name (str): The name of the gene to query.

    Returns:
        dict or None: A dictionary containing GO terms categorized by namespace 
                     (biological_process (BP), molecular_function(MF), cellular_component(CC)), 
                     or None if no GO terms are found or if there's an error.
    """

    # Step 1: Query MyGene.info to get GO term IDs
    mygene_url = f"https://mygene.info/v3/query?q={gene_name}&species=human&fields=go"
    response = requests.get(mygene_url)

    if response.status_code != 200:
        return None

    data = response.json()

    if "hits" not in data or not data["hits"]:
        return None

    # Extract GO term IDs
    go_terms = {
        "biological_process": [],
        "molecular_function": [],
        "cellular_component": []
    }
    go_ids = set()

    # Loop through hits to find GO terms
    for hit in data["hits"]:
        if "go" in hit:
            go_data = hit["go"]
            for category in ["BP", "MF", "CC"]:  # Ensure we check all three
                if category in go_data:
                    terms = go_data[category]
                    if isinstance(terms, dict):  # If only one GO term category is present
                        terms = [terms]
                    for term in terms:
                        go_id = term.get("id")
                        go_name = term.get("term")
                        if go_id:
                            go_ids.add(go_id)
                            go_entry = {
                                "id": go_id,
                                "name": go_name or "Unknown",
                                "category": category
                            }
                            if category == "BP":
                                go_terms["biological_process"].append(go_entry)
                            elif category == "MF":
                                go_terms["molecular_function"].append(go_entry)
                            elif category == "CC":
                                go_terms["cellular_component"].append(go_entry)

    if not go_ids:
        return None

    # Step 2: Query EBI QuickGO API with retrieved GO term IDs
    quickgo_url = "https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms"
    go_id_list = ",".join(go_ids)
    headers = {"Accept": "application/json"}

    try:
        quickgo_response = requests.get(f"{quickgo_url}/{go_id_list}", headers=headers, timeout=10)

        if quickgo_response.status_code == 200:
            quickgo_data = quickgo_response.json().get("results", [])

            for result in quickgo_data:
                go_id = result.get("id")
                go_description = result.get("definition", {}).get("text", "No description available")
                go_category = result.get("aspect", "Unknown")
                go_synonyms = [syn["name"] for syn in result.get("synonyms", []) if "name" in syn]

                # Match the GO ID with its category
                for category, terms in go_terms.items():
                    for term in terms:
                        if term["id"] == go_id:
                            term["description"] = go_description
                            term["synonyms"] = go_synonyms

            return go_terms

        else:
            return None

    except requests.exceptions.RequestException as e:
        return None



def load_tajima_d_results(directory):
    """
    Loads Tajima's D statistics from files in a directory into the database.

    The function expects the files to follow a specific naming convention and format:

    File Naming Convention:
      - '<PopulationName>_chr<ChromosomeNumber>.Tajima.D'
      - Example: 'BEB_chr10.Tajima.D'

    File Content Format:
      - Tab-separated values (.tsv)
      - Columns:
          - 1st column (index 0): Window Midpoint (ignored)
          - 2nd column (index 1): Start position of the genomic bin
          - 3rd column (index 2): Number of SNPs in the bin (nSNPS)
          - 4th column (index 3): Tajima's D value (TajimaD)

    Example File Content ('BEB_chr10.Tajima.D'):
    '''
    Position    Start   nSNPS   TajimaD
    15000       10000   25      -1.237
    25000       20000   30      0.983
    35000       30000   40      1.456
    45000       40000   50      -0.678
    '''

    The function performs the following steps:
    1. Deletes all existing Tajima's D records in the 'tajima_d_results' table Need to change this in the future.
    2. Iterates through each file in the specified directory.
    3. Extracts population and chromosome information from the filename.
    4. Reads the file line by line, skipping the header.
    5. Parses each line, extracting the relevant data.
    6. Filters out lines with missing or invalid Tajima's D values.
    7. If a population, chromosome, and bin combination does not already exist, it inserts and update the database.
    8. Inserts the valid data into the 'tajima_d_results' table.

    Assumptions:
    - The bin size is expected to be 10kb ('bin_end = bin_start + 10000').
    - The first line in each file is a header and is skipped.
    """
    with db.session.begin():
        # Iterate through all files in the directory
        for filename in os.listdir(directory):
            if filename.endswith(".Tajima.D"): # Process only '.Tajima.D' files
                # Extract population and chromosome number from filename
                population, chromosome = filename.split("_")[0], filename.split("_")[1].split(".")[0].replace("chr", "")

                with open(os.path.join(directory, filename), "r") as file:
                    next(file)# Skip the header line
                    for line in file:
                        columns = line.strip().split("\t") # Split columns based on tab separation
                        
                        # Skip lines with invalid or missing Tajima's D values
                        if len(columns) < 4 or columns[3] in ["", "NA", ".", "nan"]:
                            continue
                            
                        bin_start = int(columns[1])
                        bin_end = bin_start + 10000  # Assuming 10kb bin size
                        n_snps = int(columns[2])
                        tajima_d = float(columns[3])

                        # Check if this population + chromosome + bin already exists
                        existing_tajima = db.session.query(TajimaD).filter_by(
                            population=population, chromosome=chromosome, bin_start=bin_start
                        ).first()
                        
                        if not existing_tajima:
                            # If this combination does NOT exist, insert it
                            db.session.add(TajimaD(
                                population=population,
                                chromosome=chromosome,
                                bin_start=bin_start,
                                bin_end=bin_end,
                                n_snps=n_snps,
                                tajima_d=tajima_d
                            ))
                            
    db.session.commit() # Save changes to the database


def load_clr_results(directory):
    """
    Loads CLR results from a single file into the database.

    The function expects the file to follow a specific naming convention and format:

    File Naming Convention:
      - 'SweeD_Report.<Population>_chr<Chromosome>.txt'
      - Example: 'SweeD_Report.BEB_chr10.txt'

    File Content Format:
      - Tab-separated values (.tsv)
      - Columns:
          - 1st column (index 0): Position
          - 2nd column (index 1): Likelihood/ CLR value
          - 3rd column (index 2): Alpha value

    Example File Content ('SweeD_Report.BEB_chr10.txt'):
    '''
    Position   Likelihood       Alpha
    11501	  2.241451e-01	 3.902013e+01
    12838	  2.254904e-01	 2.918483e-04
    14176	  3.827930e-01	 1.359461e-04
    '''

    The function performs the following steps:
    1. Extracts population and chromosome information from the filename.
    2. Reads the file into a DataFrame, skipping the first two rows.
    3. Cleans up data by removing rows where 'Position' is not a number.
    4. Renames columns to match the table schema.
    5. Adds population and chromosome columns.
    6. Adds a column for id, assuming we want to manually assign ids starting from 1 or max(id)+1.
    7. Reorders columns to match the table schema.
    8. Inserts the valid data into the 'clr_results' table.

    Assumptions:
    - The first two lines in the TSV file are headers and are skipped.
    """
    with db.session.begin():
        max_id = db.session.query(func.max(CLRTest.id)).scalar()
        if max_id is None:
            max_id = 0  # Start ID numbering from 1 if table is empty

        for filename in os.listdir(directory):
            if filename.startswith('SweeD_Report'):
                # Extract population and chromosome from filename
                parts = filename.split('.')
                population = parts[1].split('_')[0]
                chromosome = parts[1].split('_')[1].replace('chr', '')

                file_path = os.path.join(directory, filename)

                # Read the TSV file, skipping the first two rows
                df = pd.read_csv(file_path, sep='\t', skiprows=2)

                # Clean up data by removing rows where 'Position' is not a number
                df = df[pd.to_numeric(df['Position'], errors='coerce').notnull()]

                # Rename columns
                df.columns = ['position', 'clr', 'alpha']

                # Add population and chromosome columns
                df['population'] = population
                df['chromosome'] = chromosome

                # Assign unique IDs starting from max_id + 1
                df['id'] = range(max_id + 1, max_id + 1 + len(df))
                max_id += len(df)  # Update max_id for the next batch

                # Reorder columns
                df = df[['id', 'population', 'chromosome', 'position', 'clr', 'alpha']]

                # Insert only if the entry does not already exist
                for _, row in df.iterrows():
                    existing_clr = db.session.query(CLRTest).filter_by(
                        population=row['population'],
                        chromosome=row['chromosome'],
                        position=row['position']
                    ).first()

                    if not existing_clr:
                        clr_test = CLRTest(
                            id=row['id'],
                            population=row['population'],
                            chromosome=row['chromosome'],
                            position=row['position'],
                            clr=row['clr'],
                            alpha=row['alpha']
                        )
                        db.session.add(clr_test)
    try: 
        db.session.commit()  # Save changes to the database
    except Exception as e:
        db.session.rollback()
        print(f"Error: {e}")
    finally:
        db.session.close()



def get_gene_coordinates_ensembl(gene_name):
    """
    Fetches gene information from Ensembl, including Ensembl ID, chromosome, start/end position, biotype, and description.

    Args:
        gene_name (str): The name of the gene to query.

    Returns:
        dict: A dictionary containing the Ensembl gene ID, biotype, chromosome, start/end position, description, and assembly name.
              Returns None if the gene is not found.
    """
        
    # Step 1: Query Ensembl to retrieve the Ensembl Gene ID
    url = f"https://rest.ensembl.org/xrefs/symbol/homo_sapiens/{gene_name}?content-type=application/json"
    response = requests.get(url)

    if response.status_code == 200:
        data = response.json()
        if data:
            ensembl_gene_id = data[0]["id"]  # Extract Ensembl Gene ID

            # Step 2: Fetch gene details using Ensembl Gene ID
            gene_url = f"https://rest.ensembl.org/lookup/id/{ensembl_gene_id}?content-type=application/json"
            gene_response = requests.get(gene_url)

            if gene_response.status_code == 200:
                gene_data = gene_response.json()
                
                # Return only the required fields
                return {
                    "ensembl_id": ensembl_gene_id,
                    "name": gene_data.get("display_name", "N/A"),
                    "biotype": gene_data.get("biotype", "N/A"),
                    "chromosome": gene_data.get("seq_region_name", "N/A"),
                    "start": gene_data.get("start", "N/A"),
                    "end": gene_data.get("end", "N/A"),
                    "assembly_name": gene_data.get("assembly_name", "N/A"),
                    "description": gene_data.get("description", "N/A")
                }
    
    return None  # Return None if the gene is not found


def load_fst_snp_results(directory):
    """
    Loads and adds SNP-based FST data from CSV files into the database.
    If the CSV does not have 'Ref' and 'Alt' columns, the values will be set to None.

    """
    import os
    from models import FstSNP

    processed_files = 0
    new_records = 0
    updated_records = 0

    try:
        # Loop over every file in the specified directory
        for filename in os.listdir(directory):
            if filename.endswith(".csv"):
                filepath = os.path.join(directory, filename)
                df = pd.read_csv(filepath)

                # Start a transaction for the current file
                with db.session.begin():
                    # Iterate over each row in the DataFrame
                    for _, row in df.iterrows():
                        # Use row.get() to safely retrieve 'Ref' and 'Alt' values.
                        # If the column is missing, row.get() returns an empty string (default value).
                        ref_value = row.get('Ref', "")
                        alt_value = row.get('Alt', "")
                        
                        ref_allele = str(ref_value).strip() or None
                        alt_allele = str(alt_value).strip() or None

                        existing = FstSNP.query.filter_by(snp_id=row['SNP']).first()
                        
                        if existing:
                            # Update existing record with null-safe values
                            existing.gene = row['Gene']
                            existing.chromosome = str(row['Chromosome'])
                            existing.position = int(row['Position'])
                            existing.reference_allele = ref_allele
                            existing.alternative_allele = alt_allele
                            existing.fst_beb = float(row['FST_BEB'])
                            existing.fst_gih = float(row['FST_GIH'])
                            existing.fst_itu = float(row['FST_ITU'])
                            existing.fst_pjl = float(row['FST_PJL'])
                            existing.fst_stu = float(row['FST_STU'])
                            updated_records += 1
                        else:
                            # Insert new SNP data into the database    
                            db.session.add(FstSNP(
                                snp_id=row['SNP'],
                                gene=row['Gene'],
                                chromosome=str(row['Chromosome']),
                                position=int(row['Position']),
                                reference_allele=ref_allele,
                                alternative_allele=alt_allele,
                                fst_beb=float(row['FST_BEB']),
                                fst_gih=float(row['FST_GIH']),
                                fst_itu=float(row['FST_ITU']),
                                fst_pjl=float(row['FST_PJL']),
                                fst_stu=float(row['FST_STU'])
                            ))
                            new_records += 1
                    processed_files += 1

        db.session.commit()
        print(f"""
        FST Data Load Complete!
        Processed files: {processed_files}
        New SNPs added: {new_records}
        Existing SNPs updated: {updated_records}
        Null Ref alleles handled: {df['Ref'].isna().sum() if 'Ref' in df.columns else 'Column Missing'}
        Null Alt alleles handled: {df['Alt'].isna().sum() if 'Alt' in df.columns else 'Column Missing'}
        """)

    except Exception as e:
        db.session.rollback()
        print(f"Error loading FST data: {str(e)}")
        raise






def get_fst_data(chromosome, region=None, populations=None):
    """
    Fetch FST statistics for a given chromosome (or region) and selected populations.
    
    Args:
        chromosome (str): Chromosome to filter by.
        region (tuple, optional): (start, end) positions defining the genomic region.
        populations (list, optional): List of populations (e.g., ["BEB", "GIH"]) to filter by.
        
    Returns:
        tuple: A dictionary of FST values (keyed by population) and a dictionary of summary statistics 
               (mean and standard deviation) for each population.
    """
    query = FstSNP.query.filter(FstSNP.chromosome == chromosome)
    if region:
        start, end = region
        query = query.filter(FstSNP.position >= start, FstSNP.position <= end)
    
    # Initialize dictionary for each selected population
    fst_data = {pop: [] for pop in populations} if populations else {}
    
    for fst in query.all():
        for pop in populations:
            # Map population code to corresponding FST column (e.g., "BEB" -> fst_beb)
            col = f'fst_{pop.lower()}'
            value = getattr(fst, col, None)
            fst_data[pop].append({
                "snp_id": fst.snp_id,
                "position": fst.position,
                "fst": value
            })
    
    # Calculate summary statistics for each population
    summary_stats = {}
    for pop, data in fst_data.items():
        values = [d["fst"] for d in data if d["fst"] is not None]
        mean_val = np.mean(values) if values else None
        std_val = np.std(values) if values else None
        summary_stats[pop] = {
            "mean_fst": round(mean_val, 4) if mean_val is not None else "N/A",
            "std_dev_fst": round(std_val, 4) if std_val is not None else "N/A"
        }
    
    return fst_data, summary_stats


def download_fst_data():
    """
    Generates a text file containing FST statistics for a specified genomic region.
    
    Request Parameters:
        - chromosome (str)
        - start (int, optional)
        - end (int, optional)
        - selected_population (list): Populations to filter by.
        - gene_name (str, optional): If provided, the gene's coordinates are used.
        
    Returns:
        flask.Response: A text file with FST values for each SNP in the region and summary statistics.
    """
    user_selected_chromosome = request.args.get("chromosome")
    gene_name = request.args.get("gene_name")
    start = request.args.get("start", type=int)
    end = request.args.get("end", type=int)
    populations = request.args.getlist("selected_population")
    
    # If gene name is provided and region is not defined, fetch coordinates
    if gene_name and (start is None or end is None):
        gene_info = get_gene_coordinates_ensembl(gene_name)
        if gene_info:
            # Validate that the gene is on the user-selected chromosome
            if gene_info["chromosome"] != user_selected_chromosome:
                return jsonify({
                    "error": f"Gene '{gene_name}' is not located on the selected chromosome {user_selected_chromosome}"
                }), 400
            # Use gene coordinates
            chromosome = gene_info["chromosome"]
            start = gene_info["start"]
            end = gene_info["end"]
        else:
            return jsonify({"error": f"Gene '{gene_name}' not found in Ensembl"}), 400
    else:
        # No gene provided; use the user-selected chromosome
        chromosome = user_selected_chromosome
    
    if start is None or end is None:
        return jsonify({"error": "Invalid region parameters. Must provide a valid start and end position or a gene name."}), 400
    
    try:
        fst_data, fst_summary_stats = get_fst_data(chromosome, (start, end), populations)
    except Exception as e:
        return jsonify({"error": f"Error retrieving FST data: {str(e)}"}), 500
    
    # Handle the case where no data exists
    if not fst_data or all(len(entries) == 0 for entries in fst_data.values()):
        response_text = f"No FST data found for Chromosome {chromosome}, Region {start}-{end}.\n"
        response = Response(response_text, mimetype="text/plain")
        response.headers["Content-Disposition"] = f"attachment; filename=FST_chr{chromosome}_{start}_{end}.txt"
        return response

    # Generate file content
    file_content = ["Population\tSNP ID\tPosition\tFST"]
    for pop, entries in fst_data.items():
        for entry in entries:
            file_content.append(f"{pop}\t{entry['snp_id']}\t{entry['position']}\t{entry['fst']:.4f}")
    
    # Append summary statistics
    file_content.append("\nSummary Statistics")
    file_content.append("Population\tMean FST\tStd Dev FST")
    for pop, stats in fst_summary_stats.items():
        file_content.append(f"{pop}\t{stats['mean_fst']}\t{stats['std_dev_fst']}")
    
    response = Response("\n".join(file_content), mimetype="text/plain")
    response.headers["Content-Disposition"] = f"attachment; filename=FST_chr{chromosome}_{start}_{end}.txt"
    return response















###### function to get mapped genes #####
# import requests
# import pandas as pd
# import time
# import os

# # Load your SNP data
# file_path = os.path.expanduser("~/Documents/BIO727P-GROUP_PROJECT/backup_web/converted_positions.csv")
# df = pd.read_csv(file_path)

# # Extract unique rsIDs and ensure it's a list (not NumPy array)
# rsids = df["dbSNP"].dropna().astype(str).tolist()

# # Ensembl VEP API URL
# ENSEMBL_VEP_URL = "https://rest.ensembl.org/vep/human/id"

# # Function to query Ensembl VEP API for mapped genes
# def get_ensembl_mapped_genes(rsid_batch, retries=3):
#     url = f"{ENSEMBL_VEP_URL}"
#     headers = {"Content-Type": "application/json", "Accept": "application/json"}
    
#     # Ensure rsid_batch is a list (fixing TypeError issue)
#     rsid_batch = list(rsid_batch)  

#     for attempt in range(retries):
#         try:
#             response = requests.post(url, json={"ids": rsid_batch}, headers=headers, timeout=10)
            
#             if response.status_code == 200:
#                 data = response.json()
#                 mapped_genes = {}

#                 for entry in data:
#                     rsid = entry.get("id")
#                     gene_names = set()
                    
#                     # Extract genes from transcript consequences
#                     for transcript in entry.get("transcript_consequences", []):
#                         gene_name = transcript.get("gene_symbol")
#                         if gene_name:
#                             gene_names.add(gene_name)

#                     mapped_genes[rsid] = ", ".join(gene_names) if gene_names else "No Mapped Gene"

#                 return mapped_genes

#             elif response.status_code == 429:  # Too many requests
#                 print(f"⚠️ API Rate Limit reached. Retrying in 10 seconds...")
#                 time.sleep(10)

#         except requests.RequestException as e:
#             print(f"Error fetching {rsid_batch}: {e}")

#     return {rsid: "No Mapped Gene" for rsid in rsid_batch}  # Default if API fails

# # Dictionary to store results
# mapped_genes_dict = {}

# # Process in batches of 50 to reduce API requests
# batch_size = 50
# for i in range(0, len(rsids), batch_size):
#     batch = rsids[i : i + batch_size]
#     mapped_genes_dict.update(get_ensembl_mapped_genes(batch))

#     # Print progress
#     print(f"⏳ Processed {i + len(batch)} SNPs... Sleeping for 1 second")
#     time.sleep(1)  # Pause to prevent hitting API limits

# # Convert results to DataFrame
# mapped_genes_df = pd.DataFrame(mapped_genes_dict.items(), columns=["rsID", "Mapped_Genes"])

# # Merge with original data
# df = df.merge(mapped_genes_df, left_on="dbSNP", right_on="rsID", how="left")

# # Save updated file
# output_path = "updated_snp_data_with_mapped_genes.csv"
# df.to_csv(output_path, index=False)

# print(f"Updated file saved as {output_path}")
# print(mapped_genes_df)


# mapped_genes_df




##### Function to convert GRCh37 postion to GRCh38 T2dkp postions from the T2DKP CSV file #######

# # Function to convert a position from GRCh37 to GRCh38 (returns only start position)
# def convert_grch37_to_grch38(chromosome, position):
#     url = f"https://rest.ensembl.org/map/human/GRCh37/{chromosome}:{position}..{position}/GRCh38?"
#     headers = {"Content-Type": "application/json"}
    
#     response = requests.get(url, headers=headers)
    
#     if response.status_code == 200:
#         data = response.json()
#         if "mappings" in data and len(data["mappings"]) > 0:
#             return data["mappings"][0]["mapped"]["start"]  # Only returning start position
    
#     return None  # Return None if mapping fails

# # Load your dataset (assuming it's in a CSV file)
# df = pd.read_csv("significant_snps.txt")

# # Convert each position
# converted_start_positions = []

# for index, row in df.iterrows():
#     chrom = row["chromosome"]
#     pos = row["position"]
    
#     new_start = convert_grch37_to_grch38(chrom, pos)
#     converted_start_positions.append(new_start)
    
#     # To avoid API rate limits
#     time.sleep(0.2)

# # Add new positions to the dataframe
# df["GRCh38_start"] = converted_start_positions

# # Save the updated dataset
# df.to_csv("converted_positions.csv", index=False)

# print("Conversion complete. Results saved to converted_positions.csv")


import pandas as pd
# import requests
# import time

# # Function to convert GRCh37 position to GRCh38 using Ensembl API for FST_results.csv
# # Note: This function uses the Ensembl REST API to convert positions from GRCh37 to GRCh38 make sure all column names are matching
# # The function returns the GRCh38 start position if successful, otherwise None
# def convert_grch37_to_grch38(chromosome, position):
#     url = f"https://rest.ensembl.org/map/human/GRCh37/{chromosome}:{position}..{position}/GRCh38?"
#     headers = {"Content-Type": "application/json"}
    
#     try:
#         response = requests.get(url, headers=headers)
#         if response.status_code == 200:
#             data = response.json()
#             if data.get("mappings"):
#                 return data["mappings"][0]["mapped"]["start"]  # Return GRCh38 start position
#         return None  # Return None if mapping fails
#     except Exception as e:
#         print(f"Error converting chr{chromosome}:{position} - {str(e)}")
#         return None

# # Load the CSV file
# df = pd.read_csv("FST_results.csv")

# # Convert positions and update the DataFrame
# grch38_positions = []
# total_snps = len(df)

# for index, row in df.iterrows():
#     chrom = row["Chromosome"]  # Ensure column name matches (case-sensitive)
#     pos = row["Position"]
    
#     new_pos = convert_grch37_to_grch38(chrom, pos)
#     grch38_positions.append(new_pos if new_pos is not None else pos)  # Fallback to original if conversion fails
    
#     # Progress tracking
#     print(f"Processed SNP {index + 1}/{total_snps}: {row['SNP']} -> {new_pos}")
    
#     time.sleep(0.2)  # Avoid API rate limits

# # Replace the Position column with GRCh38 coordinates
# df["Position"] = grch38_positions

# # Save to new CSV
# df.to_csv("FST_results_GRCh38.csv", index=False)
# print("Conversion complete. Updated data saved to FST_results_GRCh38.csv")
