import pandas as pd
from models import db, SNP, TajimaD, Fst, CLRTest
import requests
import numpy as np
import os
from models import db, SNP, TajimaD
from flask import jsonify, Response, request

def get_snp_info(snp_id=None, chromosome=None, start=None, end=None, gene_name=None):
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
    
    return results if results else None

def get_tajima_d_data(chromosome, region=None, populations=None):
    """Fetch Tajima's D statistics for a chromosome or a region."""
    query = TajimaD.query.filter(TajimaD.chromosome == chromosome)

    if populations:
        query = query.filter(TajimaD.population.in_(populations))

    if region:
        start, end = region
        query = query.filter(TajimaD.bin_start >= start, TajimaD.bin_end <= end)

    tajima_d_data = {pop: [] for pop in populations} if populations else {}

    for tajima in query.all():
        tajima_d_data[tajima.population].append({
            "bin_start": tajima.bin_start,
            "bin_end": tajima.bin_end,
            "tajima_d": tajima.tajima_d
        })

    summary_stats = {
        pop: {
            "mean": round(np.mean([d["tajima_d"] for d in data]), 4),
            "std_dev": round(np.std([d["tajima_d"] for d in data]), 4)
        }
        for pop, data in tajima_d_data.items()
    }

    return tajima_d_data, summary_stats


def get_t2d_snps(chromosome, start=None, end=None):
    """Fetch significant T2D SNPs within a region."""
    query = SNP.query.filter(SNP.chromosome == chromosome)

    if start and end:
        query = query.filter(SNP.grch38_start >= start, SNP.grch38_start <= end)

    return [{"snp_id": snp.snp_id, "position": snp.grch38_start} for snp in query.all()]


def download_tajima_d_data():
    """Generate a text file containing Tajima's D statistics for a region."""
    # Extract parameters
    chromosome = request.args.get("chromosome")
    gene_name = request.args.get("gene_name")  # Optional gene name
    start = request.args.get("start", type=int)
    end = request.args.get("end", type=int)
    populations = request.args.getlist("selected_population")

    # If gene name is provided, fetch its genomic coordinates
    if gene_name and (start is None or end is None):
        gene_info = get_gene_coordinates_ensembl(gene_name)
        if gene_info:
            chromosome = gene_info["chromosome"]
            start = gene_info["start"]
            end = gene_info["end"]
        else:
            return jsonify({"error": f"Gene '{gene_name}' not found in Ensembl"}), 400

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



def load_snps_from_csv(csv_file):
    df = pd.read_csv(csv_file)

    with db.session.begin():
        db.session.query(SNP).delete()

        for _, row in df.iterrows():
            # Handle Mapped_Genes column properly
            mapped_genes = None
            if pd.notna(row["Mapped_Genes"]):
                mapped_genes = row["Mapped_Genes"].strip()  # Remove surrounding spaces
                if mapped_genes.startswith('"') and mapped_genes.endswith('"'):  
                    mapped_genes = mapped_genes[1:-1]  # Remove surrounding quotes
                mapped_genes = mapped_genes.replace(", ", ",")  # Normalize spacing
                mapped_genes = " ".join(mapped_genes.split(","))  # Ensure proper formatting

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

    db.session.commit()



def get_gene_ontology_terms(gene_name):
    """
    Retrieves GO terms for a given gene name using MyGene.info and QuickGO APIs.
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
                    if isinstance(terms, dict):  # If only one GO term is present
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
    with db.session.begin():
        db.session.query(TajimaD).delete()
        for filename in os.listdir(directory):
            if filename.endswith(".Tajima.D"):
                population, chromosome = filename.split("_")[0], filename.split("_")[1].split(".")[0].replace("chr", "")
                with open(os.path.join(directory, filename), "r") as file:
                    next(file)
                    for line in file:
                        columns = line.strip().split("\t")
                        if len(columns) < 4 or columns[3] in ["", "NA", ".", "nan"]:
                            continue
                        db.session.add(TajimaD(
                            population=population,
                            chromosome=chromosome,
                            bin_start=int(columns[1]),
                            bin_end=int(columns[1]) + 10000,
                            n_snps=int(columns[2]),
                            tajima_d=float(columns[3])
                        ))
    db.session.commit()

def get_gene_coordinates_ensembl(gene_name):
    """
    Fetch chromosome, start, and end position of a gene from Ensembl.
    """
    url = f"https://rest.ensembl.org/xrefs/symbol/homo_sapiens/{gene_name}?content-type=application/json"
    response = requests.get(url)

    if response.status_code == 200:
        data = response.json()
        if data:
            ensembl_gene_id = data[0]["id"]  # Extract Ensembl Gene ID

            # Fetch gene location details
            gene_url = f"https://rest.ensembl.org/lookup/id/{ensembl_gene_id}?content-type=application/json"
            gene_response = requests.get(gene_url)

            if gene_response.status_code == 200:
                gene_data = gene_response.json()
                return {
                    "chromosome": gene_data["seq_region_name"],
                    "start": gene_data["start"],
                    "end": gene_data["end"]
                }
    
    return None  # Return None if gene is not found


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




##### Function to convert GRCh37 postion to GRCh38 T2dkp postions from the T2DKp CSV file #######

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

