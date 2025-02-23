import pandas as pd
from models import db, SNP, TajimaD, Fst
import requests
import plotly.graph_objects as go
from urllib.parse import quote_plus
import os

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





# ENSG00000147883

# import requests
# headers = {'content-type': 'application/x-www-form-urlencoded'}
# params = 'ids=ENSG00000147883,695&fields=go,symbol,refseq.rna'
# res = requests.post('http://mygene.info/v3/gene', data=params, headers=headers)

# res.json()

# def load_population_frequencies(base_directory):
#     import pandas as pd
#     from main import app
#     import os

#     populations = ['BEB', 'GIH', 'ITU', 'PJL', 'STU']

#     with app.app_context():
#         for population in populations:
#             population_dir = os.path.join(base_directory, population, f"{population}_t2dfreqs")
#             if os.path.exists(population_dir):
#                 for file_name in os.listdir(population_dir):
#                     if file_name.endswith(".txt"):
#                         file_path = os.path.join(population_dir, file_name)
#                         df = pd.read_csv(file_path, sep='\t', header=None, names=['chromosome', 'position', 'snp_id', 'reference_allele', 'alternative_allele', 'frequency'])

#                         for _, row in df.iterrows():
#                             freq_entry = PopulationFrequency(
#                                 population=population,
#                                 chromosome=row['chromosome'],
#                                 position=row['position'],
#                                 snp_id=row['snp_id'],
#                                 reference_allele=row['reference_allele'],
#                                 alternative_allele=row['alternative_allele'],
#                                 frequency=row['AF']
#                             )
#                             db.session.add(freq_entry)
#         db.session.commit()

# print("Population frequencies database model updated with foreign key relationship and MAF calculation.")
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



# import requests
# import pandas as pd
# import time

# # Load your SNP data
# file_path = "~/Documents/BIO727P-GROUP_PROJECT/backup_web/converted_positions.csv"  # Update with actual file path
# df = pd.read_csv(file_path)

# # Extract unique rsIDs
# rsids = df["dbSNP"].dropna().unique()

# # Function to query NCBI API for mapped genes
# def get_ncbi_mapped_gene(rsid):
#     url = f"https://api.ncbi.nlm.nih.gov/variation/v0/beta/refsnp/{rsid.lstrip('rs')}"
#     response = requests.get(url)

#     if response.status_code == 200:
#         data = response.json()
#         genes = set()
#         try:
#             for annotation in data.get("primary_snapshot_data", {}).get("assembly_annotation", []):
#                 for gene in annotation.get("genes", []):
#                     if "locus" in gene:
#                         genes.add(gene["locus"])
#         except KeyError:
#             return None  # No mapped genes found

#         return list(genes) if genes else None  # Return unique gene names
#     return None  # If API request fails

# # Function to query ClinicalTables API for mapped genes
# def get_clinicaltables_mapped_gene(rsid):
#     url = f"https://clinicaltables.nlm.nih.gov/api/snps/v3/search"
#     params = {"terms": rsid}
#     response = requests.get(url, params=params)

#     if response.status_code == 200:
#         data = response.json()
#         try:
#             for snp_info in data[3]:  # SNP data is in the fourth list
#                 if snp_info[0] == rsid:  # Ensure correct SNP
#                     return [snp_info[4]] if snp_info[4] else None  # Return gene if found
#         except (IndexError, TypeError):
#             return None  # No mapped gene found
#     return None  # If API request fails

# # Query both APIs in batches to avoid rate limits
# mapped_genes_dict = {}

# for i, rsid in enumerate(rsids):
#     mapped_gene = get_ncbi_mapped_gene(rsid)

#     # If NCBI fails, use ClinicalTables API
#     if not mapped_gene:
#         mapped_gene = get_clinicaltables_mapped_gene(rsid)

#     mapped_genes_dict[rsid] = mapped_gene

#     if i % 10 == 0:  # Prevent hitting API rate limits
#         time.sleep(1)

# # Convert results to DataFrame
# mapped_genes_df = pd.DataFrame(list(mapped_genes_dict.items()), columns=["rsID", "Mapped_Genes"])

# # Merge with original data
# df = df.merge(mapped_genes_df, left_on="dbSNP", right_on="rsID", how="left")

# # Save updated file
# output_path = "updated_snp_data_with_mapped_genes.csv"
# df.to_csv(output_path, index=False)

# print(f"Updated file saved as {output_path}")
# print(mapped_genes_df)
