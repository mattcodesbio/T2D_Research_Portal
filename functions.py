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
    with db.session.begin():  # No need to import app; db is initialized correctly
        db.session.query(SNP).delete()
        for _, row in df.iterrows():
            db.session.add(SNP(
                snp_id=row["dbSNP"],
                chromosome=str(row["chromosome"]),
                grch38_start=int(row["GRCh38_start"]),
                gene_name=row["gene_name"].strip("[]").replace("\"", "") if row["gene_name"] else None,
                p_value=row["pValue"],
                reference_allele=row["reference"],
                alternative_allele=row["alt"],
                consequence=row["consequence"]
            ))
    db.session.commit()



def get_gene_ontology_terms(gene_name):
    """
    Step 1: Fetch GO term IDs from MyGene.info API
    Step 2: Query QuickGO API to fetch details for each GO term ID
    """

    # Step 1: Query MyGene.info to get GO term IDs
    mygene_url = f"https://mygene.info/v3/query?q={gene_name}&species=human&fields=go"

    response = requests.get(mygene_url)

    if response.status_code != 200:
        return None  # API request failed

    data = response.json()

    if "hits" not in data or not data["hits"]:
        return None  # No results found

    # Extract GO term IDs from MyGene.info response
    go_terms = {
        "biological_process": [],
        "molecular_function": [],
        "cellular_component": []
    }

    go_categories = data["hits"][0].get("go", {})
    go_ids = set()  # Store unique GO term IDs

    for category, terms in go_categories.items():
        for term in terms:
            go_id = term["id"]
            go_ids.add(go_id)  # Collect GO IDs for QuickGO lookup

    if not go_ids:
        return None  # No GO terms found

    # Step 2: Query EBI QuickGO API with retrieved GO term IDs
    quickgo_url = "https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms"

    # Convert GO IDs to a comma-separated string
    go_id_list = ",".join(go_ids)

    quickgo_response = requests.get(f"{quickgo_url}/{go_id_list}", headers={"Accept": "application/json"})

    if quickgo_response.status_code == 200:
        quickgo_data = quickgo_response.json().get("results", [])

        for result in quickgo_data:
            go_id = result["id"]
            go_name = result["name"]
            go_category = result["aspect"]  # BP, MF, CC
            go_description = result.get("definition", {}).get("text", "No description available")

            # Extract only the synonym names (ignore type)
            go_synonyms = [syn["name"] for syn in result.get("synonyms", []) if "name" in syn]

            go_entry = {
                "id": go_id,
                "name": go_name,
                "description": go_description,
                "synonyms": go_synonyms  # Now a clean list of names
            }

            if go_category == "biological_process":
                go_terms["biological_process"].append(go_entry)
            elif go_category == "molecular_function":
                go_terms["molecular_function"].append(go_entry)
            elif go_category == "cellular_component":
                go_terms["cellular_component"].append(go_entry)

    return go_terms  # Returns structured GO terms categorized by BP, MF, and CC




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

        
        