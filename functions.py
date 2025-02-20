import pandas as pd
from models import db, SNP, TajimaD, Fst
import requests
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
    Fetch Gene Ontology terms for a given gene and retrieve detailed GO descriptions.
    """
    # Step 1: Get Ensembl Gene ID for the given gene name
    lookup_url = f"https://rest.ensembl.org/xrefs/symbol/homo_sapiens/{quote_plus(gene_name)}?content-type=application/json"
    response = requests.get(lookup_url)

    if response.status_code == 200:
        data = response.json()
        if data:
            ensembl_gene_id = data[0]["id"]  # Extract Ensembl Gene ID

            # Step 2: Fetch GO terms from MyGeneInfo API
            url = f"https://mygene.info/v3/gene/{ensembl_gene_id}"
            params = {"fields": "go"}
            response = requests.get(url, params=params, headers={"accept": "application/json"})

            if response.status_code == 200:
                go_data = response.json().get("go", {})

                # Extract top 2 GO terms per category
                biological_process = go_data.get("BP", [])[:2]
                molecular_function = go_data.get("MF", [])[:2]
                cellular_component = go_data.get("CC", [])[:2]

                # Fetch detailed GO term info
                enriched_go_terms = {}
                for category, go_terms in zip(["biological_process", "molecular_function", "cellular_component"], 
                                              [biological_process, molecular_function, cellular_component]):

                    enriched_go_terms[category] = []
                    
                    for go_term in go_terms:
                        go_id = go_term["id"]
                        go_detail_url = f"https://api.geneontology.org/api/ontology/term/{quote_plus(go_id)}"
                        go_response = requests.get(go_detail_url)

                        if go_response.status_code == 200:
                            enriched_data = go_response.json()
                            enriched_go_terms[category].append({
                                "id": go_id,
                                "name": enriched_data.get("label", "N/A"),
                                "definition": enriched_data["definition"]["text"] if isinstance(enriched_data.get("definition"), dict) else enriched_data.get("definition", "No definition available"),
                                "synonyms": enriched_data.get("synonyms", []),
                                "namespace": enriched_data.get("namespace", "N/A")
                            })

                return enriched_go_terms  # Return enriched GO terms
    
    return None  # Return None if no data is found



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

        
        