import requests
import pandas as pd
import time
from tqdm import tqdm
import logging

# Set up logging for error and status messages
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Configuration
INPUT_CSV = "DATA/SEC_DATA.csv"
OUTPUT_CSV = "TEST2_FST_results.csv"

# Define subpopulations
SAS_SUBPOPS = [
    "1000GENOMES:phase_3:BEB",  # Bengali
    "1000GENOMES:phase_3:GIH",  # Gujarati Indian
    "1000GENOMES:phase_3:ITU",  # Indian Telugu
    "1000GENOMES:phase_3:PJL",  # Punjabi
    "1000GENOMES:phase_3:STU"   # Sri Lankan Tamil
]

EUR_SUBPOPS = [
    "1000GENOMES:phase_3:CEU",  # Utah residents with Northern and Western European ancestry
    "1000GENOMES:phase_3:TSI",  # Toscani in Italy
    "1000GENOMES:phase_3:FIN",  # Finnish
    "1000GENOMES:phase_3:GBR",  # British
    "1000GENOMES:phase_3:IBS"   # Iberian populations in Spain
]

def get_ensembl_data(snp_id):
    """
    Retrieve variant data from Ensembl using the endpoint that returns a 'populations' key.
    """
    url = f"https://rest.ensembl.org/variation/human/{snp_id}?pops=1"
    headers = {"Content-Type": "application/json"}
    try:
        response = requests.get(url, headers=headers)
        response.raise_for_status()
        return response.json()
    except requests.exceptions.RequestException as e:
        logger.error(f"Error fetching data for {snp_id}: {e}")
        return None

def build_population_data(populations, csv_ref, csv_alt):
    """
    Build a dictionary keyed by population code that holds allele-specific data.
    Each value is a dictionary with alleles as keys and a dict containing frequency and allele_count.
    """
    pop_data = {}
    for entry in populations:
        pop_code = entry.get("population")
        allele = entry.get("allele", "").upper()
        try:
            freq = float(entry.get("frequency"))
        except (TypeError, ValueError):
            continue
        allele_count = entry.get("allele_count")
        if allele_count is None:
            continue
        allele_count = int(allele_count)
        if pop_code not in pop_data:
            pop_data[pop_code] = {}
        pop_data[pop_code][allele] = {"frequency": freq, "allele_count": allele_count}
    return pop_data

def get_ref_allele_data(pop_data, pop_code, csv_ref, csv_alt):
    """
    For a given population, retrieve the frequency and count for the reference allele.
    If the reference allele is missing but the alternate is available, infer the reference data.
    Returns: (p, count, total_count) or (None, None, None) if data is insufficient.
    """
    if pop_code not in pop_data:
        return None, None, None
    data = pop_data[pop_code]
    # Compute total allele count in this population (sum over available alleles)
    total = sum(info["allele_count"] for info in data.values())
    if csv_ref in data:
        count = data[csv_ref]["allele_count"]
        p = float(count) / total if total > 0 else None
        return p, count, total
    elif csv_alt in data:
        alt_count = data[csv_alt]["allele_count"]
        p = 1 - (float(alt_count) / total) if total > 0 else None
        count = total - alt_count
        return p, count, total
    else:
        return None, None, None

def calculate_hudson_fst(p1, p2):
    """
    Calculate Hudson's FST for a biallelic SNP given allele frequencies p1 and p2.
    For a single SNP, Hudson's estimator is:
        FST = (p1 - p2)^2 / (p1 + p2 - 2*p1*p2)
    """
    if p1 is None or p2 is None:
        return None
    numerator = (p1 - p2) ** 2
    denominator = p1 + p2 - 2 * p1 * p2
    if denominator == 0:
        return 0.0
    return round(numerator / denominator, 4)

def process_snp(row):
    """Process a single SNP row and calculate Hudson's FST for each SAS subpopulation compared to aggregated European data."""
    snp_id = getattr(row, "dbSNP", None)
    if not snp_id:
        logger.warning("Missing dbSNP id, skipping row")
        return None

    data = get_ensembl_data(snp_id)
    if not data:
        logger.warning(f"No data returned for {snp_id}")
        return None

    # Get allele string from mappings and derive alleles from it
    try:
        mapping = data.get("mappings", [])
        allele_string = mapping[0].get("allele_string", "") if mapping else ""
        if not allele_string:
            logger.warning(f"No allele string for {snp_id}")
            return None
        alleles = [a.upper() for a in allele_string.split("/")]
    except Exception as e:
        logger.error(f"Error processing allele string for {snp_id}: {e}")
        return None

    # Get CSV alleles
    csv_ref = getattr(row, "Reference", "").upper()
    csv_alt = getattr(row, "Alternate", "").upper()

    # For biallelic SNPs, check consistency. For multi-allelic sites, ensure CSV alleles are among the options.
    if len(alleles) != 2:
        if csv_ref in alleles and csv_alt in alleles:
            # Use CSV alleles as provided.
            pass
        else:
            logger.warning(f"CSV alleles {csv_ref}/{csv_alt} not found in allele string for {snp_id}: {allele_string}")
            return None
    else:
        # For biallelic SNPs, if the API order differs, check if it's a simple flip.
        if alleles[0] != csv_ref or alleles[1] != csv_alt:
            if alleles[0] == csv_alt and alleles[1] == csv_ref:
                # Order is flipped, but since we use counts to derive frequencies, no action is needed.
                pass
            else:
                logger.warning(f"Allele mismatch for {snp_id}: API alleles {alleles} vs CSV alleles {csv_ref}/{csv_alt}")
                return None

    # Build a dictionary for population data using both frequency and allele counts.
    pop_data = build_population_data(data.get("populations", []), csv_ref, csv_alt)

    # Aggregate European data by summing counts across the designated European populations.
    total_ref_count = 0
    total_allele_count = 0
    for pop in EUR_SUBPOPS:
        p, count, total = get_ref_allele_data(pop_data, pop, csv_ref, csv_alt)
        if p is not None:
            total_ref_count += count
            total_allele_count += total

    if total_allele_count > 0:
        p_eur = total_ref_count / total_allele_count
    else:
        p_eur = None

    # Prepare the result dictionary with SNP metadata
    results = {
        "SNP": snp_id,
        "Gene": getattr(row, "Gene", None),
        "Chromosome": getattr(row, "Chromosome", None),
        "Position": getattr(row, "Position", None),
        "Ref": csv_ref,
        "Alt": csv_alt
    }

    # For each SAS subpopulation, compute Hudson's FST against the aggregated European frequency.
    for pop in SAS_SUBPOPS:
        pop_short = pop.split(":")[-1]
        p_sas, count_sas, total_sas = get_ref_allele_data(pop_data, pop, csv_ref, csv_alt)
        if p_sas is None or p_eur is None:
            results[f"FST_{pop_short}"] = None
        else:
            results[f"FST_{pop_short}"] = calculate_hudson_fst(p_sas, p_eur)
    return results

def main():
    df = pd.read_csv(INPUT_CSV)
    results = []
    total_snps = len(df)
    logger.info(f"Processing {total_snps} SNPs from your CSV...")

    for row in tqdm(df.itertuples(index=False), total=total_snps):
        result = process_snp(row)
        if result:
            results.append(result)
        time.sleep(0.2)  # Respectful delay between API calls

    result_df = pd.DataFrame(results)
    result_df.to_csv(OUTPUT_CSV, index=False)
    logger.info(f"Done! Results saved to {OUTPUT_CSV}")
    logger.info(f"Processed {len(result_df)}/{total_snps} SNPs successfully")

if __name__ == '__main__':
    main()
