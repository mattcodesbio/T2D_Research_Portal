
# as we dont have a datbase right now the below fucntion simulates retrieving SNP data
def get_snp_info(snp_name, chromosome, start, end, gene_name):
    # Simulating SNP data retrieval
    # In practice, this would involve querying a database or an API
    if not snp_name:
        return None
    
    snp_data = [
        {
            'name': 'rs12345',
            'position': f'Chromosome {chromosome}, {start}-{end}',
            'p_value': 0.03,
            'gene': gene_name if gene_name else 'GeneA',
            "positive_selection_stat_1": "2.5",
            "positive_selection_stat_2": "3.8",
            # 'selection_summary': 'Positive selection observed in South Asian populations (p-value < 0.05).'
        },
        {
            'name': 'rs67890',
            'chromosome': '2',
            'position': '234567-234890',
            'p_value': 0.01,
            'gene': 'GENE2',
            'positive_selection_stat_1': "3.0",
            'positive_selection_stat_2': "4.5",
        },
    ]

    # Look up SNP name in the simulated database
    # Searching for the SNP by name
    for snp in snp_data:
        if snp['name'] == snp_name:
            # If SNP matches, return the SNP details
            return snp

    # If SNP not found, return None
    return None