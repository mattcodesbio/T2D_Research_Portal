# Type 2 Diabetes (T2D) Research Portal
Welcome to our T2D Research Portal!

This software is intended for researchers, bioinformaticians, and geneticists interested in exploring the genetic basis of Type 2 Diabetes (T2D) and positive selection in South Asian populations. It is free to use under the MIT license.

## Features
- User-friendly interface for genetic data exploration.
- Functional genomics and ontology terms for T2D-associated genes.
- Tools for identifying genetic variants associated with T2D.
- Statistical analysis of positive selection signals.
  
## Installation
- Clone this repository:
```sh
git clone https://github.com/your-repo/T2D-Research-Portal.git
```
- Navigate to the project directory:
 ```sh
cd T2D_Web_development_project
```
Create and activate a virtual environment:
 ```sh
python3 -m venv venv
source venv/bin/activate  # On Windows use `venv\Scripts\activate`
```
- Install the required dependencies:
```sh
pip install -r application_requirements.txt
```
- Run the Flask application:
 ```sh
python main.py
```
## Usage
- Explore T2D-associated SNPs and genes in the provided interface. 
- Perform variant analysis using built-in statistical tools.
- Visualise and compare positive selection of T2D SNPs in South Asian populations.
- Download summary statistics data for further analysis. 
- Explore information on the populations analysed.
  
## Example Usage
- To start a query in the home page, user can use any of the following terms in the search bars:
 ```sh
SNP name (e.g. rs12244851)
Chromosome (e.g. 8)
Genomic position on the chromosome (e.g. Start- 341,000; End- 441,000)
Gene name (e.g. CDKAL1)
```
## License
This project is licensed under the MIT License. See the LICENSE file for details.

## Contact
For any inquiries, please reach out to us via email or GitHub issues.
