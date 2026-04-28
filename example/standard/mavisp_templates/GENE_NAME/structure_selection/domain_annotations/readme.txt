#Step 1: Get Pfam domains from InterPro
*Last updated*: 05/03/2024
#Requirements
have mutlist in cancermuts folder in the project folder of the target protein
#run do.sh as:
bash do.sh UNIPROT_ID mutlist_file.txt

e.g., 
bash do.sh P04080 ../../cancermuts/mutlist_17032023.txt 

#Step 2: Get TED domain classifcations
This Python script queries TED (The Encyclopedia of Domains) domain annotations for a given UniProt ID using the TED API, retaining the TED domain IDs, the respective residue range and the CATH-label. 

## Features
- Fetches TED domain annotations for a given UniProt ID
- Outputs results as a CSV file

## Usage-alone
python ted_api.py -id P36896

## In this folder it is already included in do.sh
