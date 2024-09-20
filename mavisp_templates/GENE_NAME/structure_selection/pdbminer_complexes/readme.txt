# AFTER HAVING SUCCESSFULLY RUN PDBMINER YOU CAN FILTER ON COMPLEXES:

#uniprot ID starting and ending residue of the domain
bash do.sh P36896 29 110

#activate the python environment (module load python)

#to run for mavisp - replace uniprot ID with the one of your protein
python find_PDBminer_complexes.py -i P35869_all.csv --binding_interface -d 10 -s 35 -e 591 -o P35869_filtered.csv


#general info
#SIMPLE USE, just specify input and output
python find_PDBminer_complexes.py -i uniprot_all.csv -o uniprot_filtered.csv

#DOMAIN USE, define your domain of interest
python find_PDBminer_complexes.py -i uniprot_all.csv -o uniprot_filtered.csv -s start_residue -e end_residue

#INTERFACE USE:
python find_PDBminer_complexes.py -i uniprot_all.csv -o uniprot_filtered.csv --binding_interface

#INTERFACE USE WITH DEFINED DISTANCE:
python find_PDBminer_complexes.py -i uniprot_all.csv -o uniprot_filtered.csv --binding_interface -d 5 -s start_residue -e end_residue


