source /usr/local/envs/PDBminer/bin/activate
PDBminer -g $1 -u $2 -n 1 -f csv
#$1 gene_name
#$2 uniprot accession number
