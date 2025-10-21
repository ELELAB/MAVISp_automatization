#Aggregate Script

#This script aggregates outputs from the following rules:

## -string2pdb

## -mentha2pdb

## -pdbminer

## -pdbminer_complexes

#It combines the results into a single CSV file containing ranked protein–protein interactors with structural annotations.
#The command executed by the Snakefile is:

./aggregate -m {mentha2pdb_output} -s {string2pdb_output} -p {pdbminer_output} -c {pdbminer_complexes_output}

