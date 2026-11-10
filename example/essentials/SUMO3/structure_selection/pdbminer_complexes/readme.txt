# DESCRIPTION:
# This step:
# 1.Filters complexes extracted by PDBminer based on a 10 Å interaction distance, storing the results in {uniprot_ac}_filtered.csv.
# 2.Downloads the PDB files of the complexes, if available, and stores them in the dir {uniprot_ac}_pdb_complexes/.

# REQUIREMENT:
# PDBminer must have been run successfully

# INPUT FILE:
# {uniprot_ac}_all.csv 
# Must be present in the pdbminer/ dir. 


# COMMANDS OPERATED IN SNAKEFILE:
# Run find_PDBminer_complexes.py to retain filtered complexes and to download PDBs
python find_PDBminer_complexes.py\
               -i ../pdbminer/results/{wildcards.uniprot_ac}/{wildcards.uniprot_ac}_all.csv\
               -o {wildcards.uniprot_ac}_filtered.csv\
               --binding_interface -d 10


# If PDB structures are available, create directory for the PDBs files and move them there:
if ls *.pdb 1> /dev/null 2>&1; then
    mkdir -p {wildcards.uniprot_ac}_pdb_complexes/
    mv *.pdb {wildcards.uniprot_ac}_pdb_complexes/
fi
