## To run Alphafill co-factor you need to start at the Alphafill website (see guideline doc).

# 1 search your protein on uniprot accession number
# 2 Choose 30% identity (indicating < 30%)
# 3 Identify Compounds and Asym letters. 
# 3a - choose the natural compounds (metals, ATP and so forth), not drugs for example. 
#      you can always refer to the PDB-ID to investigate what
#      the compound is.
# 4 run the script
# 5 identify the residues coordinating the compund in the output inter_res.csv

#activate the environment:
source /usr/local/envs/py310/bin/activate

#run the script, {} to show where you input, should not be in the command
#you can add as many ASYM letter and compound options as you like.
python alphafill_cofactor.py -u {uniprot_accession_number} -cha {chain} -co {ASYM letter} {compund} {ASYM letter} {compound} -t 4

#Example: 
#python alphafill_cofactor.py -u Q86TW2 -cha A -co B ADP C ATP -t 4

