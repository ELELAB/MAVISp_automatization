#Requirements
- to have generated the data for the mutations of interest with cancermuts in raw_data
- the pdb file(s) of the model(s) you'll use in the folder so that the check of the WT residue can occur properly

#symbolic link to cancermuts metatable
#symbolic link to pdbs

# The script get_mutlists.py takes in input the cancermut's metatable
# and return one or more mutlist.txt files depending on the commandline flags
# containing only the mutations with a source and occurring at positions covered by the selected structure

# If you have the full-length protein or a trimmed model in which only the N-ter and/or C-ter portions have been removed
# you can follow the single domain example (i.e., it is enough to specify the range of residue after the -d flag [-d 1:756]).
# If you have trimmed the structure in multiple domains, follow the multiple domains example:
# For instance, if your structures are covering two ranges of residues (e.g., 1:341 and 501:756),
# you have to specify both the ranges (space-separated) after the -d flag (i.e., -d 1:341 501:756).

# Compatibility of the mutlist:
# one letter code (R13L) - automatically generated
# mutatex (RA13L) - generated if the '-M' flag is used (required -ch flag [default = A)]
# rosetta (A.R.13.L A) - generated if the '-R' flag is used (required -ch flag with the chain [default = A])
# hgvs (p.Arg13Leu) - generated if the '-H' flag is used 
# cabsflex (4P7A.pdb A 13 ARG LEU) - generated if the '-C' flag is used (required the -ch [default = A])
# ELM mutlist - automatically generated

# for further information see GitHub

. /usr/local/envs/cancermuts/bin/activate

#For info:
python get_mutlists.py -h

#single domain
python get_mutlists.py -m metatable_pancancer_ADCK1.csv -d 1:530 -M -R -H -C -p Q86TW2.pdb  

# multiple domains 
python get_mutlists.py -m metatable_pancancer_MLH1.csv -d 1:341 501:756 -M -R -H -C -p 4P7A.pdb 3RBN.pdb
