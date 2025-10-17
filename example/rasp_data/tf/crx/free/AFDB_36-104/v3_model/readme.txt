#The environment is activated in the script.

#Add your pdb file here.

#MAKE SURE YOU HAVE:
#	1) CHAIN (if it is not chain A, fix run_analysis -c to be your chain. 
#	2) END/TER/ENDMOL lines have been removed. 
#	3) Amino acids have natural names.

tsp -N 4 bash run_analysis.sh
