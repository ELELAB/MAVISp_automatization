saturation_mutlist.py script takes as input the uniprot_ac of a protein of interest and returns a file collecting all the possible mutations (included the phosphoryations in the FoldX format: p,y,s) of each residue.
In case a range of residue is provided the script will return only all the possible mutations associated to the residues belonging to that range. The residue range must be in the following format: starting_residue-ending_residue.


Usage:

saturation_mutlist.py -a uniprot_ac -r residue-range -o my_outputfile -es

N.B: The -r and -o flags are optional. If not specified the script will annotate all the possible mutations for all the protein residues in a file called saturation_mutlist.txt. -es flag is optional as well and used in case the self mutations need to be removed

Example:

module load python/3.10/modulefile

bash do.sh UNIPROT_ID AA_RANGE

N.B. if you have more domains remember to concatenate the saturation lists with cat and call 
the final file saturation_mutlist.txt

#python saturation_mutlist.py -a P51587 -r 3415-3418 -o BRCA2_saturation_mutlist.txt -es
