filtre_pdb.py is a script to trim pdb strucutres based on a residue range:
The script is suited for working on pdb with only a chain (in case of multiple chains
the script will generate a new pdb file containing only the chain that need to be 
trimmed).  

The script takes the following parameters as arguments:

1) file.pdb --> input pdb file with the protein to trim

2) chain_id --> chain_id of the chain on which perform the trimming

3) start    --> first residue of the trimmed protein

4) end      --> last residue of the trimmed protein

5) output   --> file name of the trimmed pdb file

The script generates a new pdb file named as specified in the 5) argument containing
the residues in the range between the specified start and end arguments (3-4).

Example:

module load python/3.10/modulefile
python filter_pdb.py 6TNYab.pdb B 1 38 6TNYab_1-38.pdb 



