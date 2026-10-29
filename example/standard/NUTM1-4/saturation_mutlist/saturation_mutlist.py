import pandas as pd 
import re
import requests
import argparse


parser=argparse.ArgumentParser(description='saturation mutlist generator')

parser.add_argument("-a","--uniprot_ac",
	                    dest="uniprot_ac",
	                    type=str,
                        help="uniprot_ac of the protein which fasta"\
                             " sequence need to be downloaded")

parser.add_argument("-r","--residue_range",
	                    dest="residue_range", 
	                    type=str,
                        help="range of residues of interest")

parser.add_argument("-es","--exclude_self",
	                    dest="exclude_self",
	                    action = "store_true",
	                    help="exclude from the output file the self mutations")

parser.add_argument("-o","--output",
	                    dest="output",
	                    default="saturation_mutlist.txt",
	                    type=str,
                        help="range of residues of interest")

args=parser.parse_args()

def download_sequence_from_uniprot(uniprot_ac):
    url = f"https://www.uniprot.org/uniprot/{uniprot_ac}.fasta"
    response = requests.get(url)
    if response.status_code == 200:
        content = response.text
        return content
    else:
        print(f"Failed to download sequence for {uniprot_ac}. Exiting...")
        exit(1)

#------------------------ Download Fasta sequence --------------------------#

uniprot_ac = args.uniprot_ac
fasta=download_sequence_from_uniprot(uniprot_ac).split("\n")
fasta.pop(0)
sequence="".join(fasta)

#-------------------  Position-Residue dictionary creation -----------------#

residue_pos_dict={}

# In case of residue range keep only the region included in the provided range
if args.residue_range:
	# Check the format of the residue range
	pattern = r"(\d+)-(\d+)"
	match = re.match(pattern,args.residue_range)
	if not match:
		print("The residue range is written in the wrong format,"\
			  "please provide"\
		      " an argument in the following format:"\
			  "starting_residue-ending_residue i.e 2345-3987")
		exit(1)

	# Parse the residue range
	starting_residue=args.residue_range.split("-")[0]
	ending_residue=args.residue_range.split("-")[1]

	#compile the dictionary with the filtered residues
	for pos,res in enumerate(sequence):
		if pos in range(int(starting_residue)-1,int(ending_residue)):
			residue_pos_dict[pos]=res

	# Check the dictionary is consistent with the residue range provided
	if int(ending_residue)>list(residue_pos_dict.items())[-1][0]+1:
		print("The residue range provided contains residues positions"\
		      " not present in the protein. Exit...")
		exit(1)

# If not residue range create a dicitonary with all the protein residue
else:
	for pos,res in enumerate(sequence):
		residue_pos_dict[pos]=res

#------------------------------- Output ------------------------------------#

# List with the resiudes to be used to generate the mutations
mutation_list=["A","C","D","E","F","G","H","I","L","M","N","P",
               "Q","R","S","T","V","Y","W","K"]

# For each position-residue in the dictionary write the all the possible mut
# ation in a output file.
with open(args.output,"w") as f:
	for key,value in residue_pos_dict.items():
		tmp_mutation_list = mutation_list.copy()
		if args.exclude_self: 
			tmp_mutation_list.remove(value)
		for residue in tmp_mutation_list:
			mutation=value+str((int(key)+1))+residue
			f.write(mutation+"\n")





