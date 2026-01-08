#the files in the folder called embedded_residues.txt and  output_ppm3.txt  are just for illustrative purpouses 
#remember to replace them with the ones for your protein

#fill in also the info.txt file since important for FAIR principles  


### INPUT
# - OPM/PPM output csv: embedded_residues.txt
# - mutlist : 'mutlist.txt' #taken from cancermuts

#to create embedded_residues.txt to copy the output from ppm server to a txt file (see output file in the this folder)
#copy only the lines with the embedded residues info and add the header
#chain_ID       tot_res intervals

### Python load

module load python/3.10/modulefile

### Command

python filter_membranes.py -i embedded_residues.txt -m mutlist.txt -o filtered_mutlist_opm.txt -p

