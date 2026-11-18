#running PDBminer:
tsp -N 1 bash do.sh GENE UNIPROT

#example
#tsp -N 1 bash do.sh SCP P22307

#NOTICE! 
#the run is setup to output a csv file. IF you need
#to use the data as more than a list of structures, 
#it is *MUCH* easier to handle the json file format,
#because json file format allows different datatypes 
#for each column, e.g. stings, lists, dictionary, 
#while the csv file format converts everything to strings, 
#which requires you to re-format everything when loading
#it into a dataframe.

#if you run with json (just remove -f csv), you import:

#import pandas as pd
#df = pd.read_json(uniprot_all.json)
