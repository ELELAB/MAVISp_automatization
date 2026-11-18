#FOR MAVISp run with mutation list
#Requirements in the folder
yaml files 
res_list.txt
aggregate2mutatex.sh
mutlist_rosetta.txt (from cancermuts folder in project folder)


#Activate virtual environment before running on tsp or run on HPC resources

####on bioinfo01, 03 or 04
#source /usr/local/envs/rosettaddgprediction/bin/activate


#load module
module load rosetta/2022.11

tsp -N 1 rosetta_ddg_run -p file.pdb  -cr relax2020_ref2015.yaml -n 1 -r /usr/local/rosetta-2022.11/ -cs rosettampi.yaml
#NB if after this step you run ddg you need to use the relaxed PDB structure and cartesian2020_ref2015.yaml

 tsp -N 8 rosetta_ddg_run -p relax/relax_P01116_1-176_0001.pdb -cr cartesian2020_ref2015.yaml -n 8 -r /usr/local/rosetta-2022.11/ -cs rosettampi.yaml -l mutlist_rosetta.txt 


#to check for crashes during the run
rosetta_ddg_check_run -cr cartesian2020_ref2015.yaml -d cartesian -mf cartesian/mutinfo.txt

#at the end of the run 
#run aggregate using tsp 
tsp -N 8 sh aggregate2mutatex.sh

#final list file name in aggregate/:
ddg_mutations_aggregate.csv
ddg_mutations_structures.csv

