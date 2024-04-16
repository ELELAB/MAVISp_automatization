pdb=$1
list=$2

. /usr/local/envs/mutatex/bin/activate

#label file from PDB
pdb2labels -p $pdb
ddg2summary -p $pdb -d final_averages -l mutation_list.txt -L $list
mv summary.txt summary_stability.txt


