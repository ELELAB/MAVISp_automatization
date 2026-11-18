#module load python/3.10/modulefile
uniprot=$1
aa=$2
#saturation list removing self mutations
python saturation_mutlist.py -a $uniprot -r $aa -o saturation_mutlist_$aa.txt -es
