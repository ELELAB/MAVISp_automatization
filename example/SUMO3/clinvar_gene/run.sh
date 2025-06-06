module load python/3.10
python3 clinvar.py -g gene.csv -o genes_output.csv
cp genes_output.csv variants_output.csv
