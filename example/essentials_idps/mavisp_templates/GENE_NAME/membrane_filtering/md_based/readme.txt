#To do the analysis:

### INPUT
# - simulation analysis  output csv: filt_res_lipid_surround.csv
# - mutlist : 'mutlist.txt' #taken from cancermuts

### Python load

module load python/3.10/modulefile

### Command

python filter_membranes.py -i filt_res_lipid_surround.csv -m mutlist.txt -o filtered_mutlist_md.txt -p

