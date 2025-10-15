source /usr/local/envs/py310/bin/activate

python get_domains.py -u $1 -m $2
python ted_api.py -u $1
#$1 uniprot accession number
#$2 mutationlist
