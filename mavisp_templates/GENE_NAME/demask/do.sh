#!/bin/bash
uniprot=$1
source /usr/local/envs/demask/demask_env/bin/activate

python3 -m demask.homologs -s $uniprot -o myquery_homologs.a2m -c config.ini
python3 -m demask.predict -i myquery_homologs.a2m -o myquery_predictions.txt

