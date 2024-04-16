#!/bin/sh
module load python/3.10

python3 clinvar.py -v variants.csv -o variants_output.csv 
