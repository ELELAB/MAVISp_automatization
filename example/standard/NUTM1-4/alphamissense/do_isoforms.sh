#!/bin/bash

# set ENST id of interest
enst=$1

# link the AlphaMissense isoform database
ln -s /data/databases/alphamissense/AlphaMissense_isoforms_aa_substitutions.tsv.gz .

# write header to output file
zcat AlphaMissense_isoforms_aa_substitutions.tsv.gz | head -n 4 | tail -n 1 > am.tsv

# find lines containing ENST id
zgrep -P "${enst}" AlphaMissense_isoforms_aa_substitutions.tsv.gz >> am.tsv

# compress output
gzip am.tsv