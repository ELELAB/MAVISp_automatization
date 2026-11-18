# set UniProt AC of interest in this variable - to be changed
upac=$1
enst=$2

# link the AF aminoacids substitution scores files in this directory
canonical_db="/data/databases/alphamissense/AlphaMissense_aa_substitutions.tsv.gz"
isoform_db="/data/databases/alphamissense/AlphaMissense_isoforms_aa_substitutions.tsv.gz"

ln -snf "$canonical_db" .
ln -snf "$isoform_db" .

# write header to output file
zcat "$canonical_db" | head -n 4 | tail -n 1 > am.tsv

# find lines containing UniProt AC of interest and save them to file
if zgrep -P "${upac}\t" "$canonical_db" >> am.tsv; then
    echo "AlphaMissense predictions found for UniProt accession ${upac}"
else
    echo "No AlphaMissense predictions found for ${upac}; trying transcript fallback"

    if [ -z "$enst" ]; then
        echo "WARNING: no Ensembl transcript ID available for fallback"
    elif zgrep -P "${enst}" "$isoform_db" >> am.tsv; then
        echo "AlphaMissense predictions found for transcript ${enst}"
    else
        echo "WARNING: no AlphaMissense predictions found for transcript ${enst}"
    fi
fi

# compress output file
gzip am.tsv
