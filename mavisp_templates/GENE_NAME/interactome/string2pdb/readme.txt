string2pdb: STRING Interactor Retrieval

This script is used in the MAVISP Snakemake pipeline to retrieve interactors from the STRING database for a given UniProt accession (UPAC).

Usage:
    ./string2pdb {uniprot_ac}

Default settings:
- STRING physical subnetwork is queried.
- Minimum interaction score: 0.15
- Only interactions supported by either curated database annotations (database score > 0) or experimental data (experimental score > 0) are included.
- Available PDBs for retrieved complexes are also identified and included.

The script is automatically called by Snakemake using these defaults. Output is saved to:
    interactome/string2pdb/{uniprot_ac}_string_interactors.csv
