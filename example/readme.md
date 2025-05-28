# Mavisp automatization pipeline

Cancer Structural Biology, Danish Cancer Society Research Center, 2100, Copenhagen, Denmark

Cancer Systems Biology, Health and Technology Department, Section for Bioinformatics, 2800, Lyngby, Denmark

## Introduction

"ToUpdate" is a Snakemake pipeline that has been designed to automate the various steps involved in the mavisp framework. With this pipeline, users are able to perform the following tasks automatically: 
- mutations aggregation, 
- retrieval and trimming of PDB models from AplphaFold ,
- retrieval information of all available PDBs in the Protein Data Bank, 
- filtering of mutation lists for further calculations,
- retrieval information about interactors through two different databases,  
- protein phosphorylation prediction,
- RaSP calculation,
- classification of mutations occurring in protein sequence through different approaches

## Requirements

### Required software

The pipeline uses the following Python packages:

- pandas
- os
- re
- shutil

The pipeline requires the following programs:
- cancermuts
- pdbminer
- rasp
- netphos
- demask

The pipeline requires the following databases:
- mentha
- AlphaFold2
- hpc-atlas
- Alphamissense
- AF_Huri

Before running the pipeline, please ensure that all the requirements are met.

### Input structure

The pipeline necessitates a CSV input file, a folder containing the scripts essential for the MAVISP framework (mavisp_templates), and a config.yaml file specifying the locations of the scripts essential for the pipeline. The input file, 'input.csv,' and the scripts directory should reside in the same directory as the Snakefile (an example is provided in the repository). The CSV file should encompass the subsequent comma-separated columns:

|protein|uniprot_ac|uniprot_id|trimmed|ref_seq|research_field|
|-------|----------|----------|-------|-------|--------------|
|BLM|P54132|BLM_HUMAN|1-359_368-1290|NP_000048|vus|

where:
"protein" is the hugo name of the protein
"uniprot_id" is the uniprot id of the protein
"uniprot_ac" is the uniprot accession number of the protein
"trimmed" represents the residues range of your trimmed model (in case of multiple models specify the ranges "_" separeted)
"ref_seq" is the ref seq code associated with the isoform used by cancermuts (usually the first one). The RefSeq associated with the first isoform can be accessed through the Uniprot database in the "Sequence & isoforms" or "Sequence" section (depending on the protein being investigated) under the field "Sequence databases". The RefSeq code required by the pipeline consists of the "NP" followed by all subsequent characters until the dot, excluding it. There may be multiple RefSeq codes available. The one associated with the correct isoform includes the isoform code, which is enclosed in square brackets, corresponding to the canonical isoform. Verify by clicking on it that it corresponds to the correct isoform.
"research field" is the project or the research field in which the protein is involved (this is the name of the folder in which the cancermuts and rasp calculation will be organized).
An additional and optional input file in txt format containing the mutations from other sources than COSMIC, cBioPortal, and Clinvar can be provided. The mutations will be aggregated along with all the mutations found in the afromentioned databases by cancermuts in the metatable. This file name must start with the hugo name of the protein in capital letters and end with .txt extension (i.e BLM_my_mutation_list.txt); the mutations must be reported in one letter code as reported below without header
A34P
C56Y
Y678L

### Configfile.yaml
The pipeline supports a configuration file in the yaml format allowing customization of the script location and the environment settings required by the various software utilized in the pipeline. Below is an example of a configuration file:

```
input:
    path: input.csv
modules:
    interactome:
        mentha2pdb:
            database_path: /data/databases/mentha-20230417/2023-04-17
            AF_Huri_path:  /data/databases/AF_Huri_HuMAP/summary/HuRI.csv
            AF_HuMAP_path: /data/databases/AF_Huri_HuMAP/summary/humap.csv
            AF_Huri_HuMAP_path: /data/databases/AF_Huri_HuMAP
        hpc_atlas:
            database_path: /data/databases/hpc_atlas/HPC-Atlas_gene.txt
    structure_selection:
        alphafold:
            dssp_exec: /usr/local/dssp-3.0.10/bin/mkdssp
        pdbminer:
            conda_activation: /usr/local/miniconda3/bin/activate
            pdbminer_conda_env: /usr/local/envs/PDBminer
    mutations_classifier:
        demask:
            source: /usr/local/envs/demask/demask_env/bin/activate
            database_path: /data/databases/uniref90demask/uniref90demask
            blastp_path: /usr/local/blast_plus-2.13.0/bin/blastp
        alphamissense:
            database_path: /data/databases/alphamissense/AlphaMissense_aa_substitutions.tsv.gz
    mutations_aggregation:
        cancermuts:
            source: /usr/local/envs/cancermuts/bin/activate
            Cosmic_database_path: /data/databases/cosmic-v96/CosmicMutantExport.tsv
            template_path: /data/raw_data/computational_data/cancermuts_testing/mavisp_templates_final/pancancer.py
            ELM_connections_per_run: 1
            folder_name: /data/raw_data/computational_data/cancermuts_testing/
    rasp:
        conda_activation: /usr/local/miniconda3/bin/activate
        rasp_conda_env: /usr/local/envs/RaSP_workflow
        output_path_folder: /data/raw_data/computational_data/rasp_testing
        folder_name: rasp
```

Before running the pipeline, please ensure that you customize the configuration file by providing the correct paths for the environment, databases and the location of the cancermuts script template.
N.B the field ELM_connections_per_run specifiy the number of cancermuts run that can be run in parallell. It's essential to query ELM database once at time, for this reason the parameter needs to be set 1. 
N.B The template files for the RasP calculations are collecte din the follwoing path: 

```
/data/raw_data/computational_data/rasp_data/mavisp_templates/free/AF2_XX-YY/model_vX/
```

## Workflow

The pipeline automates the following steps for each entry in the input CSV file, collecting the corresponding readouts and organizing them into specific folders based on the HUGO name of each entry (refer to the example folder):

- **Retrieving AlphaFold models**: The associated PDB file for each entry in the input file is downloaded from the AlphaFold database and stored in the "structure_selection/alphafold_db/" 
  path (see example folder).

- **Retrieving available PDB information through PDB miner**: Information about all possible experimental structures available in the Protein Data Bank (PDB), such as resolution, 
  experimental method, and missing residues, is collected for each entry in the input file. The readouts are stored in the "structure_selection/pdbminer/" path.

- **Trimming AlphaFold models**: The AlphaFold model is trimmed based on the specified range in the input file. Only regions with high pLDDT scores are retained. The trimmed PDB files are 
  stored in the "structure_selection/trimmed_model/" path. The residue range is used in subsequent steps to filter the mutation list for calculations.

- **Retrieving mutations from the ClinVar database**: Missense mutations reported in the ClinVar database, along with corresponding classifications, review status, and associated 
  conditions (diseases), are collected for each entry in the input file. The readouts are stored in the "clinvar_gene/" folder.

- **Cancermuts step**: The mutations obtained from ClinVar are converted into a supported format for the "cancermuts" tool and provided as input, along with any external mutation lists 
  provided by the user. The pipeline aggregates cancer mutations from cBioPortal, COSMIC, and ClinVar databases, along with structural information (SLiMs, secondary structure, Post Translational Modification sites PTMs) of the protein. The output is stored in the same folder as the "cancermuts_data" Snakefile.

- **Creation of MutateX, Rosetta, and CabsFlex compatible mutation lists**: Using the "metatable.csv" file from the cancermuts output and the residue range from the trimmed model, 
  mutation lists are created in the appropriate format for the following calculations: protein stability through the FoldX energy function (MutateX), protein stability through the Rosetta energy function Ref2015, and conformational structural ensemble generation of mutants and WT through CabsFlex (used with the ensemble mode). Additionally, a mutation list with potential phospho residues for the MutateX calculation is created. Finally, a mutation list with mutations expressed in one-letter code is created for other annotation steps. The files are stored in the "cancermuts" folder.

- **Domain annotations based on the mutation list**: The mutation list obtained in the previous step is used to retrieve information about the domains in the protein under investigation. 
  The readouts are stored in the "structure_selection/domain_annotations/" path.

- **ClinVar annotations retrieval for the mutation list**: The same mutation list used in the previous step is used to retrieve corresponding classifications, conditions, review statuses, 
  and phenotypes from the ClinVar database. The output is stored in the "clinvar" folder.

- **Retrieving interactors from the Mentha database**: Information about experimentally validated interactors collected in the Mentha and database is retrieved and organized in a CSV file, 
  along with information about possible PDB IDs and associated PMID references. Each interactor is assigned a score indicating the certainty of the interaction. The output is stored in the "interactome/mentha2pdb" path.

- **Retrieving interactors from the hpc-atlas database**: Interactor information collected in the hpc-atlas database is retrieved and organized into a txt file. The file is stored in the 
  "interactome/hpc_atlas" path.

- **Mutation classification using the Demask approach**: The uniprot fasta file, automatically retrieved by the pipeline, is provided to the Demask software to obtain a classification for 
  every possible mutation at each sequence position. The output is stored in the "demask" folder.

- **Mutation classification using the Alphamissense approach**: The uniprot_ac, automatically retrieved by the pipeline, is used to filter the Alphamissense file to obtain a classification for 
  every possible mutation at each sequence position. The output is stored in the "alphamissense" folder.

- **Rasp calculation on trimmed models**: The trimmed models obtained in the previous step are provided to the Rasp workflow to obtain stability effect predictions for every possible 
  mutation at each sequence position. The output is stored in the same folder as the "rasp_data" Snakefile.

- **Protein phosphorylation predictions using NetPhos**: The NetPhos software is used to predict phosphorylatable residues in the protein sequence, along with the corresponding kinases 
  responsible for the phosphorylation and their associated confidence scores. The output is stored in the "netphos" folder.


## Output structure

A folder is created for each entry, named after the entry's HUGO name, and it collects all the data associated with that protein. The CancerMuts data is stored separately in the following path: /data/raw_data/computational_data/cancermuts_testing/. The RaSP data is stored in a location specified in the configuration file.

The output structure is based on the first entry from the input file, but it follows a similar pattern for each entry. It's important to note that the folders containing the CancerMuts data and RaSP data are located outside the BLM folder, as mentioned earlier.
```
├── cancermuts
│   ├── P54132_368-1290.pdb
│   ├── get_mutlists.py
│   ├── mutlist_08092023.txt
│   ├── mutlist_ELM_08092023.txt
│   ├── mutlist_hgvs_08092023.txt
│   ├── mutlist_mutatex_08092023.txt
│   ├── mutlist_mutatex_P_08092023.txt
│   ├── mutlist_rosetta_08092023.txt
│   └── readme.txt
├── clinvar
│   ├── clinvar.py
│   ├── entry_not_found_variants.csv
│   ├── readme.txt
│   ├── variants.csv
│   ├── variants_output.csv
│   └── variants_to_check.csv
├── clinvar_gene
│   ├── BLM_mutation_list.txt
│   ├── clinvar.py
│   ├── genes.csv
│   ├── readme.txt
│   └── variants_output.csv
├── demask
│   ├── P54132.blast.json
│   ├── P54132.fasta
│   ├── config.ini
│   ├── myquery_homologs.a2m
│   ├── myquery_predictions.txt
│   └── readme.txt
├── external_mutation_lists
│   ├── BLM_icope_sention_input.txt
│   └── BLM_icope_sention_input2.txt
├── interactome
│   ├── hpc_atlas
│   │   ├── BLM.out
│   │   ├── HPC-Atlas_gene.txt -> /data/databases/hpc_atlas/HPC-Atlas_gene.txt
│   │   └── readme.txt
│   └── mentha2pdb
│       ├── 2023-04-17 -> /data/databases/mentha-20230417/2023-04-17
│       ├── P54132.csv
│       ├── inputs_afmulti
│       │   └── BLM
│       │       ├── AIPL1
│       │       │   └── input.fasta
│       │       ├── APBB1
│       │       │   └── input.fasta
│       │       ├── APITD1
│       │       │   └── input.fasta
│       │       ├── ATM
│       │      ...
│       ├── mentha2pdb.py
│       ├── readme.txt
│       └── target_uniprot_ac.txt
├── netphos
│   ├── P54132.fasta
│   ├── netphos.out
│   └── readme.txt
└── structure_selection
    ├── alphafold_db
    │   ├── BLM
    │   │   ├── P54132.csv
    │   │   ├── P54132.json
    │   │   ├── P54132.pdb
    │   │   ├── P54132.yaml
    │   │   └── P54132_dssp.out
    │   ├── config_alphafolddb.yaml
    │   ├── get_alphafolddb_data.py
    │   ├── readme.txt
    │   └── regions_pLDDT.csv
    ├── domain_annotations
    │   ├── domains_mutlist.csv
    │   ├── get_domains.py
    │   ├── mutlist.txt -> ../../cancermuts/mutlist_08092023.txt
    │   ├── readme.md
    │   └── summary.csv
    ├── pdbminer
    │   ├── input_file.csv
    │   ├── log.txt
    │   ├── readme.txt
    │   └── results
    │       └── P54132
    │           └── P54132_all.csv
    └── trimmed_model
        ├── P54132_368-1290.pdb
        ├── filtre_pdb.py
        └── readme.txt


```
## Usage

In order to run the pipeline (see example folder):

conda deactivate

module load python/3.10/modulefile

tsp -N 4 snakemake -s Snakefile -c 4





