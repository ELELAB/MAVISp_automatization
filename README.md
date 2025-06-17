# Mavisp automatization pipeline

Cancer Structural Biology, Danish Cancer Society Research Center, 2100, Copenhagen, Denmark

Cancer Systems Biology, Health and Technology Department, Section for Bioinformatics, 2800, Lyngby, Denmark

## Introduction

MAVISp_automatization is a Snakemake pipeline that has been designed to automate the various steps involved in the mavisp framework. With this pipeline, users are able to perform the following tasks automatically: 
- genaration of metadata relative to each protein and the curator,
- mutations aggregation, 
- retrieval and trimming of PDB models from AplphaFold ,
- retrieval information of all available PDBs in the Protein Data Bank, 
- filtering of mutation lists for further calculations,
- retrieval information about interactors through two different databases,  
- protein phosphorylation prediction,
- RaSP calculation,
- classification of mutations occurring in protein sequence through different approaches,
- identification of denovo phosphorylation sites upon mutation
- early folding regions prediction
- relevant data collection to be used by mavisp.py

## Requirements

### Required software

The pipeline uses the following Python packages:

- pandas
- os
- re
- shutil
- glob
- datetime

The pipeline requires the following programs:
- cancermuts
- pdbminer
- rasp
- netphos
- demask
- RosettaDDGprediction
- efoldmine
- procheck

The pipeline requires the following databases:
- mentha
- AlphaFold2
- hpc-atlas
- Alphamissense
- AF_Huri

It also requires to have a full mutational scan performed using MutateX.

Before running the pipeline, please ensure that all the requirements are met.

### Input structure

The pipeline necessitates a CSV input file, a folder containing the scripts essential for the MAVISP framework (mavisp_templates), and a config.yaml file specifying the locations of the scripts essential for the pipeline. The input file, 'input.csv,' and the scripts directory should reside in the same directory as the Snakefile (an example is provided in the repository). The CSV file should encompass the subsequent comma-separated columns:

|protein|uniprot_ac|uniprot_id|trimmed|ref_seq|research_field|input_pdb|structure_source|model|curator_name|affiliation|
|-------|----------|----------|-------|-------|--------------|---------|----------------|-----|------------|-----------|
|BLM|P54132|BLM_HUMAN|1-359_368-1290|NP_000048|vus||AF2|v4|Elena Papaleo_Pablo Sanchez-Izquierdo|"DTU, Denmark; DCI, Denmark_DTU, Denmark; DCI, Denmark"|

where:
"protein" is the hugo name of the protein
"uniprot_id" is the uniprot id of the protein
"uniprot_ac" is the uniprot accession number of the protein
"trimmed" represents the residues range of your trimmed model (in case of multiple models specify the ranges "_" separeted)
"ref_seq" is the ref seq code associated with the isoform used by cancermuts (usually the first one). The RefSeq associated with the first isoform can be accessed through the Uniprot database in the "Sequence & isoforms" or "Sequence" section (depending on the protein being investigated) under the field "Sequence databases". The RefSeq code required by the pipeline consists of the "NP" followed by all subsequent characters until the dot, excluding it. There may be multiple RefSeq codes available. The one associated with the correct isoform includes the isoform code, which is enclosed in square brackets, corresponding to the canonical isoform. Verify by clicking on it that it corresponds to the correct isoform.
"research field" is the project or the research field in which the protein is involved (this is the name of the folder in which the cancermuts and rasp calculation will be organized).
"input_pdb" is an optional input pdb file the user provides
"structure_source" the source of the input structure, it needs to follow the nameing of this dictionary:
AFDB: "AlphaFold database",
AF3: "AlphaFold3 webserver",
AF2: "AlphaFold2 standalone",
PDB: "Experimental PDB structure",
Mod: "Homology model (PDB template, reconstruction)
"model" is the model of the structure used (vX)
An additional and optional input file in txt format containing the mutations from other sources than COSMIC, cBioPortal, and Clinvar can be provided. The mutations will be aggregated along with all the mutations found in the afromentioned databases by cancermuts in the metatable. This file name must start with the hugo name of the protein in capital letters and end with .txt extension (i.e BLM_my_mutation_list.txt); the mutations must be reported in one letter code as reported below without header
A34P
C56Y
Y678L
"curator_name" is the names of the curators that have worked in the curation of the protein (in case of multiple curators specify the names "_" separated).
"affiliation" is the organitation/s where the curator is affiliated to (in case of multiple curators specify the affiliations "_" separated in the same order as where introduced the curators). Put this whole field between comas ("affiliation1_affiliation2") to avoid parsing issues.

### Configfile.yaml
The pipeline supports a configuration file in the yaml format allowing customization of the script location and the environment settings required by the various software utilized in the pipeline. Below is an example of a configuration file:

```
input:
    path: input.csv
modules:
    interactome:
        mentha2pdb:
            database_date: 2025-04-28
            database_path: /data/databases/mentha-20250428
            AF_Huri_path:  /data/databases/AF_Huri_HuMAP/summary/HuRI.csv
            AF_HuMAP_path: /data/databases/AF_Huri_HuMAP/summary/humap.csv
            AF_Huri_HuMAP_path: /data/databases/AF_Huri_HuMAP
        hpc_atlas:
            database_path: /data/databases/hpc_atlas/HPC-Atlas_gene.txt
    structure_selection:
        alphafold:
            dssp_exec: /usr/local/dssp-4.4/bin/mkdssp
        pdbminer:
            pdbminer_env: . /usr/local/envs/PDBminer/bin/activate
        procheck:
            procheck_env: export prodir=/usr/local/procheck && export PATH=$PATH:$prodir
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
            saturation: true
            ELM_connections_per_run: 1
            folder_name: /data/raw_data/computational_data/cancermuts_testing/
        mutlist:
            source: . /usr/local/envs/cancermuts/bin/activate
    rasp: 
        conda_activation: /usr/local/miniconda3/bin/activate
        rasp_conda_env: /usr/local/envs/RaSP_workflow
        output_path_folder: /data/raw_data/computational_data/rasp_testing/
        folder_name: rasp
    rosetta_relax:
        rosetta_env: /usr/local/envs/rosettaddgprediction/bin/activate
        rosetta_module: /usr/local/rosetta-2022.11/
        rosetta_folder:  /data/raw_data/computational_data/rosetta_testing/
        rosetta_template: /data/raw_data/computational_data/rosetta_data/mavisp_templates/stability/ref2015_cartesian2020
    efoldmine:
        environment: set +eu && . /usr/local/envs/efoldmine/bin/activate && set -eu
```

Before running the pipeline, please ensure that you customize the configuration file by providing the correct paths for the environment, databases and the location of the cancermuts script template.
N.B the field ELM_connections_per_run specifiy the number of cancermuts run that can be run in parallell. It's essential to query ELM database once at time, for this reason the parameter needs to be set 1. 
N.B The template files for the RasP calculations are collected in the following path: 

```
/data/raw_data/computational_data/rasp_data/mavisp_templates/free/AF2_XX-YY/model_vX/
```

## Workflow

The pipeline automates the following steps for each entry in the input CSV file, collecting the corresponding readouts and organizing them into specific folders based on the HUGO name of each entry (refer to the example folder):

- **Retrieving AlphaFold models**: The associated PDB file for each entry in the input file is downloaded from the AlphaFold database and stored in the "structure_selection/alphafold_db/" 
  path (see example folder).

- **Retrieving available PDB information through PDB miner**: Information about all possible experimental structures available in the Protein Data Bank (PDB), such as resolution, 
  experimental method, and missing residues, is collected for each entry in the input file. The readouts are stored in the "structure_selection/pdbminer/" path.
  
- **Filtering PDB miner results**: Complexes extracted by PDB miner are filtered based on a 10 Å interaction distance. Only complexes meeting this criterion are retained, and their binding interface residues are extracted. If available, the PDB structures of these complexes are downloaded. The readouts are stored in the "structure_selection/pdbminer_complexes/" folder. 

- **Trimming AlphaFold models**: The AlphaFold model is trimmed based on the specified range in the input file. Only regions with high pLDDT scores are retained. The trimmed PDB files are 
  stored in the "structure_selection/trimmed_model/" path. The residue range is used in subsequent steps to filter the mutation list for calculations.

- **Running Procheck for Structure Quality Assessment**: The Procheck tool is used to assess the quality of the trimmed AlphaFold models. This analysis generates summary files that provide detailed information about the stereochemical quality of the trimmed models generated in the previous step, including Ramachandran plot statistics. The readouts are stored in the "structure_selection/procheck/" folder.

- **Retrieving mutations from the ClinVar database**: Missense mutations reported in the ClinVar database, along with corresponding classifications, review status, and associated 
  conditions (diseases), are collected for each entry in the input file. The readouts are stored in the "clinvar_gene/" folder.

- **Cancermuts step**: The mutations obtained from ClinVar are converted into a supported format for the "cancermuts" tool and provided as input, along with any external mutation lists 
  provided by the user. The pipeline aggregates cancer mutations from cBioPortal, COSMIC, and ClinVar databases, along with structural information (SLiMs, secondary structure, Post Translational Modification sites PTMs) of the protein. In case multiple entrez IDs are found for one protein, no gene ID will be assigned and mutation retrieval from cBioPortal will be skipped. The output is stored in a path specified in the config.yaml file.

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

- **RosettaDDGprediction relax step**: The Cartesian2020 protocol and the ref2015 energy function from the Rosetta suite are exploited to perform the relax step required by the software using the trimmed structure obtained in the step above. The output is stored in a path specified in the config.yaml file.

- **identification of denovo phosphorylation sites**: netphos is used on wild-type or mutant sequences to identify cases in which mutations cause significant changes in the propensity of a certain site to be phosphorylated

- **Prediction of early folding sites**: EFoldMine is employed in order to predict regions with early folding propensity from the target protein's primary sequence, so that it can be subsequently identified whether the mutation sites of the investigated variants fall within the predicted early folding regions. The EFoldMine output is stored in the "efoldmine/" folder.

- **Metadata generation step**: All the relevant information of the curation of each protein for each run is stored in "metadata" folder.

- **Output collection**: The most relevant outputs from all the steps are collected in one single folder that will be further used to automatically generate the mavisp.csv. The output is stored in the "simple_mode" folder.  

## Output structure

A folder is created for each entry, named after the entry's HUGO name, and it collects all the data associated with that protein. The CancerMuts data, the Rosetta relax step output and the RaSP data are stored separately in a location specified in the configuration file.

The output structure is based on the first entry from the input file, but it follows a similar pattern for each entry. It's important to note that the folders containing the CancerMuts data, the relax step from Rosetta and the RaSP data are located outside the BLM folder, as mentioned earlier.
```
BLM
├── alphamissense
│   ├── AlphaMissense_aa_substitutions.tsv.gz
│   ├── am.tsv.gz
│   ├── am_decompressed.tsv
│   ├── am_head.tsv
│   ├── do.sh
│   └── readme.txt
│
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
│
├── clinvar
│   ├── clinvar.py
│   ├── entry_not_found_variants.csv
│   ├── readme.txt
│   ├── variants.csv
│   ├── variants_output.csv
│   └── variants_to_check.csv
│
├── clinvar_gene
│   ├── BLM_mutation_list.txt
│   ├── clinvar.py
│   ├── genes.csv
│   ├── readme.txt
│   └── variants_output.csv
│
├── demask
│   ├── P54132.blast.json
│   ├── P54132.fasta
│   ├── config.ini
│   ├── myquery_homologs.a2m
│   ├── myquery_predictions.txt
│   └── readme.txt
│
├── external_mutation_lists
│   ├── BLM_icope_sention_input.txt
│   └── BLM_icope_sention_input2.txt
│
├── interactome
│   ├── hpc_atlas
│   │   ├── BLM.out
│   │   ├── HPC-Atlas_gene.txt → /data/databases/hpc_atlas/HPC-Atlas_gene.txt
│   │   └── readme.txt
│   └── mentha2pdb
│       ├── 2025-04-28 → /data/databases/mentha-20250428/2025-04-28
│       ├── P54132.csv
│       ├── inputs_afmulti
│       │   └── BLM
│       │       ├── AIPL1
│       │       │   └── input.fasta
│       │       ├── APBB1
│       │       │   └── input.fasta
│       │       └── … (etc.)
│       ├── mentha2pdb.py
│       ├── readme.txt
│       └── target_uniprot_ac.txt
│
├── netphos
│   ├── P54132.fasta
│   ├── netphos.out
│   └── readme.txt
│
├── structure_selection
│   ├── original_model
│   │   ├── BLM
│   │   │   ├── P54132.csv
│   │   │   ├── P54132.json
│   │   │   ├── P54132.pdb
│   │   │   ├── P54132.yaml
│   │   │   └── P54132_dssp.out
│   │   ├── config_alphafolddb.yaml
│   │   ├── get_alphafolddb_data.py
│   │   ├── readme.txt
│   │   └── regions_pLDDT.csv
│   ├── domain_annotations
│   │   ├── domains_mutlist.csv
│   │   ├── get_domains.py
│   │   ├── mutlist.txt → ../../cancermuts/mutlist_08092023.txt
│   │   ├── readme.md
│   │   └── summary.csv
│   ├── efoldmine
│   │   ├── readme.md
│   │   ├── P54132.fasta
│   │   └── P54132.tabular
│   ├── pdbminer
│   │   ├── input_file.csv
│   │   ├── log.txt
│   │   ├── readme.txt
│   │   └── results
│   │       └── P54132
│   │           └── P54132_all.csv
│   ├── pdbminer_complexes
│   │   ├── readme.txt
│   │   ├── find_PDBminer_complexes.py
│   │   ├── P54132_filtered.csv
│   │   └── P54132_pdb_complexes
│   │      ├── 7XUW.pdb
│   │      └── 7XV0.pdb
│   ├── procheck
│   │   ├── P54132_368-1290.pdb
│   │   ├── anglen.log
│   │   ├── P54132_368-1290_*.ps
│   │   ├── P54132_368-1290.new
│   │   ├── P54132_368-1290.rin
│   │   ├── pplot.log
│   │   ├── bplot.log
│   │   ├── P54132_368-1290.out
│   │   ├── P54132_368-1290.sco
│   │   ├── clean.log
│   │   ├── P54132_368-1290.lan
│   │   ├── P54132_368-1290.sdh
│   │   ├── secstr.log
│   │   ├── nb.log
│   │   ├── P54132_368-1290.nb
│   │   ├── P54132_368-1290.pln
│   │   ├── P54132_368-1290.sum
│   │   ├── tplot.log
│   │   └── readme.txt
│   └── trimmed_model
│       ├── P54132_368-1290.pdb
│       ├── filtre_pdb.py
│       └── readme.txt
│
├── metadata
│   ├── metadata.yaml
│   └── importing.yaml
│
└── simple_mode
    ├── clinvar
    │   └── variants_output.csv
    │
    ├── ptm
    │   ├── summary_stability.txt
    │   ├── sasa.rsa
    │   └── metatable.csv
    │
    ├── sas
    │   └── sasa.rsa
    │
    ├── demask
    │   └── myquery_predictions.txt
    │
    ├── alphamissense
    │   └── am.tsv.gz
    │
    ├── cancermuts
    │   └── metatable_pancancer_BLM.csv
    │
    ├── efoldmine
    │   └── P54132.tabular
    │
    ├── metadata.yaml
    │
    ├── pfam
    │   └── summary.csv
    │
    ├── alphafold
    │   └── P54132.csv
    │
    ├── stability
    │   └── AF2_36-104
    │   │   └── AF2
    │   │       └── model_v3
    │   │           ├── rasp
    │   │           │   └── post_processed.csv
    │   │           └── foldx5
    │   │               └── energies.csv
    │   └── AF2_105-200
    │       └── AF2
    │           └── model_v3
    │               ├── rasp
    │               │   └── post_processed.csv
    │               └── foldx5
    │                   └── energies.csv
    │
    └── mutation_list
        └── mutations_pmid_22052025.txt

```

## Usage

In order to run the pipeline (see example folder):

conda deactivate

module load python/3.10

snakemake -c 1





