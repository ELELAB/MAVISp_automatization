# Get mutlists  
*Last updated*: 09/09/2024

## Description

The script `get_mutlists.py` takes in input the cancermut's metatable and, depending on the commandline flags, returns one or more `mutlist.txt` files that will contain the mutations of the cancermuts metatable found with a source (i.e., COSMIC, cBioPortal, and/or manual annotation) and occurring at positions covered by the selected structure. The scripts also performs a check of the compatibility of WT residues between the metatable and the pdb files selected for subsequent analyses.

## Requirements

- Python >= 3.7 
- Pandas
- Biopython

## Input
- cancermuts metatable
- get_mutlists.py

## Usage

`python get_mutlists.py [-h] -m METATABLE -d DOMAIN [DOMAIN ...] -p PDBFILE [PDBFILE ...] [-ch CHAIN] [-M] [-R] [-H] [-C] [-D DAY]`

### Required flags
- `-m` : name of the metatable
- `-d` : domain(s) residues (see below)
- `-p` : pdb file(s) 

If you have the full-length protein or a trimmed model in which only the N-ter and/or C-ter portions have been removed you can follow the single domain example (i.e., it is enough to specify the range of residue after the `-d` flag: `-d 1:756`).  
- `python get_mutlists.py -m metatable_pancancer_MLH1.csv -d 1:756` 

If you have trimmed the structure in multiple domains, follow the multiple domains example. For instance, if your structures are covering two ranges of residues (e.g., 1:341 and 501:756), you have to specify both the ranges (space-separated) after the `-d` flag: `-d 1:341 501:756`.
- `python get_mutlists.py -m metatable_pancancer_MLH1.csv -d 1:341 501:756`

With only the required flags you will get the mutations list in the one-letter format (e.g., `R13L`) and the ELM mutlist.

### Optional flags
If you wish to disable the sanity check between the PDB file and the metatable, or to generate mutation lists in formats compatible with MutateX, Rosetta, HGVS, and CABSFLEX, or filtered to include only mutations curated by the user and those reported in the COSMIC, ClinVar, and cBioPortal databases, you must specify the appropriate combination of optional flags:
- `-s` : remove the sanity check on the consistency of positions and residues between the provided PDB file and the input metatable.
- `-M` : to get MutateX mutations list format and a filtered version that only includes the phosphorylated variants (Mutatex compatible) 
- `-R` : to get Rosetta mutations list format.  
- `-H` : to get hgvs mutations list format 
- `-S` : Generate a filtered mutation list excluding those only present in the saturation mutagenesis list
- `-C` : to get CABSFLEX mutations list format 
- `-ch` : chain label (default = 'A')
- `-D` : date (DDMMYYYY) that is added to the output file names. (default = day of the cancermuts run)

Notes:  

The `-ch` flag is required to get the MutateX, Rosetta and CABSFLEX compatible formats, it has been set by default to **A**.   

Please, be aware that you need to type one pdb name (flag `-p`) for each domain you specify (see `-d`).  

If the structure is trimmed in two or more domains, the **rosetta mutlist** will automatically be split in multiple files containing only the mutations of a certain domain; additionally, the mutlist with all the mutations in the rosetta-compatible format is returned.

### Example
See the `example` directory and the `do.sh` script within.

## Output
The script returns one or multiple `.txt` file(s) containing mutations lists in the following format(s):
- `R13L`
- `RA13L` (MutateX compatible)
- `TA56C` and `TA56p` (MutateX compatible - phospho version)
- `A.R.13.L A` (Rosetta compatible)
- `p.Ala13Leu` (hgvs format)
- `4P7A.pdb A 13 ARG LEU` (CABSFLEX compatible)
- `R13L, "XXX binding motif"` (ELM mutlist)
