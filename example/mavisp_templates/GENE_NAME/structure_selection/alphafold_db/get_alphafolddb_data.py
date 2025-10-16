#!/usr/bin/env python

# Standard library
import itertools
import os
# Third-party packages
import Bio.PDB as bpdb
from Bio.PDB import PDBParser
from Bio.PDB.DSSP import DSSP

import pandas as pd
import yaml


def get_dssp_dataframe(pdb_file, dssp_location):
    
    #create empty lists
    chain_list = []
    pos_list = []
    ss_list = []
    
    #get the model
    p = PDBParser()
    structure = p.get_structure("pdb", pdb_file)
    model = structure[0]
    
    #create the dssp dictionary:
    dssp = DSSP(model, pdb_file, dssp=dssp_location)

    #get secondary structure for each residue
    for residue in list(dssp.keys()):
        chain_list.append(residue[0])
        pos_list.append(residue[1][1])
        ss = dssp[residue][2]
        ss_list.append("-" if ss == "P" else ss)
    
    #repare "P" with "-": 
    for i in range(len(ss_list)):
        if ss_list[i] == "P":
            ss_list[i] = "-"
            
    #create df:
    df = pd.DataFrame({"chain":chain_list,"resnum":pos_list,"secstruc":ss_list})
    
    return df


def get_plddt_per_residue(pdb_file):
    """Get per-residue pLDDT scores from the AlphaFold PDB file.
    """

    # Create an empty list to store the data for each
    # residue
    data = []

    # Get the structure
    structure = bpdb.PDBParser().get_structure("struct", pdb_file)

    # For each chain in the first model (only one model in
    # AlphaFold structures)
    for idx, residue in enumerate(structure[0].get_residues()):

        # Get the residue name
        resname = residue.get_resname()

        # Get only the first atom of the current residue
        # (we are only interest in residue-level properties)
        first_atom = list(residue.get_atoms())[0]

        # Get properties of interest from the ID of the first atom
        struct, mod, chain, (hetflag, resnum, icode), _ = \
            first_atom.get_full_id()
        
        # Get the B-factor for the first atom (= pLDDT score of
        # the residue)
        bfactor = first_atom.get_bfactor()

        # Add a dictionary containing the data for the current
        # residue to the final list
        data.append({"chain" : chain,
                     "resnum" : resnum,
                     "resname" : resname,
                     "pLDDT" : bfactor})

    # Create a data frame with the per-residue data and return it
    return pd.DataFrame(data)


def get_regions_by_plddt(plddt_df, plddt_cutoff):
    """Separate the protein regions having a pLDDT score greater than
    or equal to a pre-defined cut-off from those having a pLDDT score
    lower than the cut-off, and report them in ranges.
    """

    # Define a helper function to define ranges
    def get_ranges(i):
        ranges = []
        for a, b in itertools.groupby(enumerate(i), \
                                      lambda pair: pair[1] - pair[0]):
            b = list(b)
            ranges.append((b[0][1], b[-1][1]))
        return ranges

    # Keep only those residues with pLDDT score >= cut-off
    plddt_greaterequal = \
        plddt_df[plddt_df["pLDDT"] >= plddt_cutoff]["resnum"].tolist()

    # Get the protein regions covered by those residues as ranges
    plddt_greaterequal_ranges = get_ranges(plddt_greaterequal)

    # Convert the ranges whose beginning and end coincide (region
    # with only one residue) to length-1 tuples
    plddt_greaterequal_ranges = \
        [f"{b}-{e}" if b != e else str(b) for (b, e) \
         in plddt_greaterequal_ranges]

    # Keep only those residues with pLDDT score < cut-off
    plddt_lower = \
        plddt_df[plddt_df["pLDDT"] < plddt_cutoff]["resnum"].tolist()

    # Get the protein regions covered by those residues as ranges
    plddt_lower_ranges = get_ranges(plddt_lower)

    # Convert the ranges whose beginning and end coincide (region
    # with only one residue) to length-1 tuples
    plddt_lower_ranges = \
        [f"{b}-{e}" if b != e else str(b) for (b, e) \
         in plddt_lower_ranges]

    # Return the ranges
    return plddt_greaterequal_ranges, plddt_lower_ranges


def main(config,
         wd):

    # Set an empty list to store data for the final data frame
    data = []

    # Get the cut-off for the pLDDT score
    plddt_cutoff = config["plddt_cutoff"]
    
    # For each UniProt ID
    for uniprot_id, protein_data in config["uniprot_ids"].items():


        #----------------------- Local PDB ------------------------#


        # Create the directory to store the data
        # about the protein associated with the
        # UniProt ID
        protein_dir_name = protein_data["dir_name"]
        protein_dir = os.path.join(wd, protein_dir_name)
        os.makedirs(protein_dir, exist_ok = True)

        # Expect a local PDB placed beforehand at <wd>/<dir_name>/<uniprot_id>.pdb
        out_pdb_file = os.path.join(protein_dir, f"{uniprot_id}.pdb")

        if not os.path.exists(out_pdb_file):
            raise FileNotFoundError(
                f"Expected local PDB not found: {out_pdb_file}\n"
                f"Place a file named {uniprot_id}.pdb in {protein_dir} before running."
            )


        #--------------------------- DSSP ----------------------------#

        # Get DSSP df for the AlphaFold PDB file        
        df = get_dssp_dataframe(out_pdb_file, config["dssp_exec"])

        #--------------------------- pLDDT ---------------------------#


        # Get the per-residue pLDDT score
        plddt_df = get_plddt_per_residue(out_pdb_file)

        # Get the protein regions with pLDDT score greater than/equal
        # to or lower than the pre-defined pLDDT cut-off
        plddt_greaterequal_ranges, plddt_lower_ranges = \
            get_regions_by_plddt(plddt_df, plddt_cutoff)

        # Append the residue ranges representing these regions
        # to the final list
        data.append(\
                {"uniprot_id" : uniprot_id,
                 "dir_name" : protein_dir_name,
                 f"regions_pLDDT_greaterequal_{plddt_cutoff}" : \
                    ";".join(plddt_greaterequal_ranges),
                 f"regions_pLDDT_lower_{plddt_cutoff}" : \
                    ";".join(plddt_lower_ranges)})


        #------------------ Per-residue data frame -------------------#


        # Create a data frame merging per-residue information about
        # the secodnary structure assigned by DSSP and the pLDDT score
        merged_df = pd.merge(plddt_df, df, on = ["chain", "resnum"])        
        merged_df = merged_df[["chain","resnum","resname","pLDDT","secstruc"]]

        # Save the new data frame to a CSV file
        merged_csv = \
            os.path.join(protein_dir, f"{uniprot_id}.csv")
        merged_df.to_csv(merged_csv, sep = ",", index = False)


    #------------------- Protein regions data frame ------------------#


    # Generate a data frame containing the data about
    # the protein regions with pLDDT scores greater than/equal to
    # or lower than the pre-definied pLDDT cut-off for all proteins
    regions_df = pd.DataFrame(data)

    # Save the data frame in a CSV file
    regions_csv_file = \
        os.path.join(wd, "regions_pLDDT.csv")
    regions_df.to_csv(regions_csv_file, sep = ",", index = False)



if __name__ == "__main__":


    import argparse


    # Create the argument parser
    parser = argparse.ArgumentParser()

    # Add the arguments
    c_helpstr = "Configuration file for the curation."
    parser.add_argument("-c", "--config-file",
                        type = str,
                        required = True,
                        help = c_helpstr)

    d_helpstr = "Path to the working directory."
    parser.add_argument("-d", "--work-dir",
                        type = str,
                        default = os.getcwd(),
                        help = d_helpstr)

    # Parse the arguments
    args = parser.parse_args()

    # Get the working directory
    wd = args.work_dir

    # Load the configuration file
    config = yaml.safe_load(open(args.config_file, "r"))

    # Get and save AlphaFold data
    main(config = config,
         wd = wd)
