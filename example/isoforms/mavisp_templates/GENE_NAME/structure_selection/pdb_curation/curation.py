#!/usr/bin/env python
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-



# Standard library
import copy
import logging
import os
import re
import subprocess
from tempfile import NamedTemporaryFile
# Third-party packages
import Bio
from Bio import AlignIO
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB import PDBIO, PDBParser, Select
import pandas as pd
import requests
import yaml
import zipfile
import pathlib

# Set up the logging
logging.basicConfig(level = logging.INFO)

# Get the logger
log = logging.getLogger(__name__)



########################## HELPER FUNCTIONS ###########################



def get_uniprot_fasta(uniprot_id,
                      fasta_file):
    """Get the FASTA sequences associated to the protein with the
    given UniProt ID and write it to an output file.
    """

    # Generic URL for FASTA sequences in UniProt
    UNIPROT_URL = "https://www.uniprot.org/uniprot/{:s}.fasta"

    # Get the data
    data = requests.get(UNIPROT_URL.format(uniprot_id)).text

    # Write the data to the FASTA file
    with open(fasta_file, "w") as out:
        out.write(data)


def get_pdb_structure(pdb_id,
                      pdb_file):
    """Get the structure (in PDB format) associated with a
    given PDB ID and write it to an output file.
    """

    # Generic URL for PDB files in the PDB
    PDB_URL = "https://files.rcsb.org/view/{:s}.pdb"

    # Get the data
    data = requests.get(PDB_URL.format(pdb_id)).text

    # Write the data to the PDB file
    with open(pdb_file, "w") as out:
        out.write(data)


def get_mmcif_structure(pdb_id,
                        mmcif_file):
    """Get the structure (in mmCIF format) associated with a
    given PDB ID and write it to an output file.
    """

    # Generic URL for mmCIF files in the PDB
    MMCIF_URL = "https://files.rcsb.org/view/{:s}.cif"

    # Get the data
    data = requests.get(MMCIF_URL.format(pdb_id)).text

    # Write the data to the mmCIF file
    with open(mmcif_file, "w") as out:
        out.write(data)


def get_pdbredo_files(pdb_id,
                      pdbredo_dir):
    """Get the PDBredo files associated with a given PDB ID
    and store them in an output directory.
    """

    # Generic URL for PDBredo ZIP files
    PDBREDO_URL = "https://pdb-redo.eu/db/{:s}/zipped"

    # Get the data (make sure the PDB ID is given in
    # all lowercase letters)
    data = \
        requests.get(PDBREDO_URL.format(pdb_id.lower())).content

    # Create a temporary file (it is deleted as soon
    # as the 'with' block ends)
    with NamedTemporaryFile("wb") as out:

        # Write the data to a temporary ZIP file
        out.write(data)

        # Unzip the file to a directory
        with zipfile.ZipFile(out.name) as zip_content:
            for zip_info in zip_content.infolist():
                path = pathlib.Path(zip_info.filename)
                zip_info.filename = str(pathlib.Path(*path.parts[1:]))
                zip_content.extract(zip_info, path=pdbredo_dir)

def parse_clustalo_output(out_clustalo):
    """Get the sequence coverage of a PDB chain with respect
    to the UniProt sequence of the protein of interest and
    whether there are mutations in the PDB chain.

    The input alignment (output from Clustal Omega) is supposed
    to be an alignment in Clustal format containing two sequences
    aligned, the first one being the UniProt sequence of the protein
    of interest and second one being the sequence of a PDB chain
    representing (a) region(s) of the protein of interest.
    """

    # Read the alignment in Clustal format
    al = AlignIO.read(out_clustalo, format = "clustal")
    
    # Create a list to store the starting and ending points of
    # the protein UniProt sequence covered by the PDB sequence
    regions_covered = []

    # Create a list to store the mutations found in the PDB
    # sequence with respect to the UniProt sequence
    mutations = []

    # Initialize the residue number and the starting residue
    # number to zero
    res_num = 0
    res_start = 0

    # Inform the user about the UniProt and PDB IDs of
    # the sequences present in the alignment
    log.info(f"UniProt sequence ID: {al[0].id}")
    log.info(f"PDB sequence ID: {al[1].id}")

    # Initialize the variable storing the previous residue
    # in the PDB sequence to an empty string
    prev_res_pdb = ""
    
    # For each residue in the UniProt sequence and aligned
    # residue in the PDB sequence
    for res_uniprot, res_pdb in zip(al[0].seq, al[1].seq):

        # If the residue in the UniProt sequence is not
        # a gap residue (there should be no gaps in the
        # UniProt sequence)
        if res_uniprot != "-":

            # Update the counter for the number of the
            # residue under consideration
            res_num += 1

            # If the UniProt residue and the PDB residue
            # coincide
            if res_uniprot == res_pdb:

                # If the previous residue was a gap and
                # this is the first region covered
                # by the PDB sequence
                if prev_res_pdb == "-" or res_start == 0:

                    # The number of the starting residue of
                    # this region will be the current
                    # residue number
                    res_start = res_num

            # If the residue found in the PDB sequence differs
            # from the residue found in the UniProt sequence
            else:
                
                # If the residue is a gap
                if res_pdb == "-":

                    # If the previous residue was also a gap, keep
                    # on going (the region covered before the gap
                    # region has been already saved when the
                    # first gap symbol was encountered)
                    if prev_res_pdb == "-":
                        continue
                    
                    # If this is a gap region in the middle of
                    # the protein or at the end (but not at
                    # the beginning)
                    if res_start != 0:
                        
                        # Update the regions covered by the PDB
                        # sequence with the region covered so far
                        regions_covered.append((res_start, res_num-1))
                
                # If the residue is not a gap
                else:

                    # Set the start of the first region covered
                    # by the PDB sequence to the current residue
                    # number
                    res_start = res_num
                    
                    # Update the list of mutations
                    mutations.append(\
                        f"{res_uniprot}.{res_num}.{res_pdb}")

        # Update the variable storing the previous residue in the
        # PDB sequence
        prev_res_pdb = res_pdb

    # If at least one region was found (if the PDB sequence
    # covers the entire UniProt sequence or only one region
    # of the sequence)
    if regions_covered and res_start != regions_covered[-1][0]:
        
        # Update the regions covered by the PDB
        # sequence with the region covered so far
        regions_covered.append((res_start, res_num))

    # Convert the regions covered into strings for easier
    # output formatting
    regions_covered = \
        [f"{start}-{end}" for start, end in regions_covered]

    # Inform the user about the regions covered and
    # mutations found
    log.info(\
        f"Regions of the Uniprot sequence covered by the " \
        f"PDB sequence: {', '.join(regions_covered)}.")

    if mutations:
        log.info(\
            f"Mutations in the PDB sequence with respect to " \
            f"the UniProt sequence: {', '.join(mutations)}.")
    else:
        log.info("No mutations found.")

    # Return regions covered and mutations
    return regions_covered, mutations


def clean_pdb(res_to_remove,
              pdb_file,
              pdb_file_cleaned):
    """Remove unwanted residues/cofactors/water/etc. from
    a PDB file.
    """

    # Create a list to store the residues removed and a set
    # to store the IDs of the residues to remove
    removed_res = []
    removed_res_ids = set()
    
    # Create the PDB parser
    parser = PDBParser()

    # Get the structure 
    structure = parser.get_structure("struct", pdb_file)

    # For each model in the structure
    for model in structure:

        # For each chain in the model
        for chain in model:

            # For each residue in the chain
            for residue in chain:

                # Get the residue name
                res_name = residue.get_resname()

                # If the residue name is in the list
                # of residues to remove
                if res_name in res_to_remove:

                    # Store the chain ID, the residue number,
                    # and the residue name in the list
                    # of removed residues
                    s, m, chain_id, res_id = residue.get_full_id()
                    removed_res.append((chain_id, res_id[1], res_name))

                    # Store the full residue ID to remove it from
                    # the chain later
                    removed_res_ids.add(residue.id)

    # Create a class that defines a residue selector (only
    # residues respecting the criteria set will be included
    # in the final PDB)
    class ResSelect(Select):
        
        def accept_residue(self, res):
            # Ignore the residue if the ID is in the list of
            # IDs of residues to remove
            if res.id in removed_res_ids:
                return False
            else:
                return True


    # Create the PDB I/O object
    io = PDBIO()

    # Set the structure
    io.set_structure(structure)

    # Save the cleaned PDB file using the residue selector
    io.save(pdb_file_cleaned, ResSelect())

    # Return the list of removed residues
    return removed_res


def extract_first_model(pdb_file,
                        pdb_file_extracted):
    """Extract only the first model in a NMR PDB structure.
    """

    # Helper class to select only the model of interest
    class ModelSelect(Select):

        def __init__(self, model_number):
            self.model_number = model_number
        
        def accept_model(self, model):
            
            # Discard all models not having the correct model
            # number.
            # Using serial numbers instead of ids because
            # models are numbered from 1 but ids start from 0.
            if model.serial_num != self.model_number:
                return 0
            
            return 1

    # Create a PDB parser
    parser = PDBParser()

    # Parse the structure
    name = pdb_file.replace(".pdb", "")
    structure = parser.get_structure(name, pdb_file)

    # Save the processed structure
    w = PDBIO()
    w.set_structure(structure)
    w.save(pdb_file_extracted, ModelSelect(model_number = 1))


def run_pdb_fromcif(pdb_fromcif_exec,
                    mmcif_file,
                    pdb_file):
    """Run the 'pdb_fromcif' tool from pdb-tools to
    convert a file in mmCIF format to a file in
    PDB format.
    """

    # Convert the file
    pdb_fromcif_proc = \
        subprocess.run([pdb_fromcif_exec, mmcif_file],
                       capture_output = True,
                       encoding = "utf-8")

    # Save the standard output to the output PDB file
    with open(pdb_file, "w") as out:
        out.write(pdb_fromcif_proc.stdout)


def run_pdb_splitchain(pdb_splitchain_exec,
                       pdb_file,
                       pdbs_dir):
    """Run the 'pdb_splitchain' tool from pdb-tools to
    split a PDB file into as many PDB files as its
    constituent chains.
    """

    # Get the current working directory
    curr_wd = os.getcwd()

    # Go to the directory where the output PDB files should
    # be saved
    os.chdir(pdbs_dir)

    # Run 'pdb_splitchain' to split the PDB file into multiple
    # PDB files, each one containing one of the constituent
    # chains
    pdb_splitchain_proc = \
        subprocess.run([pdb_splitchain_exec,
                        pdb_file])

    # Go back to the working directory
    os.chdir(curr_wd)


def run_pdb2fasta(pdb2fasta_exec,
                  pdb_file,
                  fasta_file):
    """Run 'pdb2fasta' to extract the FASTA sequence from the PDB.
    """

    # Extract the FASTA sequence from the PDB
    pdb2fasta_proc = \
        subprocess.run([pdb2fasta_exec, pdb_file],
                       capture_output = True,
                       encoding = "utf-8")

    # Write out the captured output (FASTA sequence) to a new
    # FASTA file
    with open(fasta_file, "w") as out:
        out.write(pdb2fasta_proc.stdout)


def run_pdb4amber(pdb4amber_exec,
                  pdb_file,
                  out_dir):
    """Run 'pdb4amber' (AMBER tools) on a given PDB file.
    """

    # Get the name of the PDB file (without extension)
    pdb_name = os.path.splitext(os.path.basename(pdb_file))[0]

    # Set the path for the output PDB file
    out_pdb = os.path.join(out_dir, f"{pdb_name}_amber.pdb")

    # Set the path for the output log file
    out_log = os.path.join(out_dir, f"{pdb_name}_amber.log")

    # Run
    pdb4amber_proc = \
        subprocess.run([pdb4amber_exec,
                        "-i", pdb_file,
                        "-o", out_pdb,
                        "-l", out_log])


def run_cat(in_files,
            out_file):
    """Merge multiple input files into a single output
    file with 'cat'.
    """

    # Merge the input files with 'cat'
    cat_proc = \
        subprocess.run(["cat", *in_files],
                        capture_output = True,
                        encoding = "utf-8")

    # Save the captured output to the output file
    with open(out_file, "w") as out:
        out.write(cat_proc.stdout)


def run_clustalo(clustalo_exec,
                 in_align,
                 out_align):
    """Run a sequence alignment with Clustal Omega.
    """

    # Run Clustal Omega
    clustalo_proc = \
        subprocess.run([clustalo_exec,
                        "-i", in_align,
                        "-o", out_align,
                        "--outfmt=clustal",
                        "--resno"])


def main(config,
         wd,
         curation_dirs,
         seqcov_mut_csv):

    

    # Inform the user about the versions of the third-party packages
    # used
    log.info(f"Using 'biopython' version {Bio.__version__}.")
    log.info(f"Using 'pandas' version {pd.__version__}.")
    log.info(f"Using 'requests' version {requests.__version__}.")
    log.info(f"Using 'pyyaml' version {yaml.__version__}.")


    #---------------------- UniProt FASTA files ----------------------#


    # Create an empty dictionary to map each UniProt IDs to the path
    # of the file where the FASTA sequence associated to the protein
    # represented by the ID is stored
    uniprot_fasta_files = {}

    # Get the FASTA files with the UniProt sequences of the proteins
    # of interest
    for uniprot_id in config["uniprot_ids"]:

        # Set the path to the FASTA file where data from UniProt
        # will be written
        uniprot_fasta_file = \
            os.path.join(curation_dirs["uniprot"],
                         f"{uniprot_id}.fasta")

        # Get the data and write the FASTA file
        get_uniprot_fasta(uniprot_id, uniprot_fasta_file)

        # Save the path to the FASTA file to the dictionary
        uniprot_fasta_files[uniprot_id] = uniprot_fasta_file

    # Create an empty list to store information about the
    # sequence converage of the single PDB chains with respect
    # to the UniProt sequences of the proteins they represent
    seqcov_mut_list = []


    #------------------------ Structure files ------------------------#


    # For each PDB ID and the data associated to it
    for pdb_id, pdb_data in config["pdb_ids"].items():

        # Get the data associated to the PDB ID
        file_types, exp_method, chain2uniprot, res_to_remove = \
            pdb_data.get("file_types", []), \
            pdb_data.get("exp_method", ""), \
            pdb_data.get("chain2uniprot", ""), \
            pdb_data.get("res_to_remove", [])


        #--------------------- Get the structure ---------------------#


        # If the structure in PDB format is requested
        if "pdb" in file_types:

            # Set the path to the file that will contain the structure
            pdb_file = \
                os.path.join(curation_dirs["pdbs"], f"{pdb_id}.pdb")

            # Get the PDB structure
            get_pdb_structure(pdb_id, pdb_file)

        # If the structure in mmCIF format is requested
        if "mmcif" in file_types:

            # Set the path to the file that will contain the
            # mmCIF structure
            mmcif_file = \
                os.path.join(curation_dirs["pdbs"], f"{pdb_id}.cif")

            # Get the mmCIF structure
            get_mmcif_structure(pdb_id, mmcif_file)

            # If the structure was not requested in PDB format
            # (possibly because the PDB is not available in the
            # Protein Data Bank)
            if "pdb" not in file_types:

                # Set the path to the file that will contain the
                # PDB file generated from the mmCIF file
                pdb_file = \
                    os.path.join(curation_dirs["pdbs"], f"{pdb_id}.pdb")

                # Convert the mmCIF file to a PDB file
                run_pdb_fromcif(config["pdb_fromcif_exec"],
                                mmcif_file,
                                pdb_file)


        #---- Check for alt states, missing atoms/residues, etc. -----#


        # Create a directory for each PDB inside the one dedicated
        # to PDB4amber results
        pdb4amber_dir = \
            os.path.join(curation_dirs["p4a"], pdb_id)
        os.makedirs(pdb4amber_dir, exist_ok = True)

        # Run pdb4amber to check for alternative states, missing atoms,
        # missing residues, etc.
        run_pdb4amber(config["pdb4amber_exec"],
                      pdb_file,
                      pdb4amber_dir)


        #------------------------ NMR model 1 ------------------------#

        
        # If the structure is a NMR structure, take only the first
        # model
        if exp_method == "nmr":

            # Set the path to the processed PDB file
            pdb_file_extracted = \
                os.path.join(curation_dirs["pdbs_proc"],
                             f"{pdb_id}_model1.pdb")

            # Extract only the first model
            extract_first_model(pdb_file,
                                pdb_file_extracted)

            # The structure extracted will be the one processed
            # from now on
            pdb_file = pdb_file_extracted


        #-------------------- Clean the structure --------------------#


        # Set the path to the cleaned PDB file
        pdb_file_cleaned = \
            os.path.join(curation_dirs["pdbs_proc"],
                         f"{pdb_id}.pdb")

        # Clean the PDB file
        removed_res = \
            clean_pdb(res_to_remove,
                      pdb_file,
                      pdb_file_cleaned)

        # If residues have been removed
        if removed_res:

            # Inform the user about the residues removed
            removed_res_str = \
                ", ".join(["{:s}-{:d}{:s}".format(c, rnum, rname) \
                           for c, rnum, rname in removed_res])
            logstr = \
                f"The following residues have been removed from " \
                f"{pdb_file_cleaned}: {removed_res_str}"
            log.info(logstr)


        #--------------- Get FASTA from the structure ----------------#


        # Set the path to the FASTA file that will contain
        # the sequence extracted from the PDB structure
        pdb_fasta_file = \
            os.path.join(curation_dirs["p2f"],
                         f"{pdb_id}.fasta")

        # Run pdb2fasta and get the FASTA sequence from
        # the PDB structure
        run_pdb2fasta(config["pdb2fasta_exec"],
                      pdb_file_cleaned,
                      pdb_fasta_file)


        #-------------------- Split the structure --------------------#


        # Split the PDB file into its constituent chains
        run_pdb_splitchain(config["pdb_splitchain_exec"],
                           pdb_file_cleaned,
                           curation_dirs["pdbs_dec"])

        # For the PDB files with the splitted chains
        for _pdb_file in os.listdir(curation_dirs["pdbs_dec"]):


            #---------- Get the structure's FASTA sequence -----------#


            # Consider only the PDBs of interest
            if _pdb_file.startswith(pdb_id):

                # Get the PDB name
                _pdb_name = \
                    os.path.splitext(os.path.basename(_pdb_file))[0]

                # Set the path to the FASTA file that will contain
                # the sequence extracted from the PDB structure
                _pdb_fasta_file = \
                    os.path.join(curation_dirs["p2f_dec"],
                                 f"{_pdb_name}.fasta")

                # Get the full path to the PDB file
                _pdb_file = \
                    os.path.join(curation_dirs["pdbs_dec"], _pdb_file)

                # Run pdb2fasta and get the FASTA sequence from
                # the PDB structure
                run_pdb2fasta(config["pdb2fasta_exec"],
                              _pdb_file,
                              _pdb_fasta_file)

                # Get the chain contained in the current PDB file
                chain = _pdb_name.split("_")[1]

                # Try to get the UniProt ID of the protein covered
                # by that chain
                try:
                    uniprot_id_chain = chain2uniprot[chain]
                
                # If the chain was not found in the configuration,
                # do not consider it
                except KeyError:

                    # Warn the user and go to the next structure
                    logstr = \
                        f"Chain '{chain}' of '{pdb_id}' is not " \
                        f"specified in the configuration file. " \
                        f"It will be ignored for alignment " \
                        f"purposes."
                    log.warning(logstr)
                    continue

                # Get the UniProt sequence the chain should be
                # aligned to
                uniprot_fasta_to_align = \
                    uniprot_fasta_files[uniprot_id_chain]

                # Set the path to the file containing the input
                # sequences
                in_align = \
                    os.path.join(\
                        curation_dirs["align"],
                        f"{uniprot_id_chain}_and_{_pdb_name}.fasta")

                # Set the path to the file containing the output
                # alignment
                out_align = \
                    os.path.join(\
                        curation_dirs["align"],
                        f"{uniprot_id_chain}_and_{_pdb_name}_aligned.fasta")

                # Merge the input sequences into a single file
                run_cat([uniprot_fasta_to_align, _pdb_fasta_file], 
                         in_align)
                    
                # Run a pairwise alignment with Clustal Omega
                run_clustalo(config["clustalo_exec"],
                             in_align,
                             out_align)

                # Parse the output alignment file
                regions_covered, mutations = \
                    parse_clustalo_output(out_align)

                # Update the list summarizing the info about
                # the sequence converage of the PDB files
                seqcov_mut_list.append(\
                    {"uniprot_id" : uniprot_id_chain,
                     "pdb_id" : pdb_id,
                     "pdb_chain" : chain,
                     "sequence_coverage" : ";".join(regions_covered),
                     "mutations" : ";".join(mutations)})

        # If the structure comes from x-ray crystallography
        if exp_method == "xray":

            # Set the path to the directory where the PDBredo
            # files will be stored
            pdbredo_dir = \
                os.path.join(curation_dirs["pdbs_redo"], pdb_id)

            # Get the associated PDBredo files
            get_pdbredo_files(pdb_id, pdbredo_dir)

            # Get the PDBredo structure in PDB format
            pdb_redo_file = \
                os.path.join(pdbredo_dir, f"{pdb_id.lower()}_final.pdb")


            #-------------------- Clean the structure --------------------#


            # Set the path to the cleaned PDB file
            pdb_redo_file_cleaned = \
                os.path.join(curation_dirs["pdbs_redo_proc"],
                             f"{pdb_id.lower()}_final.pdb")

            # Clean the PDB file
            removed_res = \
                clean_pdb(res_to_remove,
                          pdb_redo_file,
                          pdb_redo_file_cleaned)

            # If residues have been removed
            if removed_res:

                # Inform the user about the residues removed
                removed_res_str = \
                    ", ".join(["{:s}-{:d}{:s}".format(c, rnum, rname) \
                               for c, rnum, rname in removed_res])
                logstr = \
                    f"The following residues have been removed from " \
                    f"{pdb_redo_file_cleaned}: {removed_res_str}"
                log.info(logstr)


            #------------------ Split the structure ------------------#


            # Split the PDB file into its constituent chains
            run_pdb_splitchain(config["pdb_splitchain_exec"],
                               pdb_redo_file_cleaned,
                               curation_dirs["pdbs_redo_dec"])


    #------------------- Sequence coverage CSV file ------------------#


    # Create a Pandas data frame storing the sequence coverage
    # and mutations associated to each PDB chain
    seqcov_mut_df = pd.DataFrame(seqcov_mut_list)

    # Save the data frame to the designated output CSV file
    seqcov_mut_df.to_csv(seqcov_mut_csv, sep = ",")



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

    # Names of the directories needed for the curation of the PDBs
    curation_dirs = \
        {"align" : os.path.join(wd, "alignment"), 
         "p2f" : os.path.join(wd, "pdb2fasta"), 
         "p2f_dec" : os.path.join(wd, "pdb2fasta_decomposed"),
         "p4a" : os.path.join(wd, "pdb4amber"),
         "pdbs" : os.path.join(wd, "pdbs"), 
         "pdbs_dec" : os.path.join(wd, "pdbs_decomposed"),
         "pdbs_redo" : os.path.join(wd, "pdbs_redo"),
         "pdbs_redo_proc" : os.path.join(wd, "pdbs_redo_processed"), 
         "pdbs_redo_dec" : os.path.join(wd, "pdbs_redo_decomposed"),
         "pdbs_proc" : os.path.join(wd, "pdbs_processed"),
         "uniprot" : os.path.join(wd, "uniprot")}

    # For each directory
    for curation_dir in curation_dirs.values():

        # Do not create it if it already exists
        os.makedirs(curation_dir, exist_ok = True)

    # Name of the output CSV file that will contain information
    # about the sequence coverage of each PDB chain
    seqcov_mut_csv = \
        os.path.join(curation_dirs["align"], "seqcov_mutations.csv")

    # Run the curation
    main(config = config,
         wd = wd,
         curation_dirs = curation_dirs,
         seqcov_mut_csv = seqcov_mut_csv)
