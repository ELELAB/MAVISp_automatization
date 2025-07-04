#!/usr/bin/env python

# Copyright (C) 2023 Ludovica Beltrame <beltrameludo@gmail.com>
# Danish Cancer Society

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import os
import sys
import time
import argparse
import pandas as pd
import logging as log
from Bio.PDB import *
from Bio.SeqUtils import seq3

def get_day(metatable):
    '''The function returns the day of the cancermuts
    run in the DD-MM-YYYY format'''

    # Time creation in seconds
    ti_c = os.path.getmtime(metatable)

    # Convert to standard time
    s_ti = time.ctime(ti_c)

    # Create a time object/structure
    t_obj = time.strptime(s_ti)
    
    return time.strftime("%d%m%Y", t_obj)

def keep_curated_mutations(row_value):
        # Split the string by commas and strip whitespace
        sources = [s.strip() for s in row_value.split(',')]
        # Keep the row if at least one substring does not contain 'saturation'
        return any(source and 'saturation' not in source for source in sources)

def pdb_dict(pdb_file):
    '''The functions reads the pdb file and returns
    a dictionary {resnum : resname}.'''

    # Parse target structure
    parser = PDBParser(QUIET=True)
    structure_id = 'target'
    filename = pdb_file
    structure = parser.get_structure(structure_id, filename)

    # Make dictionary {residue_position : residue_name}
    pdb_dict = {}
    for model in structure.get_list():
        for chain in model.get_list():
            for residue in chain.get_list():
                resseq = residue.get_id()[1]
                resname = residue.get_resname()
                pdb_dict.update({resseq: resname})
    
    return pdb_dict

def check_WT(metatable, pdb_dict, error_file,exclude_sanity_check):
    '''The function compares the wt residue of the metatable
    and of the pdb and returns a filtered metatable not including
    the rows with unmatched wt residues. Furthermore, it creates
    an error file which will be including the wrong wt.'''

    # Map in the metatable the residue name based on position
    # and store the information in a new column
    metatable['wt_pdb'] = metatable['aa_position'].map(pdb_dict)
    
    # Check pdb range and cancermuts metatable range
    # to be correspondent
    if metatable['wt_pdb'].isnull().values.any():
        if not exclude_sanity_check:
            log.error('Please check the residue range in the pdb(s) and ' \
                       'the range used in cancermuts correspond. \n' \
                       'Inconsistencies found. Exiting...')
            sys.exit(1)
        else:
            log.error('The residue range in the pdb(s) and ' \
                      'the range used in cancermuts does not'\
                      'correspond. This is probably due to the' \
                      'presence of missing residues in the pdb file \n' \
                      'The sanity check on the pdb file has been removed,'\
                      'so the mutation lists will be generated anyway.')


    # Create dataframe containing the unmatched residues
    error_df = metatable.loc[~(metatable['wt_3L'].str.upper() == metatable['wt_pdb'])].drop_duplicates(subset = ['wt_3L', 'aa_position', 'wt_pdb']).copy()
    
    # If the error dataframe contains unmatches
    if not exclude_sanity_check and not error_df.empty:

        # Print warning
        number = error_df.shape[0]
        
        # Rename columns
        dict_names = {'wt_3L': 'wt_metatable', 'aa_position': 'resnum'}
        error_df.rename(columns = dict_names, inplace=True)

        # Capitalize first letter of wt_pdb
        error_df['wt_pdb'] = error_df['wt_pdb'].apply(lambda x: x.capitalize())

        # Save error file (cancermuts_wt, position, pdb_wt)
        error_df.to_csv(error_file, index = False, columns = ['wt_metatable', 'resnum', 'wt_pdb'])

        log.error(f'{number} unmatched WT residue(s) found. ' \
                  f'See {error_file} for details. Exiting...')
        sys.exit(1)
    else:
        # Print info that no unmatched residues have been found
        if not exclude_sanity_check:
            log.info('All the WT residues match correctly.')

    checked_df = metatable.loc[(metatable['wt_3L'].str.upper() == metatable['wt_pdb'])].copy()
    checkef_df = checked_df.drop(columns = 'wt_pdb')

    return checked_df

def mutatex_format(dataframe, chain, day):
    '''Mutation format: wt chain position mut (e.g., AA75C)'''

    #Add mutation column
    dataframe['mutatex_mutation'] = dataframe['ref_aa'] + \
                                    chain + \
                                    dataframe["aa_position"].astype(str) + \
                                    dataframe["alt_aa"]
    
    return dataframe['mutatex_mutation'].to_csv(f'mutlist_mutatex_{day}.txt', header=None, index=None)

def mutatexP_format(dataframe, chain, day):
    '''Mutation format: wt chain position mut (e.g. TA75C)
    and wt chain position phospho_type (e.g., TA75p)'''

    #Define dictionary of phosphorylations
    phospho_dic = {'T': 'p', 'S': 's', 'Y': 'y'}
    
    #Filter dataframe for positions supporting phosphorylation
    phospho_dataframe = dataframe[dataframe['phosphorylation_site'] == 'P'].copy()
    
    #Add mutations column
    phospho_dataframe['mutatex_mutation'] = phospho_dataframe['ref_aa'] + \
                                            chain + \
                                            phospho_dataframe["aa_position"].astype(str) + \
                                            phospho_dataframe["alt_aa"]
    
    phospho_dataframe['phospho_type'] = phospho_dataframe['ref_aa'].apply(lambda x: phospho_dic.get(x))
    
    phospho_dataframe['mutatex_P_mutation'] = phospho_dataframe['ref_aa'] + \
                                              chain + \
                                              phospho_dataframe["aa_position"].astype(str) + \
                                              phospho_dataframe['phospho_type']
    
    return pd.concat([phospho_dataframe["mutatex_mutation"], phospho_dataframe["mutatex_P_mutation"].drop_duplicates()], \
                axis = 0).to_csv(f'mutlist_mutatex_P_{day}.txt', header=None, index=None)

def rosetta_format(dataframe, chain, domain, day, multi):
    '''Mutation format: chain.wt.position.mut  chain (e.g., A.A.75.C A).
    The function returns both the full length mutlist and the '''

    #Add mutation column
    dataframe['rosetta_mutation'] = chain + '.' + \
                                    dataframe['ref_aa'] + '.' + \
                                    dataframe["aa_position"].astype(str) + '.' + \
                                    dataframe["alt_aa"] + ' ' + \
                                    chain
    
    # Filter dataframe if multiple domains provided
    if multi == True:
        subset_df, n, c = subset_dataframe(dataframe, domain)

    return dataframe['rosetta_mutation'].to_csv(f'mutlist_rosetta_{day}.txt', header=None, index=None), \
            subset_df['rosetta_mutation'].to_csv(f'mutlist_rosetta_{n}-{c}_{day}.txt', header=None, index=None) \
                if multi == True else dataframe['rosetta_mutation'].to_csv(f'mutlist_rosetta_{day}.txt', header=None, index=None)

def hgvs_format(dataframe, day):
    '''Mutation format: p. wt3L position mut3L (e.g., p.Ala75Cys)'''

    #Add mutation column
    dataframe['mutation_hgvs'] = 'p.' + \
                                 dataframe['wt_3L'] + \
                                 dataframe["aa_position"].astype(str) + \
                                 dataframe["mut_3L"]

    return dataframe['mutation_hgvs'].to_csv(f'mutlist_hgvs_{day}.txt', header=None, index=None)

def cabsflex_format(dataframe, pdb, chain, domain, day):
    '''Mutation format: pdb chain position wt mut (e.g., 2XWR.pdb A 75 ALA CYS)'''
    
    subset_df, n, c = subset_dataframe(dataframe, domain)
    
    # Add mutation column
    subset_df['mutation_cabsflex'] = pdb + ' ' + \
                                     chain + ' ' + \
                                     subset_df['aa_position'].astype(str) + ' ' + \
                                     subset_df['wt_3L'].str.upper() + ' ' + \
                                     subset_df['mut_3L'].str.upper()

    return subset_df['mutation_cabsflex'].to_csv(f'mutlist_cabsflex_{n}-{c}_{day}.txt', header=None, index=None)

def elm_mutlist(dataframe, day):
    '''Mutation formate: one letter. The output is a dataframe
    containing the mutation in the first column and the linear 
    motif annotation in the second one.'''

    # Filter dataframe to remove empty linear motif rows
    dataframe = dataframe[pd.notnull(dataframe['linear_motif'])].copy()

    #Add mutation column
    dataframe['mutation'] = dataframe['ref_aa'] + \
                            dataframe["aa_position"].astype(str) + \
                            dataframe["alt_aa"]
    
    return dataframe[['mutation', 'linear_motif']].to_csv(f'mutlist_ELM_{day}.txt', index=None)


def subset_dataframe(dataframe, domain):
    '''The function filters the metatable to keep only
    the mutations within a certain range of residues.'''

    # Derive Nter and Cter
    nter, cter = domain.split(':')

    # Derive a subset of the dataframe only
    # containing the residues of selected domain
    subset_df = dataframe[dataframe.aa_position.between(left=int(nter), 
                                                        right=int(cter), 
                                                        inclusive='both')].copy()

    return subset_df, nter, cter

def main():   

    # Basic logging configuration
    log.basicConfig(level=log.INFO, \
                    format='%(levelname)s - %(message)s')

    # Add arguments required to the script to a argparse.ArgumentParser instance.
    description = "Get different mutation lists from cancermuts metatable"
    parser = argparse.ArgumentParser(description = description)

    m_helpstr = "Input cancermuts metatable"
    parser.add_argument("-m", "--metatable",
                        action = "store",
                        type = str,
                        help = m_helpstr,
                        required = True)

    d_helpstr = "Range of each domain.\n If multiple domains are given, " \
                "please provide them space-separated (e.g., 1:150 300-750)."
    parser.add_argument("-d", "--domain",
                        nargs = '+',
                        type = str,
                        help = d_helpstr,
                        required = True)

    ch_default = 'A'
    ch_helpstr = f"Chain of interest. (deafult = {ch_default})"
    parser.add_argument("-ch", "--chain",
                        type = str,
                        default = ch_default,
                        help = ch_helpstr)

    p_helpstr = "pdb(s) of interest.\n The number of provided pdbs " \
                "must match the number of domain ranges."
    parser.add_argument("-p", "--pdbfile",
                        nargs = '+',
                        type = str,
                        required = True,
                        help = p_helpstr)

    s_helpstr = "Remove the sanity check on the consistency of"\
                " positions and residues between the provided"\
                " PDB file and the input metatable."

    parser.add_argument("-s", "--exclude_sanity_check",
                        action = "store_true",
                        default = False,
                        required = False,
                        help = s_helpstr)

    mutatex_helpstr = "Generate mutatex mutlist"
    parser.add_argument("-M", "--mutatex",
                        action = 'store_true',
                        help = mutatex_helpstr)

    rosetta_helpstr = "Generate rosetta mutlist"
    parser.add_argument("-R", "--rosetta",
                        action = 'store_true',
                        help = rosetta_helpstr)

    hgvs_helpstr = "Generate HGVS mutlist"
    parser.add_argument("-H", "--hgvs",
                        action = 'store_true',
                        help = hgvs_helpstr)

    cabsflex_helpstr = "Generate cabsflex mutlist"
    parser.add_argument("-C", "--cabsflex",
                        action = 'store_true',
                        help = cabsflex_helpstr)

    exclude_saturation_helpstr = "Generate a filtered"\
                      " mutation list excluding "\
                      " those only present in the"\
                      " saturation mutagenesis list"

    parser.add_argument("-S", "--exclude_saturation",
                        action = 'store_true',
                        help = exclude_saturation_helpstr)

    date_helpstr = "Date in DDMMYYYY added to the output file names " \
                   "(default: day of the cancermuts run)"
    parser.add_argument("-D", "--day",
                        type = str,
                        help = date_helpstr)

    args = parser.parse_args()

    # Make dictionary from pdb(s)
    pdb_dictionary = {}
    for pdb in args.pdbfile:
        pdb_dictionary.update(pdb_dict(pdb))

    # Read cancermuts metatable
    df = pd.read_csv(args.metatable, index_col=0)
    
    # If a date is not provided by the user,
    # derive it from the creation date of the metatable
    if args.day is None:
        date = get_day(args.metatable)
    else:
        date = args.day

    # Filter out lines without source
    df = df[pd.notnull(df['sources'])]

    # Filter out lines in which the mutant
    # residue is not a standard amino acid
    aa1 = list("ACDEFGHIKLMNPQRSTVWY")
    df = df[df['alt_aa'].isin(aa1)]

    # Filter out lines with aa_position not covered by the structure
    filtered_df = pd.DataFrame()
    for r in args.domain:
        n, c = r.split(':')
        filtered_df = filtered_df.append(df[df.aa_position.between(left=int(n), 
                                                                   right=int(c), 
                                                                   inclusive='both')])

    # Convert letter code of wt and mutant residues
    filtered_df['wt_3L'] = filtered_df['ref_aa'].apply(seq3)
    filtered_df['mut_3L'] = filtered_df['alt_aa'].apply(seq3) 
    
    # Check WT residues
    input_df = check_WT(filtered_df, pdb_dictionary, 'wt.err',args.exclude_sanity_check)
    
    # Generate one letter mutlist 
    output = f'mutlist_{date}.txt'

    input_df['mutation'] = input_df['ref_aa'] + \
                           input_df["aa_position"].astype(str) + \
                           input_df["alt_aa"]
    input_df['mutation'].to_csv(output, header=None, index=None)

    if args.exclude_saturation:
        output = f'curated_mutlist_{date}.txt'
        curated_input_df = input_df[input_df['sources'].apply(keep_curated_mutations)]
        curated_input_df['mutation'].to_csv(output, header=None, index=None)

    # Generate mutatex mutlist
    if args.mutatex:
        mutatex_format(input_df, args.chain, date)
        mutatexP_format(input_df, args.chain, date)

    # Generate HGV mutlist
    if args.hgvs:
        hgvs_format(input_df, date)

    # Generate rosetta mutlist
    if args.rosetta:
        if len(args.domain) > 1:
            for domain in args.domain:
                rosetta_format(input_df, args.chain, domain, date, multi=True)
        else:
            rosetta_format(input_df, args.chain, args.domain[0], date, multi=False)

    
    # Generate cabsflex mutlist
    if args.cabsflex:
        if not len(args.pdbfile) == len(args.domain):
            log.error('A pdb ID for each domain must be provided ' \
                      'to get the CABSFLEX mutlist. Exiting...')
            sys.exit(1)
        else:
            for i in range(0, len(args.pdbfile), 1):
                pdb, domain = args.pdbfile[i], args.domain[i]
                cabsflex_format(input_df, pdb, args.chain, domain, date)

    # Generate ELM mutlist
    elm_mutlist(input_df, date)
    
    log.info('Mutlists done!')

main()
