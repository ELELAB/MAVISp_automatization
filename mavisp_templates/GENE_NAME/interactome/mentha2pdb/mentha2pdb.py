# -*- coding: utf-8 -*-
"""
Created on Thu May 12 09:14:52 2022
updated December 2023 v 1.4

@author: Matteo Lambrughi
"""

import json
import os
from os.path import isfile, join
from pathlib import Path
import shutil
from concurrent.futures import ThreadPoolExecutor, as_completed
import sys
import argparse
import time
from decimal import Decimal
import numpy as np
import pandas as pd
import re
import requests
import warnings
import csv


def get_pdb_entries_for_uniprot(uniprot_id):
    """
    Queries PDB for entries based on UniProt Accession Code (AC) and human taxonomy ID (9606).
    Returns list of PDB IDs, or [] if none found.
    """
    url = "https://search.rcsb.org/rcsbsearch/v2/query"
    headers = {'Content-Type': 'application/json'}

    payload = {
        "query": {
            "type": "group",
            "logical_operator": "and",
            "nodes": [
                {
                    "type": "group",
                    "logical_operator": "and",
                    "nodes": [
                        {
                            "type": "terminal",
                            "service": "text",
                            "parameters": {
                                "attribute": "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession",
                                "operator": "in",
                                "negation": False,
                                "value": [uniprot_id]
                            }
                        },
                        {
                            "type": "terminal",
                            "service": "text",
                            "parameters": {
                                "attribute": "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_name",
                                "operator": "exact_match",
                                "value": "UniProt",
                                "negation": False
                            }
                        }
                    ]
                },
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "rcsb_entity_source_organism.taxonomy_lineage.id",
                        "operator": "exact_match",
                        "negation": False,
                        "value": "9606"
                    }
                }
            ]
        },
        "return_type": "entry",
        "request_options": {
            "return_all_hits": True
        }
    }

    try:
        response = requests.post(url, headers=headers, json=payload)
        response.raise_for_status()
        result_data = response.json()
        return [entry["identifier"] for entry in result_data.get("result_set", [])]
    except requests.exceptions.RequestException as e:
        print(f"Error querying PDB for UniProt ID {uniprot_id}: {e}")
        return []


def make_target_interactor_sequence_files(dataframe_out):

    # get target list so we cover -x option (splitted outs) and normal (with all the targets in the same dataframe
    target_list = list(dict.fromkeys(dataframe_out['target uniprot id'].tolist()))

    Path("inputs_afmulti").mkdir(parents=True, exist_ok=True)

    url = 'https://rest.uniprot.org/uniref/search?query=uniprot_id:'

    for target in target_list:
        print('>>Making folders/files for target {}                   '.format(target))

        # filter dataframe
        target_data = dataframe_out[(dataframe_out['target uniprot id'] == target)]
        interactor_uniprot_ids = target_data['interactor uniprot id'].to_list()
        interactor_genes = target_data['interactor uniprot gene'].to_list()

        # get first gene value -> same target = all the same
        target_uniprot_gene = target_data['target uniprot gene'].values[0]

        # get target sequence
        result = make_request(url, 'get', target)

        target_sequence = ''

        #CHOSE RIGHT RESULT : Homo sapiens (Human) in organism name and target in id
        for res in result['results']:
            if target in res['id'] and res['representativeMember']['organismName'] == 'Homo sapiens (Human)':
                target_sequence = res['representativeMember']['sequence']['value']
                break
            else:
                # print(f"result discarded cause {target} not in {res['id']} \n\t "
                #       f"or \n\t {res['representativeMember']['organismName']} is not Homo sapiens (Human)")
                pass

        # fix for uniprot genes of type U2AF1L5 {ECO:0000312|HGNC:HGNC:51830} -> error creating folder
        # covering no space case U2AF1L5{ECO:0000312|HGNC:HGNC:51830} and space case U2AF1L5 {ECO:0000312|HGNC:HGNC:51830}
        if ' ' in target_uniprot_gene:
            target_uniprot_gene = target_uniprot_gene.split(' ')[0].rstrip()
        if '{' in target_uniprot_gene:
            target_uniprot_gene = target_uniprot_gene.split('{')[0].rstrip()

        # make dir for target
        Path("inputs_afmulti/" + target_uniprot_gene).mkdir(parents=True, exist_ok=True)

        for interactor_id, interactor_gene in zip(interactor_uniprot_ids, interactor_genes):
            # for every interactor make request make dir and then build file
            result = make_request(url, 'get', interactor_id)
            interactor_sequence = ''
            if result['results'] != []:
                ##########CHOSE RIGHT RESULT : Homo sapiens (Human) in organism name and target in id
                interactor_sequence = ''
                for res in result['results']:
                    if interactor_id in res['id'] and res['representativeMember']['organismName'] == 'Homo sapiens (Human)':
                        interactor_sequence = res['representativeMember']['sequence']['value']
                        break
                    else:
                        # print(f"result discarded cause "
                        #       f"{interactor_id} not in {res['id']} \n\t "
                        #       f"or \n\t "
                        #       f"{res['representativeMember']['organismName']} is not Homo sapiens (Human)")
                        pass
            else:
                print('***INTERACTOR {} of target {} returned NO results, skipping folder/sequence creation'.format(
                    interactor_id, target))
                continue

            if ' ' in interactor_gene:
                interactor_gene = interactor_gene.split(' ')[0].rstrip()
            if '{' in interactor_gene:
                interactor_gene = interactor_gene.split('{')[0].rstrip()

            # make dir for interactor
            Path("inputs_afmulti/" + target_uniprot_gene + '/' + interactor_gene).mkdir(parents=True, exist_ok=True)

            print('>>Made folder {}                     '.format(
                "inputs_afmulti/" + target_uniprot_gene + '/' + interactor_gene), end='\r')

            with open("inputs_afmulti/" + target_uniprot_gene + '/' + interactor_gene + '/input.fasta',
                      'w+') as alpha_file:
                alpha_file.write('>' + target_uniprot_gene + '\n')
                alpha_file.write(target_sequence + '\n')
                alpha_file.write('>' + interactor_gene + '\n')
                alpha_file.write(interactor_sequence + '\n')

    return 0


def pmid_adder(data, dataframe_out):
    # if p option selected -> PMID search and add
    # we have
    # dataframe data -> mentha db
    # dataframe dataframeOut
    df_out = dataframe_out.copy(deep=True)
    pmid_list = []
    for index, row in dataframe_out.iterrows():
        target_protein = row["target uniprot id"]
        interactor_protein = row["interactor uniprot id"]

        data_direct = data[(data['Protein A'] == target_protein) & (data['Protein B'] == interactor_protein)]
        data_invers = data[(data['Protein A'] == interactor_protein) & (data['Protein B'] == target_protein)]
        direct_pmid = ''
        invers_pmid = ''

        direct_pmid = data_direct['PMID']
        invers_pmid = data_invers['PMID']

        pmid = []
        pmid += direct_pmid.to_list()
        pmid += invers_pmid.to_list()
        pmid = ' '.join(x for x in list(dict.fromkeys(pmid)))

        pmid_list.append(pmid)

    df_out['PMID'] = pmid_list
    return df_out


def get_experiment(pdb):
    """
    This function retrieves PDB > experiment

    :param pdb: String,
    :return: resolution: String
    """

    url = 'https://www.ebi.ac.uk/pdbe/api/pdb/entry/experiment/'
    data = make_request(url, "get", pdb)
    resolution = ''

    if data is None:
        resolution = 'none'
        # print("########")
        # print("EXPERIMENT CALL ERROR -> no data")
        # print("########")
    else:
        if 'resolution' in data[pdb.lower()][0].keys():
            resolution = Decimal(str(data[pdb.lower()][0]['resolution']))
        else:
            # print('resolution not in data')
            resolution = 'na'

    return resolution


def get_summary(pdb):
    """
    This function retrieves PDB > summary

    :param pdb: String,
    :return: fused: String ('yes'/'')
             dna: String
             ligands: String
             method: String
    """

    url = 'https://www.ebi.ac.uk/pdbe/api/pdb/entry/summary/'
    data = make_request(url, "get", pdb)

    dna = ''
    ligands = ''
    method = ''
    title = ''
    fused = ''

    # data None -> means that we have no data from the request
    # no data -> no title
    # no title -> fused = ''
    if data is None:
        title = 'NO TITLE'
        dna = 'none'
        ligands = 'none'
        method = 'none'
        fused = 'none'
        # print("########")
        # print("PDB ", pdb, " has no title :(")
        # print("########")
    else:
        # get INFOs
        title = data[pdb.lower()][0]['title']
        dna = data[pdb.lower()][0]['number_of_entities']['dna']
        ligands = data[pdb.lower()][0]['number_of_entities']['ligand']
        method = data[pdb.lower()][0]['experimental_method'][0]

        if 'fused' in title or 'fusion' in title:
            fused = 'yes'
        else:
            fused = 'na'

    return fused, dna, ligands, method


def get_mappings_data(pdb, targetProtein, interactorProtein):
    """
    This function will GET the mappings data from
    the PDBe API using the make_request() function

    :param pdb: String
    :return: targetChainIds: String
             targetStart: String
             targetEnd: String
             interactorChainIds: String
             interactorStart: String
             interactorEnd: String
             otherInteractors: String
    """

    base_url = "https://www.ebi.ac.uk/pdbe/"
    api_base = base_url + "api/"
    uniprot_mapping_url = api_base + 'mappings/uniprot/'
    # Check if the provided PDB id is valid
    # There is no point in making an API call
    # with bad PDB ids
    if not re.match("[0-9][A-Za-z][A-Za-z0-9]{2}", pdb):
        # print("Invalid PDB id")
        return 'none', 'none', 'none', 'none', 'none', 'none', 'none'

    # GET the mappings data
    mappings_data = make_request(uniprot_mapping_url, "get", pdb)

    targetChainIds = []
    interactorChainIds = []
    targetStart = []
    targetEnd = []
    interactorStart = []
    interactorEnd = []
    otherInteractors = []

    # Check if there is data
    if not mappings_data:
        # print("NA")
        mappings_data = ['NA']
        return 'none', 'none', 'none', 'none', 'none', 'none', 'none'
    else:
        # extract chains data
        # dict that contains uniprots
        uniprot = mappings_data[pdb.lower()]['UniProt']
        for uID in uniprot.keys():
            if uID != targetProtein and uID != interactorProtein:
                otherInteractors.append(uID)
            else:
                if uID == targetProtein:
                    uniprotTarget = uniprot[targetProtein]['mappings']
                    for mapping in uniprotTarget:
                        targetChainIds.append(mapping['chain_id'])
                        targetStart.append(mapping['unp_start'])
                        targetEnd.append(mapping['unp_end'])
                if uID == interactorProtein:
                    uniprotInteractor = uniprot[interactorProtein]['mappings']
                    for mapping in uniprotInteractor:
                        interactorChainIds.append(mapping['chain_id'])
                        interactorStart.append(mapping['unp_start'])
                        interactorEnd.append(mapping['unp_end'])

    # convert lists to strings
    targetChainIds = 'na' if targetChainIds == [] else ';'.join([str(x) for x in targetChainIds])
    targetStart = 'na' if targetStart == [] else ';'.join([str(x) for x in targetStart])
    targetEnd = 'na' if targetEnd == [] else ';'.join([str(x) for x in targetEnd])

    interactorChainIds = 'na' if interactorChainIds == [] else ';'.join([str(x) for x in interactorChainIds])
    interactorStart = 'na' if interactorStart == [] else ';'.join([str(x) for x in interactorStart])
    interactorEnd = 'na' if interactorEnd == [] else ';'.join([str(x) for x in interactorEnd])

    otherInteractors = 'na' if otherInteractors == [] else ';'.join([str(x) for x in otherInteractors])

    return targetChainIds, targetStart, targetEnd, interactorChainIds, interactorStart, interactorEnd, otherInteractors


def make_request(url, mode, pdb_id):
    """
    This function can make GET and POST requests to
    the PDBe API

    :param url: String,
    :param mode: String,
    :param pdb_id: String
    :return: JSON or None
    """
    if mode == "get":
        time.sleep(0.01)
        response = requests.get(url=url + pdb_id)
    elif mode == "post":
        time.sleep(0.01)
        response = requests.post(url, data=pdb_id)

    if response.status_code == 200:
        return response.json()
    else:
        print("NA from ", url, " for pdb ", pdb_id)

    return None

def normal_run(args):
    datasets = []
    filterSameProteinInteraction = False

    if args.filter:
        filterSameProteinInteraction = True

    # read data with pandas
    data = pd.read_csv(args.i, sep=';', converters={'Score': Decimal})

    # filtering for taxon.A = 9606 AND taxon.B = 9606 AND score >= cutoff (args.s)
    data = data[(data['Taxon A'] == 9606) & (data['Taxon B'] == 9606) & (data['Score'] >= args.s)]

    dataframeOut = pd.DataFrame(columns=['target uniprot id', 'target uniprot gene',  # 2 -> from csv
                                         'interactor uniprot id', 'interactor uniprot gene',  # 2 -> from csv
                                         'mentha score',  # 1 -> from csv
                                         'PDB id',  # 1 -> from RCSB API
                                         'fusion',  # 1 -> from summary request
                                         'target chain id', 'target starting residue', 'target ending residue',
                                         # 3 -> from mappings request
                                         'interactor chain id', 'interactor starting residue',
                                         'interactor ending residue',  # 3 -> from mappings request
                                         'other interactors',  # 1 -> from mappings request
                                         'method',  # 1 -> from summary request
                                         'resolution',  # 1 -> from experiment request
                                         'dna chains', 'num ligands'])  # 2 -> from summary request
    # -----------------------------
    # 18 columns total

    # open uniprot target file and get lines
    targets = []
    with open(args.t, 'r') as uniprotTargets:
        for uniprot in uniprotTargets:
            targets.append(uniprot)
            if args.x:
                dataframeOut = dataframeOut[0:0]
            uniprot = uniprot.rstrip()

            # get data where protein A or protein b matches uniprot selected
            uniprotData = data[(data['Protein A'] == uniprot) | (data['Protein B'] == uniprot)]

            uniprotData = uniprotData.reset_index()  # make sure indexes pair with number of rows

            targetQueryResult = get_pdb_entries_for_uniprot(uniprot)
            print('Target {}                                          '.format(uniprot))
            for index, row in uniprotData.iterrows():
                targetProtein = ''
                interactorProtein = ''
                targetGene = ''
                interactorGene = ''
                score = Decimal('0')
                # list that will later added to dataframe out
                outRow = []

                if row['Protein A'] == uniprot:
                    targetProtein = row['Protein A']
                    interactorProtein = row['Protein B']
                    targetGene = row['Gene A']
                    interactorGene = row['Gene B']
                else:
                    targetProtein = row['Protein B']
                    interactorProtein = row['Protein A']
                    targetGene = row['Gene B']
                    interactorGene = row['Gene A']

                # filter protein interaction with self
                if filterSameProteinInteraction and row['Protein A'] == row['Protein B']:
                    print('skipped protein interaction with self \n \t target {} interactor'.format(targetProtein,
                                                                                                    interactorProtein))
                    continue

                score = row['Score']

                # setup first 5 of outRow
                outRow.extend([targetProtein, targetGene, interactorProtein, interactorGene, score])

                # sending RCSB API requests
                interactorQueryResult = get_pdb_entries_for_uniprot(interactorProtein)

                # check if something went wrong in RCSB API -> set na and go next
                if not interactorQueryResult or not targetQueryResult:
                    # set output row to na (13 cause we had 5 set and 13 missing positions)
                    outRow.extend(['na'] * 13)
                    print(f'\t No PDB entries found via RCSB API for interactor {interactorProtein}         ', end='\r')
                    # append row to dataframe Out
                    dataframeOut.loc[len(dataframeOut)] = outRow
                else:
                    # get common pdbs to both proteins
                    commonPdbs = set(targetQueryResult).intersection(set(interactorQueryResult))

                    # if intersection is not empty
                    if commonPdbs != set():
                        print('\t protein interactor {} share pdbs -> {}'.format(interactorProtein, commonPdbs))
                        # intersection not empty -> run requests to get other columns
                        # do requests
                        # !! can be multiple pdbs !!
                        for pdb in commonPdbs:
                            fused = ''
                            dna = ''
                            ligands = ''
                            method = ''
                            targetChainIds = ''
                            targetStart = ''
                            targetEnd = ''
                            interactorChainIds = ''
                            interactorStart = ''
                            interactorEnd = ''
                            otherInteractors = ''
                            resolution = ''

                            # do requests on pdb chosen
                            # request summary -> from summary we get fusion, method,dna chains, num ligands
                            fused, dna, ligands, method = get_summary(pdb)
                            # request mappings -> from mappings we get chain infos (id, start, stop) and other interactors
                            targetChainIds, targetStart, targetEnd, interactorChainIds, interactorStart, interactorEnd, otherInteractors = get_mappings_data(
                                pdb, targetProtein, interactorProtein)
                            # request experiment -> from experiment we get resolution
                            resolution = get_experiment(pdb)

                            # add data to output row
                            outRow.extend([pdb, fused, targetChainIds, targetStart, targetEnd, interactorChainIds,
                                           interactorStart, interactorEnd, otherInteractors, method, resolution, dna,
                                           ligands])

                            # add out row to dataframe
                            # rowSerie = pd.Series(outRow, index = dataframeOut.columns)
                            # dataframeOut = dataframeOut.append(rowSerie, ignore_index=True)
                            dataframeOut.loc[len(dataframeOut)] = outRow

                            # reset outrow to first five values -> first five are fixed until we don't change target - interactor pair
                            # other values change basing on the pdb selected
                            outRow = outRow[:5]
                    else:
                        print('\t interactor {} ->  NO COMMON PDBS'.format(interactorProtein), end='\r')
                        # intersection empty -> set na and go on
                        # empty  pdb list, no requests set na
                        outRow.extend(['na'] * 13)
                        dataframeOut.loc[len(dataframeOut)] = outRow


            # args.x -> 1 csv per target
            if args.x:
                if args.p:
                    dataframeOutx = pmid_adder(data, dataframeOut)
                    # replace chars that will break to_csv
                    dataframeOutx.replace({',': '_'}, regex=True, inplace=True)

                    dataframeOutx.sort_values(['target uniprot id', 'mentha score'], ascending=False, inplace=True)

                    dataframeOutx['normal_or_cfg'] = 0

                    datasets.append(dataframeOutx)
                    print('>> Out for uniprot {} -> {}'.format(uniprot, 'dataframe_' + uniprot + '.csv'))
                else:
                    # replace chars that will break to_csv
                    dataframeOut.replace({',': '_'}, regex=True, inplace=True)

                    dataframeOut.sort_values(['target uniprot id', 'mentha score'], ascending=False, inplace=True)

                    dataframeOutx['normal_or_cfg'] = 0

                    datasets.append(dataframeOutx)

                # if option -a is selected we have to create folder and subfolders for input.fasta files
                if args.a:
                    make_target_interactor_sequence_files(dataframeOut)
            else:
                continue


    if not args.x:
        if args.p:
            dataframeOut = pmid_adder(data, dataframeOut)
        # replace chars that will break to_csv
        dataframeOut.replace({',': '_'}, regex=True, inplace=True)

        dataframeOut.sort_values(['target uniprot id', 'mentha score'], ascending=False, inplace=True)

        dataframeOut['normal_or_cfg'] = 0

        datasets.append(dataframeOut)

        if args.a:
            make_target_interactor_sequence_files(dataframeOut)

    return datasets, targets

def cfg_run(args):
    print('CFG')
    datasets = []

    config_file = args.c

    # if have config, read config
    if config_file != '':
        config_dict = {}
        import configparser

        config = configparser.ConfigParser()
        with open(config_file, 'r') as cfg_file:
            config.read_file(cfg_file)

        for each_section in config.sections():
            l = []
            for (each_key, each_val) in config.items(each_section):
                l.append(each_val.strip().split(','))
            config_dict[each_section] = l
    else:
        return None

    filterSameProteinInteraction = False

    if args.filter:
        filterSameProteinInteraction = True

    # read data with pandas
    data = pd.read_csv(args.i, sep=';', converters={'Score': Decimal})

    # filtering for taxon.A = 9606 AND taxon.B = 9606 AND score >= cutoff (args.s)
    data = data[(data['Taxon A'] == 9606) & (data['Taxon B'] == 9606) & (data['Score'] >= args.s)]

    dataframeOut = pd.DataFrame(columns=['target uniprot id', 'target uniprot gene',  # 2 -> from csv
                                         'interactor uniprot id', 'interactor uniprot gene',  # 2 -> from csv
                                         'mentha score',  # 1 -> from csv
                                         'PDB id',  # 1 -> from RCSB API
                                         'fusion',  # 1 -> from summary request
                                         'target chain id', 'target starting residue', 'target ending residue',
                                         # 3 -> from mappings request
                                         'interactor chain id', 'interactor starting residue',
                                         'interactor ending residue',  # 3 -> from mappings request
                                         'other interactors',  # 1 -> from mappings request
                                         'method',  # 1 -> from summary request
                                         'resolution',  # 1 -> from experiment request
                                         'dna chains', 'num ligands'])  # 2 -> from summary request
    # -----------------------------
    # 18 columns total

    # open uniprot target file and get lines
    with open(args.t, 'r') as uniprotTargets:
        pmids = []
        for uniprot in uniprotTargets:


            if args.x:
                dataframeOut = dataframeOut[0:0]
                pmids = []

            uniprot = uniprot.rstrip()

            # get data where protein A or protein b matches uniprot selected
            uniprotData = data[(data['Protein A'] == uniprot) | (data['Protein B'] == uniprot)]

            uniprotData = uniprotData.reset_index()  # make sure indexes pair with number of rows

            # cfg_data = config_dict[uniprot]
            cfg_data = config_dict.get(uniprot, [])
            if cfg_data == []:
                print(f'no config for target {uniprot} ')

            for cfg in cfg_data:
                int_id = cfg[0]
                int_gene = cfg[1]
                int_pdb = cfg[2]
                int_pmid = cfg[3]

                pmids.append(int_pmid)

                targetProtein = ''
                interactorProtein = ''
                targetGene = ''
                interactorGene = ''
                score = 'na'

                outRow = []

                for index, row in uniprotData.iterrows():
                    if row['Protein A'] == uniprot:
                        targetProtein = row['Protein A']
                        interactorProtein = int_id
                        targetGene = row['Gene A']
                        interactorGene = int_gene
                        break
                    else:
                        targetProtein = row['Protein B']
                        interactorProtein = int_id
                        targetGene = row['Gene B']
                        interactorGene = int_gene
                        break

                #try to get score from db for config lines
                for index, row in uniprotData.iterrows():
                    if row['Protein A'] == uniprot and row['Protein B'] == int_id:
                        score = row['Score']
                    elif row['Protein A'] == int_id and row['Protein B'] == uniprot:
                        score = row['Score']

                # setup first 5 of outRow
                outRow.extend([targetProtein, targetGene, interactorProtein, interactorGene, score])

                if int_pdb == '':
                    print(f'>>PDB is not present \n\t {cfg} \n row is na')
                    #setting score to na -> we have no score coming from mentha db -> target interactor not in db
                    outRow[4] = 'na'
                    #setting all columns to na -> no pdb to use for requests
                    ext = ['na' for x in range(13)]
                    outRow.extend(ext)
                    #save row
                    dataframeOut.loc[len(dataframeOut)] = outRow
                    #reset row for next config row
                    outRow = outRow[:5]
                    continue

                #we reach this part only if we have a pdb to use
                fused = ''
                dna = ''
                ligands = ''
                method = ''
                targetChainIds = ''
                targetStart = ''
                targetEnd = ''
                interactorChainIds = ''
                interactorStart = ''
                interactorEnd = ''
                otherInteractors = ''
                resolution = ''

                # do requests on pdb chosen
                # request summary -> from summary we get fusion, method,dna chains, num ligands
                fused, dna, ligands, method = get_summary(int_pdb)
                # request mappings -> from mappings we get chain infos (id, start, stop) and other interactors
                targetChainIds, targetStart, targetEnd, interactorChainIds, interactorStart, interactorEnd, otherInteractors = get_mappings_data(
                    int_pdb, targetProtein, interactorProtein)
                # request experiment -> from experiment we get resolution
                resolution = get_experiment(int_pdb)

                # add data to output row
                outRow.extend([int_pdb, fused, targetChainIds, targetStart, targetEnd, interactorChainIds,
                               interactorStart, interactorEnd, otherInteractors, method, resolution, dna,
                               ligands])

                # add out row to dataframe
                dataframeOut.loc[len(dataframeOut)] = outRow
                outRow = outRow[:5]

            if args.x:
                if args.p:
                    #add pmid col
                    df_out = dataframeOut.copy(deep=True)
                    df_out['PMID'] = pmids
                    df_out.replace({',': '_'}, regex=True, inplace=True)
                    df_out.sort_values(['target uniprot id', 'mentha score'], ascending=False, inplace=True)
                    df_out['normal_or_cfg'] = 1
                    datasets.append(df_out)
                else:
                    df_out = dataframeOut.copy(deep=True)
                    df_out.replace({',': '_'}, regex=True, inplace=True)
                    df_out.sort_values(['target uniprot id', 'mentha score'], ascending=False, inplace=True)
                    df_out['normal_or_cfg'] = 1
                    datasets.append(df_out)
                if args.a:
                    make_target_interactor_sequence_files(df_out)

    if not args.x:
        if args.p:

            dataframeOut['PMID'] = pmids
        dataframeOut.replace({',': '_'}, regex=True, inplace=True)
        dataframeOut.sort_values(['target uniprot id', 'mentha score'], ascending=False, inplace=True)
        dataframeOut['normal_or_cfg'] = 1
        datasets.append(dataframeOut)
        if args.a:
            make_target_interactor_sequence_files(dataframeOut)

    return datasets

def extract_genes(data, edf_list, target_list):
    ol = []

    for target, edf in zip(target_list, edf_list):
        #target gene
        tg = extract_helper(data, target)
        edf['target uniprot gene'] = tg

        #interactor gene
        edf_interactors = edf['interactor uniprot id'].tolist()
        edf_interactor_gene_list = []
        for interactor in edf_interactors:
            ig = extract_helper(data, interactor)
            edf_interactor_gene_list.append(ig)

        edf['interactor uniprot gene'] = edf_interactor_gene_list

        ol.append(edf)

    return ol

def extract_helper(data, id):
    gene = ''
    targetdata = data[(data['Protein A'] == id)]
    if not targetdata.empty:
        gene = targetdata['Gene A'].iloc[0]
    else:
        targetdata = data[(data['Protein B'] == id)]
        if not targetdata.empty:
            gene = targetdata['Gene B'].iloc[0]
        else:
            #last chance make request to uniprot.org
            gene = extract_gene_fromrequest(id)

    return gene

def extract_gene_fromrequest(id):

    url = 'https://rest.uniprot.org/uniprotkb/search?query='
    res = make_request(url,'get',id)
    gene = res['results'][0]['genes'][0]['geneName']['value']

    return gene

def copy_folder(ex, id1, id2, af_folder_path):
    #ex -> extra file name
    #id1 id2 -> pair components

    #clean path before extra file name and get only the name no extension
    ex = os.path.basename(ex)
    ex = ex.split('.')[0]

    folder = id1 + '-' + id2

    from_path = af_folder_path
    to_path = 'AF_Huri_HuMAP'
    if 'huri' in ex.lower():
        from_path = Path(from_path).joinpath('Huri_dimers').joinpath('HuRI').joinpath(folder)
        to_path = Path(to_path).joinpath('Huri_dimers').joinpath(folder)
    elif 'humap' in ex.lower():
        from_path = Path(from_path).joinpath('HuMAP_dimers').joinpath('pdb').joinpath(folder)
        to_path = Path(to_path).joinpath('HuMAP_dimers').joinpath(folder)

    if from_path.exists() and not to_path.exists():
        ignore_func = lambda d, files: [f for f in files if isfile(join(d, f)) and not f.endswith('.pdb')]
        shutil.copytree(from_path, to_path, ignore=ignore_func)
        print(f'>>>copy from path \n {from_path} \n to \n {to_path}')
    else:
        s = f'destination path already exists \n {to_path}' \
            if to_path.exists() else\
            f'source folder does not exist \n {from_path}'
        print(s)

def rename_pair_folder_direct(ensg1, ensg2, up1, up2, base_path='AF_Huri_HuMAP'):
    
    huri_path = Path(base_path, 'Huri_dimers')
    old_folder = huri_path / f"{ensg1}-{ensg2}"
    new_folder = huri_path / f"{up1}-{up2}"
    
    try:
        if not old_folder.exists():
            raise FileNotFoundError(f"Expected old folder missing: {old_folder}")
        if new_folder.exists():
            raise FileExistsError(f"Target already exists: {new_folder}")
        
        old_folder.rename(new_folder)
        print(f"Renamed {old_folder} → {new_folder}")

    except FileNotFoundError as e:
        print(f"ERROR: required folder not found.\n{e}", file=sys.stderr)
        sys.exit(1)
    except FileExistsError as e:
        print(f"ERROR: output file(s) already exist — please remove them and try again.\n{e}", file=sys.stderr)
        sys.exit(1)

def process_extra_files(args, extra_files):

    datasets = []
    extra_df = pd.DataFrame(columns=['target uniprot id', 'target uniprot gene',  # 2 -> from csv
                                         'interactor uniprot id', 'interactor uniprot gene',  # 2 -> from csv
                                         'mentha score',  # 1 -> from csv
                                         'PDB id',  # 1 -> from RCSB API
                                         'fusion',  # 1 -> from summary request
                                         'target chain id', 'target starting residue', 'target ending residue',
                                         # 3 -> from mappings request
                                         'interactor chain id', 'interactor starting residue',
                                         'interactor ending residue',  # 3 -> from mappings request
                                         'other interactors',  # 1 -> from mappings request
                                         'method',  # 1 -> from summary request
                                         'resolution',  # 1 -> from experiment request
                                         'dna chains', 'num ligands','PMID'])  # 2 ->



    # read data with pandas
    data = pd.read_csv(args.i, sep=';', converters={'Score': Decimal})
    # filtering for taxon.A = 9606 AND taxon.B = 9606 AND score >= cutoff (args.s)
    data = data[(data['Taxon A'] == 9606) & (data['Taxon B'] == 9606) & (data['Score'] >= args.s)]

    targets = []
    file = open(args.t, 'r')
    targets = file.readlines()
    targets = [t.strip() for t in targets]
    file.close()

    if args.extra == [] or args.extra == None:
        print('No extra files given, skipping extra files processing')
    else:

        for e in extra_files:
            extra_df[e] = []

        pairs_scores = []
        for extra_file in extra_files:
            filename = os.path.basename(extra_file).lower()
            pair_score=[]
            extra_file_data = pd.read_csv(extra_file, sep=',')
            #cut all scores under cutoff
            extra_file_data = extra_file_data[extra_file_data.pDockQ >= args.extra_cutoff]

            for _, r in extra_file_data.iterrows():
                up1, up2 = r['NameUPAC'].split('-', 1)
                score = r['pDockQ']
                if "huri" in filename:
                    ensg1, ensg2 = r['Name'].split('-', 1)
                    pair_score.append([ensg1, ensg2, up1, up2, score])
                else:
                    pair_score.append([None, None, up1, up2, score])
            pairs_scores.append([extra_file, pair_score])

        edf = []
        for target in targets:
            df_t = []
            for i, e_ps in enumerate(pairs_scores):
                ex = e_ps[0]
                ex_name = os.path.basename(ex).lower()
                ps = e_ps[1]
                for p in ps:
                    ensg1, ensg2, up1, up2, score = p
                    g1 = 'extra gene'
                    g2 = 'extra gene'
                    if up1 == target:
                        row = [up1, g1, up2, g2, 'na'] + ['na']*14

                        row = row + ['na']*i +[score]+['na']*(len(args.extra) -i -1)

                        extra_df.loc[len(extra_df)] = row
                        if "huri" in ex_name:
                            copy_folder(ex, ensg1, ensg2, args.af)
                            rename_pair_folder_direct(ensg1, ensg2, up1, up2)
                        elif "humap" in ex_name:
                            copy_folder(ex, up1, up2, args.af)

                    elif up2 == target:
                        row = [up2, g2, up1, g1, 'na'] + ['na']*14

                        row = row + ['na']*i +[score]+['na']*(len(args.extra) -i -1)

                        extra_df.loc[len(extra_df)] = row
                        if "huri" in ex_name:
                            copy_folder(ex, ensg1, ensg2, args.af)
                            rename_pair_folder_direct(ensg1, ensg2, up1, up2)
                        elif "humap" in ex_name:
                            copy_folder(ex, up1, up2, args.af)
                df_t.append(extra_df)
                extra_df = extra_df[0:0]

            df = pd.concat([d for d in df_t])
            edf.append(df)

        # grab gene from bs for extra files
        datasets = extract_genes(data, edf, targets)

        if not args.x:
            datasets = [pd.concat([d for d in datasets])]

    return datasets



def grab_result(url):
    response = session.get(url)
    #logging.info("request was completed in %s seconds [%s]", response.elapsed.total_seconds(), response.url)
    if response.status_code != 200:
        pass
        #logging.error("request failed, error code %s [%s]", response.status_code, response.url)
    if 500 <= response.status_code < 600:
        # server is overloaded? give it a break
        time.sleep(5)
    return response

def download(urls, d):
    with ThreadPoolExecutor(max_workers=THREAD_POOL) as executor:
        # wrap in a list() to wait for all requests to complete
        for response, url in zip(list(executor.map(grab_result, urls)), urls):
            _, ensg_number = url.split('=', 1)
            if response.status_code == 200:
                try:
                    r = json.loads(response.text)
                except Exception as e:
                    pass
                if r['results'] != []:
                    d[ensg_number] = r['results'][0]['primaryAccession']
                else:
                    d[ensg_number] = None
            else:
                d[ensg_number] = None
    return d


def main(argv):
    warnings.filterwarnings("ignore")
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--i', help='mentha database file')
    parser.add_argument('-t', '--t', help='File with target uniprots')
    parser.add_argument('-s', '--s', type=Decimal, help='Cutoff score')
    parser.add_argument('-o', '--o', nargs='?', const='dataframe.csv', default='dataframe.csv', help='Output name')
    parser.add_argument('-f', '--filter', action='store_true')
    parser.add_argument('-p', '--p', action='store_true', help='option to add PMID column to output')
    parser.add_argument('-x', '--x', action='store_true', help='option to have 1 csv output file per target uniprot ID')
    parser.add_argument('-a', '--a', action='store_true', help='option to have inputs_afmulti folder with subfolders and input.fasta files')
    parser.add_argument('-c', '--c', default='', help='Config file containing rows to insert into mentha db')
    parser.add_argument('-extra', '--extra-files', dest='extra', nargs='*', required=False, default=None, help='list of extra files to process')
    parser.add_argument('-ec','--extra-cutoff', dest='extra_cutoff', default=0.5, type=float, help='Cutoff on extra files pair pDockQ scores')
    parser.add_argument('-af','--af-folder', dest='af', help='AF_Huri_HuMAP folder location')

    args = parser.parse_args()

    if args.extra != None and args.af == None:
        print('Detected extra files but no AF_Huri_HuMAP folder path, use the -af parameter')
        print('quitting.')
        sys.exit(0)

    datasets = []
    config_datasets = []
    extra_datasets = []
    targets = []

    datasets, targets = normal_run(args)
    config_datasets = cfg_run(args)
    extra_datasets = process_extra_files(args, args.extra)

    if config_datasets == [] or config_datasets == None:
        config_datasets = []
        for d in datasets:
            config_datasets.append(pd.DataFrame(columns=d.columns))

    if extra_datasets == [] or extra_datasets == None:
        extra_datasets = []
        for d in datasets:
            extra_datasets.append(pd.DataFrame(columns=d.columns))

    for ds, ds_cfg, ds_extra, target in zip(datasets, config_datasets, extra_datasets, targets):
        result = pd.merge(ds, ds_cfg, how='outer',
                          left_on=['target uniprot id', 'target uniprot gene', 'interactor uniprot id', 'interactor uniprot gene', 'PDB id'],
                          right_on=['target uniprot id', 'target uniprot gene', 'interactor uniprot id', 'interactor uniprot gene', 'PDB id'])

        result.replace('na', np.nan, inplace=True)

        l = []
        dfxF = pd.DataFrame(columns=ds.columns)
        for i, row in result.iterrows():
            if row['normal_or_cfg_x'] == 0.0 and pd.isna(row['normal_or_cfg_y']):
                # row in normal run but not in config
                lfix = row.iloc[0:20].tolist()
                l.append(lfix)
                dfxF.loc[len(dfxF)] = lfix
            elif row['normal_or_cfg_x'] == 0.0 and row['normal_or_cfg_y'] == 1.0:
                # row in normal and in config
                row['pmid'] = str(row['PMID_x']) + ' ' + str(row['PMID_y'])
                lfix = row.iloc[0:18].tolist() + [row['pmid'], 2]
                l.append(lfix)
                dfxF.loc[len(dfxF)] = lfix
            elif pd.isna(row['normal_or_cfg_x']) and row['normal_or_cfg_y'] == 1.0:
                # row in cfg but not in normal
                lfix = row.iloc[0:6].tolist() + row.iloc[21:].tolist()
                l.append(lfix)
                dfxF.loc[len(dfxF)] = lfix


        dfxF.drop(['normal_or_cfg'], axis=1, inplace=True)

        if args.extra != [] and args.extra != None:
            for e in args.extra:
                dfxF[e] = 'na'

        for i,r in ds_extra.iterrows():
            index_list = []
            index_list = dfxF[(dfxF['target uniprot id'] == r['target uniprot id']) &
                              (dfxF['target uniprot gene'] == r['target uniprot gene']) &
                              (dfxF['interactor uniprot id'] == r['interactor uniprot id']) &
                              (dfxF['interactor uniprot gene'] == r['interactor uniprot gene'])
                            ].index.tolist()
            if index_list  != []:
                #ds extra row already in dataframe
                for e in args.extra:
                    dfxF.loc[index_list, e] = r[e]

            else:
                dfxF.loc[len(dfxF)] = r

        if args.extra == None:
            # for e in args.extra:
            #     dfxF.drop([e], axis=1, inplace=True)
            pass
        else:
            # fix col names
            columns = list(dfxF.columns)
            columns_to_fix = columns[-len(args.extra):]
            col_fix = [os.path.basename(c).split('.', 1)[0] for c in columns_to_fix]
            columns = columns[:-len(args.extra)] + col_fix
            dfxF.columns = columns
        # Define new column names
        new_last_column_name = "pDockQ HuRI"
        new_second_last_column_name = "pDockQ HuMap"

        # Rename the last two columns
        dfxF.columns.values[-1] = new_last_column_name
        dfxF.columns.values[-2] = new_second_last_column_name

        dfxF.sort_values(['target uniprot id', 'mentha score', 'interactor uniprot id', 'PDB id'], ascending=False, inplace=True)
        dfxF.replace(np.nan, 'na', inplace=True)

        csv_outname = args.o
        if len(datasets) == 1:
            pass
        else:
            splitted_o = args.o.split('.')
            #example
            #args.o = out.csv
            #csv_outname = out_<target>.csv
            csv_outname = f'{splitted_o[0]}_{target.strip()}.csv'
            #csv_outname = f'dataframe_{target.strip()}.csv'

        if not args.x:
            print(f'>>writing full dataframe (no splitted option selected -x) -> {csv_outname}')
            dfxF.to_csv(csv_outname, index=False, quoting=csv.QUOTE_NONE, sep=',')
        else:
            print(f'>>writing dataframe for target {target} -> {csv_outname}')
            dfxF.to_csv(csv_outname, index=False, quoting=csv.QUOTE_NONE, sep=',')


THREAD_POOL = 16

# This is how to create a reusable connection pool with python requests.
session = requests.Session()
session.mount(
    'https://rest.uniprot.org/uniprotkb/search?query=',
    requests.adapters.HTTPAdapter(pool_maxsize=THREAD_POOL,
                                  max_retries=3,
                                  pool_block=True)
)


if __name__ == "__main__":
    main(sys.argv[1:])


