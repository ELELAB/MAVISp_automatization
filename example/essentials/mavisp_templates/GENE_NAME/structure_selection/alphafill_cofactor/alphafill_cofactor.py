#!/usr/bin/env python3
# Copyright (C) 2023 Katrine Meldgård <katrine@meldgaard.dk>
# Cancer Structural Biology, Danish Cancer Society Research Center, 2100, Copenhagen, Denmark
# Cancer Systems Biology, Health and Technology Department, Section for Bioinformatics, 2800, Lyngby, Denmark

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

import argparse
import MDAnalysis as mda
import requests as rq
import sys
from cif2pdb import cif2pdb
import os
import shutil

#----------------------------------------------------
#            Creating parsing arguments
#----------------------------------------------------
parser=argparse.ArgumentParser(description='Find interacting residues for co-factors in alphafill entries.',prefix_chars='-')

required_arguments=parser.add_argument_group('required arguments')
required_arguments.add_argument('-u',type=str,required=True,help='UniProt ID')
parser.add_argument('-t', required=False,default=5,help='The A threshold for how finding the interacting residues')
parser.add_argument('-co',required=False,nargs='+',default=None,help='Co-factors to find interacting residues for. Please write them in pairs with the alphafil id and cofactor name (mind capital letters): B FAD D DTC')
parser.add_argument('-cha',required=False,nargs='+',default=None,help='The name of the chain(s), as encoded in alphafill. Should be written with a space between each chain (mind capital letters): A B')

#Execute parse_args()
args = parser.parse_args()

if args.co is not None:
    if len(args.co)%2 != 0:
        parser.print_help()
        parser.error('Argument -co is not in the correct format.')


#----------------------------------------------------
#         Collect and prepare alphafill entry
#----------------------------------------------------
#Get entry and metadata from alphafill
upid=args.u
entry_url = f'https://alphafill.eu/v1/aff/{upid}'
meta_url = f'https://alphafill.eu/v1/aff/{upid}/json'

entry_response = rq.get(entry_url)
if entry_response.status_code != 200:
    print(f'Alphafill did not respond properly to the request. Status code {entry_response.status_code}')
    sys.exit()

#Convert cif response to pdb
with open(f'{upid}.mmcif','w') as outfile:
    outfile.write(entry_response.text)
#Function saves pdb and returns chain mapping
new_to_org_map = cif2pdb(f'{upid}.mmcif',f'{upid}.pdb')
org_to_new_map = {y:x for x,y in new_to_org_map.items()}

meta_response = rq.get(meta_url)
if meta_response.status_code != 200:
    print(f'Alphafill did not respond properly to the request. Status code {meta_response.status_code}')
    sys.exit()


#----------------------------------------------------
#            Prepare universe and argumnets
#----------------------------------------------------
#Create universe from pdb
u = mda.Universe(f'{upid}.pdb')

#Find chains if none are specified
cha_flag=False
if args.cha is None:
    cha_flag=True
    args.cha = list()

id_co = dict()
cofactor = list()
for seg in u.segments:
    #Finds cofactor if none are specified
    if  len(seg.residues) == 1:
        id_co[seg.segid] = seg.residues[0].resname
        if args.co is None:
            cofactor.append(seg.segid)
    elif cha_flag and len(seg.residues) > 1:
        args.cha.append(new_to_org_map[seg.segid])


if args.co is not None: 
    for i in range(0,len(args.co),2):
        #Check that provided alphafill id matches provided cofactor name
        if id_co[org_to_new_map[args.co[i]]] != args.co[i+1]:
            parser.print_help
            raise KeyError('The alphafill id does not match the paired cofactor name.')
        else:
            cofactor.append(org_to_new_map[args.co[i]])



#----------------------------------------------------
#              Find interacting residues
#----------------------------------------------------
inter_res = dict()
for chain in args.cha:
    for seg in u.segments:
        #if the segment is not a chain and is in the list of cofactor
        if seg.segid not in args.cha and seg.segid in cofactor:
            atoms = u.select_atoms(f'segid {chain} and around {args.t} segid {seg.segid}')
            if seg.segid not in inter_res:
                inter_res[seg.segid] = str(set(atoms.resids))
            else:
                inter_res[seg.segid] += ', '+str(set(atoms.resids))


#----------------------------------------------------
#             Find interacting cofactor
#----------------------------------------------------
inter_co = dict()
for co in cofactor:
    for seg in u.segments:
        if seg.segid not in args.cha:
            atoms = u.select_atoms(f'segid {co} and around {args.t} segid {seg.segid}')
            if len(atoms) > 0:
                if co not in inter_co:
                    inter_co[co] = [[new_to_org_map[seg.segid], id_co[seg.segid]]]
                else:
                    inter_co[co].append([new_to_org_map[seg.segid], id_co[seg.segid]])


#----------------------------------------------------
#           Write results to output files
#----------------------------------------------------
header_string = 'cofactor_af_id, cofactor_name'
for chain in args.cha:
    header_string = header_string + f', interacting_residues_chain{chain}'
with open('inter_res.csv','w') as outfile:
    outfile.write(header_string + '\n')
    for key in inter_res.keys():
        inter_res[key] = inter_res[key].replace('set()','NA')
        outfile.write(f'{new_to_org_map[key]}, {id_co[key]}, {inter_res[key]}\n')

header_string = 'cofactor_af_id, cofactor_name, interacting_cofactor(s)([[af_id, cofactor_name]])'
with open('inter_cofactors.csv','w') as outfile:
    outfile.write(header_string + '\n')
    for key in sorted(inter_co.keys()):
        outfile.write(f'{new_to_org_map[key]}, {id_co[key]}, {inter_co[key]}\n')

if os.path.exists("__pycache__"):
    # Remove "__pycache__" directory
    shutil.rmtree("__pycache__")
