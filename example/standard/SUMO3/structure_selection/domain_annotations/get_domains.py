#!/usr/bin/env python3

# Copyright (C) 2022 Ludovica Beltrame
# Copyright (C) 2024 Ludovica Beltrame, Kristine Degn, Lorenzo Favaro

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

# standard library modules
import sys, errno, re, json, ssl
from urllib import request
from urllib.error import HTTPError
from time import sleep
import argparse
import pandas as pd
import numpy as np
import sys

parser = argparse.ArgumentParser(description="Retrieve pfam protein domains")

u_helpstr = 'uniprotID'
parser.add_argument("-u", "--uniprot", type=str, required=True, help=u_helpstr)

m_helpstr = 'Retrieve domains of point mutation sites'
# Optional argument for mutation file
parser.add_argument("-m", "--mutation-file", help=m_helpstr)

# Check if the -m flag is used and if the mutation file is provided
args = parser.parse_args()

def output_lists():
    #disable SSL verification to avoid config issues
    context = ssl._create_unverified_context() 
  
    next = BASE_URL 
    last_page = False
    
    attempts = 0
    while next:
      try:
        req = request.Request(next, headers={"Accept": "application/json"})
        res = request.urlopen(req, context=context)
        # If the API times out due a long running query
        if res.status == 408:
          # wait just over a minute
          sleep(61)
          # then continue this loop with the same URL
          continue
        elif res.status == 204:
          #no data so leave loop
          break
        payload = json.loads(res.read().decode())
        next = payload["next"]
        attempts = 0
        if not next:
          last_page = True
      except HTTPError as e:
        if e.code == 408:
          sleep(61)
          continue
        else:
          # If there is a different HTTP error, it wil re-try 3 times before failing
          if attempts < 3:
            attempts += 1
            sleep(61)
            continue
          else:
            sys.stderr.write("LAST URL: " + next)
            raise e
      
      # Don't overload the server, give it time before asking for more
      if next:
        sleep(1)
  
    #Store information in lists
    start = []
    end = []
    pfam_domain = []
    accession = []
    
    # Early return if 'results' is not in payload
    if not 'results' in payload:
        return start, end, pfam_domain, accession

    for item in payload['results']:
        # Check if the expected keys are present
        if 'proteins' in item and item['proteins'] and \
           'entry_protein_locations' in item['proteins'][0] and item['proteins'][0]['entry_protein_locations'] and \
           'fragments' in item['proteins'][0]['entry_protein_locations'][0] and item['proteins'][0]['entry_protein_locations'][0]['fragments']:
            
            for fragment in item['proteins'][0]['entry_protein_locations'][0]['fragments']:
                start.append(fragment['start'])
                end.append(fragment['end'])
                pfam_domain.append(item['metadata']['name'])
                accession.append(item['metadata']['accession'])

    return start, end, pfam_domain, accession

def domains_to_mutations(x):
    '''Assign domain to point mutation'''
    
    #If a domain and an accession number are available,
    #store them in the domain and accession variables
    try:
        domain = data.loc[(data['start'] <= x ) & (data['end'] >= x)]['pfam_domain'].values[0]
        accession = data.loc[(data['start'] <= x ) & (data['end'] >= x)]['accession'].values[0]
    
    #If the residue doesn't belong to any domain
    #return not found and NA
    except IndexError:
        domain = 'Not found'
        accession = 'NA'
    
    return pd.Series([domain, accession])

if __name__ == "__main__":
  
    uniprotID = args.uniprot
    
    BASE_URL = "https://www.ebi.ac.uk:443/interpro/api/entry/all/pfam/protein/UniProt/" + uniprotID + "/?page_size=200"
    
    #Interact with API, store output in a dictionary and
    #extract information in lists (start and end positions, 
    #domain name and accession number)
    try:
        start, end, pfam_domain, accession = output_lists()
    except UnboundLocalError:
        print('WARNING: no domains found annotated in pfam. Exiting...')
        start, end, pfam_domain, accession = [], [], [], []
    
    # Convert lists to dataframe
    if start and end and pfam_domain and accession:
        data = pd.DataFrame(list(zip(start, end, pfam_domain, accession)),
                     columns =['start', 'end', 'pfam_domain', 'accession'])
        data = data.sort_values(by=['start'], ascending=True)
    else:
        # Create an empty dataframe with the expected columns
        data = pd.DataFrame(columns=['start', 'end', 'pfam_domain', 'accession'])
    
    # Save summary.csv containing a summary of the protein domains
    data.to_csv('summary.csv', sep=';', index=False, columns=['start', 'end','pfam_domain', 'accession'])
     
    # If mutlist is provided do the following
    if args.mutation_file:
        # Process mutation file
        mutation_file_path = args.mutation_file
        # Read mutlist file in dataframe and add position column
        out = pd.read_csv(mutation_file_path, names=['mutation'])
        out['pos'] = out['mutation'].str[1:-1]  
        
        data["start"] = pd.to_numeric(data["start"], errors='coerce')
        data['end'] = pd.to_numeric(data['end'], errors='coerce')
        out['pos'] = pd.to_numeric(out['pos'], errors='coerce')
        
        # Get domains of point mutations            
        out[['pfam_domain', 'accession']] = out['pos'].apply(domains_to_mutations)
    
        # Save mut_domains.csv
        out.to_csv('domains_mutlist.csv', index=False, sep=';', columns=['mutation', 'pfam_domain', 'accession'])
    else:
        # Create an empty domains_mutlist.csv if no mutation file is provided
        empty_out = pd.DataFrame(columns=['mutation', 'pfam_domain', 'accession'])
        empty_out.to_csv('domains_mutlist.csv', index=False, sep=';')
