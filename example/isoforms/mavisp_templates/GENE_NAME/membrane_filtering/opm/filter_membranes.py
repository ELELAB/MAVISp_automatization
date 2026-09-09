#!/usr/bin/env python3

# Copyright (C) 2023 Mattia Utichi <biomatt90@gmail.com>
# Danish Cancer Society and Technical University of Denmark

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
import pandas as pd
import logging as log
import sys
import re

# Argument parser
desc = 'The script fiters a mutation list based on ' \
	   'a list of positions of residues in contact with the membranes. '\
	   'This list can be the DREAMM output (option -d), the OPM/PPM '\
	   'output (option -p) or the simulation analysis output (option -s)'
parser = argparse.ArgumentParser(description= desc)

parser.add_argument('-i', '--input', 
					type = str,
					required = True,
					help = 'input file with the positions to be filtered out')
parser.add_argument('-m', '--mutlist', 
					type = str,
					required = True,
					help = 'mutation list to be filtered')
parser.add_argument('-o', '--output',
					type = str,
					dest = 'out',
					default = 'filtered_mutlist.txt')

#add mutually exclusive options for the input source
source = parser.add_mutually_exclusive_group()
source.add_argument('-s','--sim-input', 
					dest = 'sim',
					action = 'store_true',
					default = False,
					help = 'if the input file comes from an '\
						   'MD simulation analysis')
source.add_argument('-p','--opm-ppm-input',
					dest = 'opm',
					action = 'store_true',
					default = False,
					help = 'if the input file comes from the OPM/PPM analysis')
source.add_argument('-d','--dreamm-input',
					dest = 'dreamm',
					action = 'store_true',
					default = False,
					help = 'if the input file comes from the DREAMM analysis')
args = parser.parse_args()

#If the input file comes from a simulation analysis
if args.sim is not False:

	# Read the simulation file
	try:
		data = pd.read_csv(args.input, usecols=['columns20percent'])
	except FileNotFoundError:
		log.error(f'File {args.input} not found. Exiting...')
		sys.exit(1)

	# Create the control list with the positions that need to be fitered out
	# as numerical strings
	control_list = data['columns20percent'].to_list()
	control_list = [x[:-2] for x in control_list]
	control_list = set(control_list)


#If the input file comes from DREAMM server
if args.dreamm is not False:

	# Read the DREAMM output csv
	try:
		data = pd.read_csv(args.input)
	except FileNotFoundError:
		log.error(f'File {args.input} not found. Exiting...')
		sys.exit(1)

	# Create the control list with the positions that need to be fitered out
	# as numerical strings
	control_list = data['resnum'].to_list()
	control_list = [str(x) for x in control_list]

if args.opm is not False:
	# Read the OPM/PPM output
	try:
		with open(args.input,'r') as f:
			control_list = []

			#read every line and get the residue list (positions and ranges)
			for line in f:

				#skip header
				if re.match(r'#',line):
					continue
				numbers = re.findall(r'[0-9\-]+',line)

				#remove the first number which is the total number
				numbers.pop(0)

				#generate all positions from ranges and append to control_list
				for number in numbers:
					res_range =  re.match(r'(\d+)-(\d+)',number)
					if res_range is not None:
						for i in range(int(res_range.group(1)),
									   int(res_range.group(2))+1):
							control_list.append(i)

					#append single positions
					else:
						control_list.append(number)

			#get the unique position list
			control_list = set(control_list)
			control_list = [str(x) for x in control_list]

	#if the input is missing or with the wrong name
	except FileNotFoundError:
		log.error(f'File {args.input} not found. Exiting...')
		sys.exit(1)

# Read the provided mutlist and initialize the filtered list
try:
	with open(args.mutlist,'r') as f:
		filtered = []

		# Fo every row of the mutlist check if the position is in the string 
		# control list. If it is not append it to the filtered list.
		for line in f:
			new_line = line.strip('\n')
			match = re.findall(r'\d+',new_line)
			if not match[0] in control_list:
				filtered.append(new_line)
except FileNotFoundError:
	log.error(f'File {args.mutlist} not found. Exiting...')
	sys.exit(1)

# Write the output filtered mutlist 
with open(args.out, 'a') as f:
	for i in filtered:
		f.write(i+'\n')
