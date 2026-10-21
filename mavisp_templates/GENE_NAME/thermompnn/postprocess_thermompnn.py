#!/usr/bin/env python3

# Copyright (C) 2026 Eszter Toldi, Technical University of Denmark
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""
Postprocess ThermoMPNN inference CSV

Applies a residue numbering offset and renames columns for ingestion in MAVISp. 

Usage:
    python postprocess_thermompnn.py --input <csv> --output <csv> --offset <int>
"""

import argparse
import pandas as pd
import sys

def postprocess(input_path, output_path, offset):

    required_cols = ['position', 'wildtype', 'mutation', 'ddG_pred']

    try:
        df = pd.read_csv(input_path, usecols=required_cols)
    except Exception as e:
        print(f"Error reading {input_path}: {type(e).__name__}: {e.args}")
        sys.exit(1)

    missing = [c for c in required_cols if c not in df.columns]
    if missing:
        print(f"Missing expected columns in {input_path}: {missing}")
        sys.exit(1)

    df['variant'] = (
        df['wildtype'].astype(str)
        + (df['position'] + offset).astype(str)
        + df['mutation'].astype(str)
    )

    if df['variant'].duplicated().any():
        dupes = df[df['variant'].duplicated(keep=False)]['variant'].unique().tolist()
        print(f"Duplicate variant labels after applying offset {offset}: {dupes}")
        print("Check that the offset value is correct.")
        sys.exit(1)

    out = df[['variant', 'ddG_pred']].rename(columns={'ddG_pred': 'ThermoMPNN_ddG'})
    out = out.set_index('variant')

    try:
        out.to_csv(output_path)
    except Exception as e:
        print(f"Error writing {output_path}: {type(e).__name__}: {e.args}")
        sys.exit(1)

    print(f"Written {len(out)} variants to {output_path}")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Postprocess ThermoMPNN inference CSV for MAVISp ingestion")
    parser.add_argument('--input',  required=True, help="Raw ThermoMPNN inference CSV")
    parser.add_argument('--output', required=True, help="Output normalized CSV path")
    parser.add_argument('--offset', required=True, type=int,
                        help="Integer offset added to zero-based position column")
    args = parser.parse_args()

    postprocess(args.input, args.output, args.offset)
