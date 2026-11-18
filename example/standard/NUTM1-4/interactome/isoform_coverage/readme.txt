The output of this script is a copy of the input table with additional DIGGER annotation columns.
It reports:
- the reference and isoform sequence lengths
- the mapped ENST used for the target
- whether target domains are unchanged, lost, gained, or rewired
- which DIGGER DDIs are retained or missing
- the percentage of missing DDIs
- the predicted PPI effect (Retained, Affected, or Uncertain)
- the DIGGER confidence label

Run:
python annotate_isoforms_digger.py <aggregated_csv> <target_isoform_uniprot> <output_csv>

Example:
module load python

tsp -N 1 python annotate_isoforms_digger.py ../aggregate/Q9NQX3-2_aggregated.csv Q9NQX3-2 Q9NQX3-2_aggregated_isoform.csv
