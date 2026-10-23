# Copy custom_inference_standalone.py, do.sh and readme.txt from
# mavisp_templates/GENE_NAME/thermompnn into the protein folder with the trimmed PDB.
# Run from that folder on a GPU server (either bioinfo03 or bioinfor04).

# Activate the ThermoMPNN environment:
conda activate /usr/local/envs/thermompnn

# Check GPU availability:
nvidia-smi

# Run:
tsp bash do.sh P60484_1-351.pdb --chain A --gpu 0 --thermompnn-dir /data/raw_data/computational_data/thermoMPNN_data/templates --offset 1

# --thermompnn-dir points to shared helper scripts and model weights.
# Replace the PDB, chain, GPU and starting residue (--offset) for your input.
# Both raw and postprocessed CSVs are saved in the protein folder.
