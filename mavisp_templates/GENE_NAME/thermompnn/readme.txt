#have do.sh and a trimmed model named ACCESSION_START-END.pdb in your protein folder

#run on bioinfo03 with conda available
#do.sh runs in /usr/local/envs/thermompnn automatically; no manual activation needed

#examples below use P60484_1-351.pdb; replace it with your PDB filename

#check paths without running the model
bash do.sh P60484_1-351.pdb --chain A --gpu 0 --dry-run

#check GPU availability; --gpu selects a visible GPU, not a reservation
#use --gpu 1 for the second visible GPU; omit --gpu to keep the CUDA default
nvidia-smi

#example of run (direct execution, not a queued job)
bash do.sh P60484_1-351.pdb --chain A --gpu 0

#this test do.sh uses Emma's central server checkout; no export needed
#helper modules, weights and local.yaml stay beside central inference; no YAML edits needed
#requires the updated inference with integrated postprocessing, --dry-run and --gpu

#before sharing, set the Conda initialization, environment and inference paths in do.sh for that installation

#Postprocessing for MAVISp ingestion is automatic
#both raw and _postprocessed.csv files are saved in the current folder
#existing files with the chosen output names are overwritten
#offset is read from the PDB filename; use --offset for other filenames

#for reproducibility, record the inference revision, Python environment and model checksums

