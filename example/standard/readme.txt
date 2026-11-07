#To run

Have all the specified needed inputs in the README.md ready in their respective locations

conda deactivate

module load python/3.10/modulefile

Run from this folder on bioinfo01 using the Snakemake command under SLURM in ../../README.md.
Snakemake submits ThermoMPNN to a GPU node through SLURM; the other steps run locally.
Processed ThermoMPNN predictions are collected under each protein's simple_mode/stability/.../thermompnn/ folder.

To queue that command with Task Spooler, prefix it with tsp -N 1 and keep all Snakemake options.
