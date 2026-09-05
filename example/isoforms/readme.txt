# TO RUN

Make sure that all specified needed inputs defined in ../../README.md are ready in their respective locations.

conda deactivate

module load python/3.10/modulefile

tsp -N 4 snakemake -c 4 isoforms
