#!/usr/bin/env bash
set -e
source /usr/local/anaconda3/etc/profile.d/conda.sh
conda activate /usr/local/envs/thermompnn
exec python "/data/user/emmaqingjie/repos/MAVISp_automatization/mavisp_templates/GENE_NAME/thermompnn/custom_inference_standalone.py" --pdb "$@"

