# Activate the appropriate Python environment
. /usr/local/envs/py310/bin/activate

# Source AMBER tools
. /usr/local/amber-20/amber.sh

# Run the curation
python curation.py -c config_curation.yaml

