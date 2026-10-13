
# In order to reproduce the run you need to:
# 1. activate the cancermuts virtual environment

. /usr/local/envs/cancermuts/bin/activate

# 1. Prepare saturation mutlist

cp /data/user/shared_projects/mavisp/SUMO3-2/saturation_mutlist/saturation_mutlist.txt .
python input_csv.py saturation_mutlist.txt
mv input.csv saturation_mutlist.csv

#In case of manually curated mutation list prepare in the same way also that file

# 3. The pancancer.py script support mutations from external sources from
    manually curated lists. In order to provide mutations list from
    an external mutation list, specify the file with the following flags:
    -e in case of manually curated mutation lists (the script supports several lists so in case
       of multiple lists, specify the files tab separated)


tsp -N 1 python pancancer.py -p SUMO3 --isoform P55854-2 -r NP_065857 -e saturation_mutlist.csv

## IMPORTANT the script will generate a file with the ENST id for that isoform, you should keep it.
# it will be usefull for alphamissense module
