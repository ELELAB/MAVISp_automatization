pdb=$1
list=$2

. /usr/local/envs/mutatex/bin/activate

#label file from PDB
pdb2labels -p $pdb

set +e
ddg2summary -p $pdb -d final_averages -l mutation_list.txt -L $list &> log
exit_code=$?
echo $exit_code

if [[ $exit_code -eq 1 ]]; then
    line=$(grep 'no usable mutations were found in mutations file; exiting...' log)
    if [[ ! -z $line ]]; then
        touch summary.txt
    else
        exit 1
    fi
fi
set -e
rm log
mv summary.txt summary_stability.txt

