
# In order to reproduce the run you need to:
# 1. activate the cancermuts virtual environment

. /usr/local/envs/cancermuts/bin/activate

# 2. Prepare input.csv
cp /data/user/shared_projects/mavisp/MAP1LC3B/clinvar_gene/MAP1LC3B_mutation_list.txt .
python input_csv.py MAP1LC3B_mutation_list.txt
mv input.csv clinvar.csv

#In case of manually curated mutation list prepare in the same way also that file

# 3. The pancancer.py script support mutations from external sources such as ClinVar or
    manually curated lists. In order to provide mutations list from Clinvar and/or from
    an external mutation list, specify the file with the following flags:
    -c in case of mutation list from clinvar
    -e in case of manually curated mutation lists (the script supports several lists so in case
       of multiple lists, specify the files tab separated)

# 4. run the script pancancer.py providing the protein name and the protein uniprotID as arguments e.g.
#     python pancancer.py -p $protein -i $uniprot_id
# notice that in some instances, the cancermuts run might fail because it is not able
# to convert the Uniprot ID to a corresponding Uniprot AC. This will show up as a crash
# related to Uniprot, with a trace like:
#   File "pancancer_clinvar.py", line 21, in <module>
#    seq = up.get_sequence(prt, upid=uniprotID)
#  File "/data/user/teo/devel/cancermuts/cancermuts/datasources.py", line 130, in get_sequence
#    this_upac = self._get_aliases(upid, ['UniProtKB_primaryAccession'])['UniProtKB_primaryAccession']
#    ...
# If this happens, the Uniprot AC can be supplied manually as
# a third argument of the script:
#     python pancancer.py -p $protein -i $uniprot_id -a $uniprot_ac
# for instance:
#     tsp -N 1 python pancancer_clinvar.py -p MAP1LC3B -i MLP3B_HUMAN -a Q9GZQ8 -c mutations_clinvar.csv -e LC3B_icope.csv LC3B_marinara.csv

# In some cases, multiple entrez IDs can be found for one protein. This will not cause an error, but
# no gene ID will be assigned and mutation retrieval from cBioPortal will be skipped.

tsp -N 1 python pancancer_clinvar.py MAP1LC3B MLP3B_HUMAN
