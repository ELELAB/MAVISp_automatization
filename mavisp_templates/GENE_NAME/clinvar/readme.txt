#create the input variant.csv file
#search the isoform name from Clinvar database in the HGVS field
#prepare input file
ln -s ../cancermuts/mutlist.txt mutlist.txt
cat mutlist.txt | awk '{print "ACVR1B;"$1";NP_004293"}' > tmp.csv
#add header to input file
cat <(echo "gene;variant_name;iso") tmp.csv >  variants.csv
rm tmp.csv

bash run.sh
