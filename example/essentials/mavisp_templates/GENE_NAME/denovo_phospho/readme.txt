Files required:
 -  mutlist with this format: wt position mut (ex: D57A) 
	cp ../cancermuts/mutlist_XXXXXX.txt mutlist.txt
 -  wt fasta file (you can copy from demask folder)
	cp ../demask/*.fasta

Modify the config file with your own desired files: 
    	mutlist: "mutlist.txt"
	fasta_file: "XXXXX.fasta"
	directory: "results/"

Run the scrip:
module load python/3.10/modulefile
snakemake -c4
