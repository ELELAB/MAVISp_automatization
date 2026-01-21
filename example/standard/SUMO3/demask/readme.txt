#have a fasta file with your full protein sequence 

wget https://www.uniprot.org/uniprot/$uniprot.fasta

#run with tsp

tsp -N 4 bash do.sh file.fasta
