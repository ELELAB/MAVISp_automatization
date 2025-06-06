#have a csv with the gene name and the isoform id corresponding to the main isoform as input.
#To get the main isoform id:
#After you have searched for your protein,  click on Sequence and isoforms section, check the accession numbers of the isoform and search for the code followed by 1). The string  with the “NP” (isoform_code without the dot and the number after it) associated with that accession number is the one to put in your input file.

tsp -N 1 bash run.sh

#check on status of the run using tail -f on the corresponding tsp log file

