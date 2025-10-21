import argparse
from cancermuts.datasources import ManualAnnotation
from cancermuts.datasources import UniProt
from cancermuts.datasources import cBioPortal, COSMIC
from cancermuts.datasources import MyVariant
from cancermuts.datasources import gnomAD
from cancermuts.datasources import PhosphoSite, MobiDB
from cancermuts.datasources import ggetELMPredictions
from cancermuts.table import Table

parser=argparse.ArgumentParser(description='Pancancer: run Cancermuts"\
                               "software')
parser.add_argument("-p", "--prt", dest="prt", help="hugo name of your protein")
parser.add_argument("-i", "--uniprotID", dest="uniprotID", help="uniprot_id of your protein (the one ending with _HUMAN)")
parser.add_argument("-a", "--uniprotAC", dest="uniprotAC", default= None, help="uniprot_ac of your protein")
parser.add_argument("-c", "--clinvar", dest="clinvar", default= None, help="csv mutation file from clinvar database")
parser.add_argument("-e", "--external_mutations", dest="external_mutations", nargs='+', default= None, help="csv file containing the external mutations to study")

args=parser.parse_args()
# create the corresponding uniprot object
up = UniProt()

# get the sequence for the protein
seq = up.get_sequence(args.prt, upid=args.uniprotID, upac=args.uniprotAC)
print(seq.sequence)

# add mutations from cBioPortal
cb = cBioPortal()
try:
    cb.add_mutations(seq, metadata=['cancer_type', 'cancer_study', 'genomic_mutations'])
except TypeError:
    print("WARNING: Skipping cBioPortal due to missing Entrez ID.")


# add mutations from COSMIC
cosmic = COSMIC(targeted_database_file='/data/databases/cosmic-v102/Cosmic_CompleteTargetedScreensMutant_v102_GRCh38.tsv',
				screen_mutant_database_file='/data/databases/cosmic-v102/Cosmic_GenomeScreensMutant_v102_GRCh38.tsv',
				classification_database_file='/data/databases/cosmic-v102/Cosmic_Classification_v102_GRCh38.tsv',
				database_encoding='latin1', lazy_load_db=True,
                )
cosmic.add_mutations(seq, genome_assembly_version='GRCh38', metadata=['genomic_coordinates', 'genomic_mutations',
                                                'cancer_site', 'cancer_histology'])
# add mutations from other sources
if args.clinvar:
    ma = ManualAnnotation(args.clinvar)

    # add mutations to the seq object
    ma.add_mutations(seq, metadata=['genomic_mutations'])

    # add PTM annotations to the sequence object
    ma.add_position_properties(seq)

    # add structure or linear motif annotation to the sequence object
    ma.add_sequence_properties(seq)


if args.external_mutations:
    for external_list in args.external_mutations:
        ma = ManualAnnotation(external_list)

        # add mutations to the seq object
        ma.add_mutations(seq)

        # add PTM annotations to the sequence object
        ma.add_position_properties(seq)

        # add structure or linear motif annotation to the sequence object
        ma.add_sequence_properties(seq)


# add annotations from MyVariant (REVEL)
mv = MyVariant()
mv.add_metadata(seq)

# add annotations from gnomAD
gnomad = gnomAD(version='2.1')
gnomad.add_metadata(seq, md_type=['gnomad_exome_allele_frequency',
	                              'gnomad_genome_allele_frequency'])

# add annotations from PhosphoSite
ps = PhosphoSite('/data/databases/phosphosite/')
ps.add_position_properties(seq)

# add annotations from MobiDB
#mdb = MobiDB()
#mdb.add_position_properties(seq)

# add annotations from ELM
elm = ggetELMPredictions()
elm.add_sequence_properties(seq,
			    exclude_elm_classes="MOD_.")


# save table
tbl = Table()
df = tbl.to_dataframe(seq)
df.to_csv("metatable_pancancer_"+args.prt+".csv")
tbl.plot_metatable(df, fname="plot_metatable_pancancer_"+args.prt+".pdf", section_size=50)
