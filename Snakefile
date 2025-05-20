import pandas as pd
import os
import re
import shutil
import time
import configparser
import yaml
import urllib
from Bio.Data import IUPACData

def mutation_converter(x):
    return f"p.{IUPACData.protein_letters_1to3.get(x[0])}{x[1:-1]}\
                        {IUPACData.protein_letters_1to3.get(x[-1])}"


configfile: "config.yaml"

######################### Script and readme paths ###########################


modules = config['modules']

# the script and the readme files for each module are added in the config
# dictionary

#------------------------------ Interactome --------------------------------#

mentha_script = f"mavisp_templates/GENE_NAME/interactome/"\
                f"mentha2pdb/mentha2pdb.py"
mentha_readme = f"mavisp_templates/GENE_NAME/interactome/"\
                f"mentha2pdb/readme.txt"
mentha_bash =   f"mavisp_templates/GENE_NAME/interactome/"\
                f"mentha2pdb/do.sh"
hpc_atlas_readme = f"mavisp_templates/GENE_NAME/interactome/"\
                   f"hpc_atlas/readme.txt"

modules['interactome']['mentha2pdb'].update({"script":mentha_script,
                                             "readme":mentha_readme,
                                             "bash":mentha_bash})
modules['interactome']['hpc_atlas'].update({"readme":hpc_atlas_readme})

#--------------------------- Structure Selection ---------------------------#

alphafold_script = f"mavisp_templates/GENE_NAME/"\
                   f"structure_selection/alphafold_db/get_alphafolddb_data.py"
alphafold_readme = f"mavisp_templates/GENE_NAME/"\
                   f"structure_selection/alphafold_db/readme.txt"

filter_pdb_script = f"mavisp_templates/GENE_NAME/"\
                    f"structure_selection/trimmed_models/filter_pdb.py"
filter_pdb_readme = f"mavisp_templates/GENE_NAME/"\
                    f"structure_selection/trimmed_models/readme.txt"

pdbminer_readme = f"/mavisp_templates/GENE_NAME/"\
                  f"structure_selection/pdbminer/readme.txt"

pdbminer_complexes_readme = f"mavisp_templates/GENE_NAME/"\
                            f"structure_selection/pdbminer_complexes/readme.txt"
pdbminer_complexes_script = f"mavisp_templates/GENE_NAME/"\
                            f"structure_selection/pdbminer_complexes/"\
                            f"find_PDBminer_complexes.py"

procheck_readme = f"/mavisp_templates/GENE_NAME/"\
                  f"structure_selection/procheck/readme.txt"


modules['structure_selection']["alphafold"].update({"script":alphafold_script,
                                                   "readme":alphafold_readme})

modules['structure_selection'].update({"trimmed_models":\
                                      {"script":filter_pdb_script,
                                       "readme":filter_pdb_readme}})

modules['structure_selection']['pdbminer'].update({"readme":pdbminer_readme})

modules['structure_selection'].update({"pdbminer_complexes":\
                                      {"script": pdbminer_complexes_script,
                                       "readme": pdbminer_complexes_readme}})

modules['structure_selection']['procheck'].update({"readme":procheck_readme})


#---------------------------- Aggregation step -----------------------------#


cancermuts_readme = "mavisp_templates/GENE_NAME/cancermuts_metatable/readme.txt"

cancermuts_input_script = "mavisp_templates/GENE_NAME/cancermuts_metatable/input_csv.py"

cancermuts_pancancer_script = "mavisp_templates/GENE_NAME/cancermuts_metatable/pancancer.py"

cancermuts_plot_script = "mavisp_templates/GENE_NAME/cancermuts_metatable/plot.py"


modules['mutations_aggregation']['cancermuts'].update({"readme" : cancermuts_readme,
                                                       "script_inputs" : cancermuts_input_script,
                                                       "script_pancancer" : cancermuts_pancancer_script,
                                                       "script_plot" : cancermuts_plot_script})

#--------------------------------- ClinVar ---------------------------------#

clinvar_gene_script = f"mavisp_templates/GENE_NAME/"\
                      f"clinvar_gene/clinvar.py"
clinvar_gene_readme = f"mavisp_templates/GENE_NAME/"\
                      f"clinvar_gene/readme.txt"
clinvar_gene_bash =   f"mavisp_templates/GENE_NAME/"\
                      f"clinvar_gene/run.sh"

modules.update({"ClinVar_database":{"clinvar_gene":\
                                         {"script" : clinvar_gene_script,
                                          "readme" : clinvar_gene_readme,
                                          "bash"   : clinvar_gene_bash}}})

column_names = ['name',
                'site',
                'type',
                'function',
                'reference',
                'genomic_mutations']

HGVSp_re = re.compile('\((p.[A-Z][a-z]{2}[0-9]+[A-Z][a-z]{2})\)')

GRCh_to_hg = {'GRCh37' : 'hg19',
              'GRCh38' : 'hg38'}

refseq_to_UCSC = {'hg38': {
                        "NC_000001.11" : "1",
                        "NC_000002.12" : "2",
                        "NC_000003.12" : "3",
                        "NC_000004.12" : "4",
                        "NC_000005.10" : "5",
                        "NC_000006.12" : "6",
                        "NC_000007.14" : "7",
                        "NC_000008.11" : "8",
                        "NC_000009.12" : "9",
                        "NC_000010.11" : "10",
                        "NC_000011.10" : "11",
                        "NC_000012.12" : "12",
                        "NC_000013.11" : "13",
                        "NC_000014.9"  : "14",
                        "NC_000015.10" : "15",
                        "NC_000016.10" : "16",
                        "NC_000017.11" : "17",
                        "NC_000018.10" : "18",
                        "NC_000019.10" : "19",
                        "NC_000020.11" : "20",
                        "NC_000021.9"  : "21",
                        "NC_000022.11" : "22",
                        "NC_000023.11" : "X",
                        "NC_000024.10" : "Y"
                }, 'hg19' : {
                        "NC_000001.10" : "1",
                        "NC_000002.11" : "2",
                        "NC_000003.11" : "3",
                        "NC_000004.11" : "4",
                        "NC_000005.9"  : "5",
                        "NC_000006.11" : "6",
                        "NC_000007.13" : "7",
                        "NC_000008.10" : "8",
                        "NC_000009.11" : "9",
                        "NC_000010.10" : "10",
                        "NC_000011.9"  : "11",
                        "NC_000012.11" : "12",
                        "NC_000013.10" : "13",
                        "NC_000014.8"  : "14",
                        "NC_000015.9"  : "15",
                        "NC_000016.9"  : "16",
                        "NC_000017.10" : "17",
                        "NC_000018.9"  : "18",
                        "NC_000019.9"  : "19",
                        "NC_000020.10" : "20",
                        "NC_000021.8"  : "21",
                        "NC_000022.10" : "22",
                        "NC_000023.10" : "X",
                        "NC_000024.9"  : "Y"
                }
}

def variant_name_to_HGVSp(row):
    hgvsps = re.findall(HGVSp_re, row['variant_name'])
    if len(hgvsps) != 1:
        print(f"WARNING: skipping {row.variant_id}, {row.variant_name}")
        return 'unexpected'
    if hgvsps[0].endswith('Ter'):
        print(f"WARNING: skipping {row.variant_id}, {row.variant_name}")
        return 'unexpected'

    return hgvsps[0]

def HGVSg_to_cancermuts(row):
    hgvsgs = row['genomic_annotation'].split()
    out = []

    for hgvsg in hgvsgs:
        build, variant = hgvsg.split(',')
        ref_seq, coords = variant.split(':')

        try:
            build = GRCh_to_hg[build]
        except KeyError:
            print(f"WARNING: unexpected genome build for {row['variant_id']} ({ref_seq}). Genomic mutation won't be annotated. ")
            continue

        ref_seq = refseq_to_UCSC[build][ref_seq]

        out.append(f"{build},{ref_seq}:{coords}")

    return " ".join(out)

#------------------------------- Mutlist -----------------------------------#

mutlist_script = f"mavisp_templates/GENE_NAME/"\
                 f"cancermuts/get_mutlists.py"
mutlist_readme = f"mavisp_templates/GENE_NAME/"\
                 f"cancermuts/readme.txt"

modules.update({"mutlist_generation":{"script":mutlist_script,
                                      "readme":mutlist_readme}})

saturation_mutlist_sh = "mavisp_templates/GENE_NAME/"\
                        "saturation_list/saturation_mutlist.py"
saturation_mutlist_readme = "mavisp_templates/GENE_NAME/"\
                            "saturation_list/readme.txt"
saturation_mutlist_py = "mavisp_templates/GENE_NAME/"\
                        "saturation_list/do.sh"

modules.update({"saturation_mutlist_generation":{"py":saturation_mutlist_py,
                                                 "sh":saturation_mutlist_sh,
                                                 "readme":mutlist_readme}})


#--------------------------- Protein annotations ---------------------------#

domains_script = f"mavisp_templates/GENE_NAME/"\
                 f"structure_selection/domain_annotations/get_domains.py"
domains_readme = f"mavisp_templates/GENE_NAME/"\
                 f"structure_selection/domain_annotations/readme.txt"

modules.update({"domain_annotations":{"script":domains_script,
                                      "readme":domains_readme}})

netphos_readme = f"mavisp_templates/GENE_NAME/"\
                 f"netphos/readme.txt"

modules.update({"netphos":netphos_readme})

ptm_mutatex_script  = "mavisp_templates/GENE_NAME/ptm/mutatex/run_ddgs.sh"
ptm_mutatex_readme  = "mavisp_templates/GENE_NAME/ptm/mutatex/readme.txt"
ptm_mutatex_mutlist = "mavisp_templates/GENE_NAME/ptm/mutatex/mutation_list.txt"
ptm_naccess_readme  = "mavisp_templates/GENE_NAME/ptm/naccess/readme.txt"

modules['mutations_aggregation']['ptm'] = {'mutatex' : {'script'  : ptm_mutatex_script,
                                                                 'readme'  : ptm_mutatex_readme,
                                                                 'mutlist' : ptm_mutatex_mutlist},
                                          'naccess' : {'readme'  : ptm_naccess_readme}}

#--------------------------- Mutations Classifiers -------------------------#

demask_readme = f"mavisp_templates/GENE_NAME/demask/readme.txt"
demask_script = f"mavisp_templates/GENE_NAME/demask/do.sh"


modules['mutations_classifier']["demask"].update({"readme":demask_readme,
                                                  "script":demask_script})


alphamissense_script = f"mavisp_templates/GENE_NAME/alphamissense/do.sh"
alphamissense_readme = f"mavisp_templates/GENE_NAME/alphamissense/readme.txt"

modules['mutations_classifier']\
       ["alphamissense"].update({"readme":alphamissense_readme,
                                  "script":alphamissense_script})

#------------------------------ Calculations -------------------------------#

rasp_script = f"/data/raw_data/computational_data/rasp_data/"\
              f"mavisp_templates/saturation/AF2_XX-YY/model_X/run_analysis.sh"

rasp_readme = f"/data/raw_data/computational_data/rasp_data/"\
              f"mavisp_templates/saturation/AF2_XX-YY/model_X/readme.txt"

modules["rasp"].update({"script":rasp_script,
                        "readme":rasp_readme})

rosetta_relax_path = f"/data/raw_data/computational_data/rosetta_data/"\
                     f"mavisp_templates/stability/ref2015_cartesian2020/"\

rosetta_relax_yaml = f"{rosetta_relax_path}relax2020_ref2015.yaml"
rosetta_relax_cartesian_yaml = f"{rosetta_relax_path}\
                                 cartesian2020_ref2015.yaml"
rosetta_relax_cartdgg_yaml = f"{rosetta_relax_path}cartddg2020_ref2015.yaml"
rosettampi_yaml = f"{rosetta_relax_path}rosettampi.yaml"
rosetta_relax_bash_script = f"{rosetta_relax_path}aggregate2mutatex.sh"
rosetta_relax_aggregate_yaml = f"{rosetta_relax_path}aggregate.yaml"
rosetta_relax_readme = f"{rosetta_relax_path}readme.txt"

modules['rosetta_relax'].update({"relax_yaml":rosetta_relax_yaml,
                                 "cart_yaml":rosetta_relax_cartesian_yaml,
                                 "cartdgg_yaml":rosetta_relax_cartdgg_yaml,
                                 "script":rosetta_relax_bash_script,
                                 "agg_yaml":rosetta_relax_aggregate_yaml,
                                 "rosettampi_yaml":rosettampi_yaml,
                                 "readme":rosetta_relax_readme})


#--------------------------------- Allosigma -------------------------------#

allosigma_aminoacid = f"mavisp_templates/GENE_NAME/long_range/"\
                      f"allosigma2/2.allosigma_classify/aminoacids.dat"
allosigma1_readme   = f"/mavisp_templates/GENE_NAME/long_range/"\
                      f"allosigma2/1.allosteric_signalling_map/readme.txt"

allosigma_classify  = f"mavisp_templates/GENE_NAME/long_range/"\
                      f"allosigma2/2.allosigma_classify/allosigma-classify"
allosigma2_readme   = f"mavisp_templates/GENE_NAME/long_range/"\
                      f"allosigma2/2.allosigma_classify/readme.txt"

allosigma_heatmap   = f"mavisp_templates/GENE_NAME/long_range/"\
                      f"allosigma2/3.allosigma_heatmap/all/allosigma-heatmap"
allosigma3_readme   = f"mavisp_templates/GENE_NAME/long_range/"\
                      f"allosigma2/3.allosigma_heatmap/all/readme.txt"

allosigma_filtering = f"mavisp_templates/GENE_NAME/long_range/"\
                      f"allosigma2/4.allosigma.filtering/allosigma-filtering"
allosigma4_readme   = f"mavisp_templates/GENE_NAME/long_range/"\
                      f"allosigma2/4.allosigma.filtering/readme.txt"

modules.update({"allosigma":{"allosigma1":{"aminoacids":allosigma_aminoacid,
                                           "readme":allosigma1_readme},
                             "allosigma2":{"script":allosigma_classify,
                                           "readme":allosigma2_readme},
                             "allosigma3":{"script":allosigma_heatmap,
                                           "readme":allosigma3_readme},
                             "allosigma4":{"script":allosigma_filtering,
                                           "readme":allosigma4_readme}}})

df = pd.read_csv(config['input']['path'])

# Rasp
rasp_path=modules['rasp']['output_path_folder']
df['output_path_folder'] = rasp_path

#Rosetta
rosetta_path=modules['rosetta_relax']['rosetta_folder']
df['output_rosetta_folder'] = rosetta_path
df["trimmed"] = df["trimmed"].str.split("_")
df_exploded = df.explode("trimmed")
#-------------------------------- Denovo phospho ------------------------------#


snakefile = "mavisp_templates/GENE_NAME/denovo_phospho/snakefile"
readme = "mavisp_templates/GENE_NAME/denovo_phospho/readme.txt"

modules.update({"denovo_phospho":{"snakefile":snakefile,
                                  "readme":readme}})

#---------------------------------- Efoldmine ---------------------------------#

efoldmine_readme = f"mavisp_templates/GENE_NAME/efoldmine/readme.md"
modules['efoldmine'].update({"readme": efoldmine_readme})


##############################################################################

                                   #Snakefile#

##############################################################################

rule all:
    input:
        expand("{hugo_name}/structure_selection/"\
                "original_model/",
                zip,
                hugo_name = df['protein'].str.upper(),
                structure_source = df['structure_source']),

        expand("{hugo_name}/structure_selection/"\
               "domain_annotations/"\
               "domains_mutlist.csv",
                zip,
                hugo_name = df['protein'].str.upper(),
                structure_source = df['structure_source']),

        expand("{hugo_name}/netphos/"\
               "netphos.out",
              hugo_name = df['protein'].str.upper()),

        expand("{hugo_name}/interactome/"\
               "hpc_atlas/{hugo_name}.out",
              hugo_name = df['protein'].str.upper()),

        expand("{hugo_name}/"\
               "structure_selection/"\
               "pdbminer/results/"\
               "{uniprot_ac}/"\
               "{uniprot_ac}_all.csv",
              zip,
              hugo_name = df['protein'].str.upper(),
              uniprot_ac = df['uniprot_ac'].str.upper()),

        expand("{hugo_name}/demask/"\
               "myquery_predictions.txt",
               hugo_name = df['protein'].str.upper()),

        expand("{hugo_name}/alphamissense/"\
               "am.tsv.gz",
               hugo_name = df['protein'].str.upper()),

        expand("{hugo_name}/"\
               "structure_selection/"\
               "trimmed_model/",
               zip,
               hugo_name = df['protein'].str.upper(),
               structure_source = df['structure_source']),

        expand("{path}/"\
               "{research_field}/"\
               "{hugo_name}/free/{structure_source}_{resrange}/{model}_model/",
               zip,
               resrange = df_exploded['trimmed'],
               hugo_name = df_exploded['protein'].str.lower(),
               path = df_exploded['output_path_folder'],
               research_field = df_exploded['research_field'],
               structure_source = df_exploded['structure_source'],
               model = df_exploded['model']),

        expand("{hugo_name}/interactome/"\
               "mentha2pdb/"\
               "{uniprot_ac}.csv",
               zip,
               hugo_name = df['protein'].str.upper(),
               uniprot_ac = df['uniprot_ac'].str.upper()),

        expand("{path}/"\
               "{research_field}/"\
               "{hugo_name}/free/"\
               "{structure_source}_{resrange}/"\
               "{model}_model/ref2015_cartesian2020/relax/relax_{uniprot_ac}_{resrange}_0001.pdb",
               zip,
               uniprot_ac = df_exploded['uniprot_ac'].str.upper(),
               hugo_name = df_exploded['protein'],
               path = df_exploded['output_rosetta_folder'],
               resrange = df_exploded['trimmed'],
               research_field = df_exploded['research_field'],
               structure_source = df['structure_source'],
               model = df['model']),

        expand("{hugo_name}/ptm/{structure_source}_{resrange}/mutatex/summary_stability.txt",
            zip,
            hugo_name = df_exploded['protein'],
            resrange = df_exploded['trimmed'],
            structure_source = df_exploded['structure_source']),

        expand("{hugo_name}/ptm/{structure_source}_{resrange}/naccess/{uniprot_ac}_trimmed_model0_checked.rsa",
            zip,
            hugo_name = df_exploded['protein'],
            resrange = df_exploded['trimmed'],
            uniprot_ac = df_exploded['uniprot_ac'].str.upper(),
            structure_source = df_exploded['structure_source']),

        expand("{hugo_name}/efoldmine/{uniprot_ac}.tabular",
            zip,
            hugo_name=df['protein'].str.upper(),
            uniprot_ac=df['uniprot_ac'].str.upper()),

        expand("{hugo_name}/structure_selection/pdbminer_complexes/{uniprot_ac}_filtered.csv",
            zip,
            hugo_name = df['protein'].str.upper(),
            uniprot_ac = df['uniprot_ac'].str.upper()),

        expand("{hugo_name}/structure_selection/procheck/",
               zip,
               hugo_name=df['protein'].str.upper()),

        expand(["{hugo_name}/metadata/metadata.yaml",
                "{hugo_name}/metadata/importing.yaml"],
                hugo_name= df_exploded['protein'].str.upper()),

        expand("{hugo_name}/simple_mode/collection_{research_field}_{structure_source}_{resrange}_{uniprot_ac}_{model}.done",
                zip,
                hugo_name = df_exploded['protein'].str.upper(),
                research_field = df_exploded['research_field'],
                structure_source = df_exploded['structure_source'],
                resrange = df_exploded['trimmed'],
                uniprot_ac = df_exploded['uniprot_ac'].str.upper(),
                model = df_exploded['model'])

###################### Structure selection and trimming ######################

rule structure_selection:
    output:
        directory("{hugo_name}/structure_selection/original_model")

    run:

        pdb = df.loc[
            (df['protein'] == wildcards.hugo_name),
            'input_pdb'].iloc[0]

        if pd.isna(pdb):

            uniprot_ac = df.loc[df["protein"] == wildcards.hugo_name,
                                                "uniprot_ac"].iloc[0]
            dssp_exec=modules['structure_selection']['alphafold']['dssp_exec']

            # Alphafold: Create the config.yaml file for the analysis

            if not os.path.exists(str(output)):
                os.makedirs(str(output))
            data = {
                "dssp_exec": dssp_exec,
                "plddt_cutoff": 70,
                "uniprot_ids": {
                    uniprot_ac: {
                        "dir_name": wildcards.hugo_name.lower(),
                        "version": "latest"
                    }
                }
            }

            # Specify the output YAML file path
            file_path = "output.yaml"

            with open(f"{output}/config_alphafolddb.yaml", "w") as f:
                yaml.dump(data, f)

            # run alphafold module

            readme=modules['structure_selection']['alphafold']['readme']
            script=modules['structure_selection']['alphafold']['script']
            shell("cd {output} && "
                "cp ../../../{readme} . && "
                "cp ../../../{script} . && "
                "python get_alphafolddb_data.py -c config_alphafolddb.yaml")
        else:
            structure_folder = f'{output}/{wildcards.hugo_name.lower()}'
            shell("mkdir -p {structure_folder}")
            shell("cp {pdb} {structure_folder}")


rule trim_model:
    input:
        "{hugo_name}/structure_selection/original_model"
    output:
        directory("{hugo_name}/structure_selection/trimmed_model")
    run:

        # list with all the residue ranges specified in the input file
        trimmed_pdb_list=df.loc[df['protein'] == wildcards.hugo_name,\
                                                'trimmed'].iloc[0]

        uniprot_ac  = df.loc[df['protein'] == wildcards.hugo_name,\
                                                'uniprot_ac'].iloc[0]
        input_files = f'{input}/{wildcards.hugo_name.lower()}' # path of the pdb files

        # list containing all the files in the folder
        files = os.listdir(input_files)

        # keep only the pdb file
        file_pdb_list = [f for f in files if f.endswith(".pdb")]

        # for every range specified in the input file create
        # the corresponding trimmed pdb file

        script=modules['structure_selection']['trimmed_models']['script']
        readme=modules['structure_selection']['trimmed_models']['readme']
        for i in trimmed_pdb_list:
            number_list = i.split("-")
            start = number_list[0]
            end = number_list[1]
            trimmed_pdb = f'{uniprot_ac}_{start}-{end}.pdb'
            shell("mkdir -p {output} &&"\
                  " python {script}"\
                      " {input_files}/{file_pdb_list[0]}"
                      " A"\
                      " {start}"\
                      " {end}"\
                      " {wildcards.hugo_name}/structure_selection/"\
                      "trimmed_model/{trimmed_pdb} && "\
                      "cp {script} {readme} {output} ")

rule pdbminer:
    output:
        "{hugo_name}/structure_selection/pdbminer/results/"\
        "{uniprot_ac}/{uniprot_ac}_all.csv"
    shell:
        '''
        readme={modules[structure_selection][pdbminer][readme]}
        set +eu; {config[modules][structure_selection][pdbminer][pdbminer_env]}; set -eu
        mkdir -p {wildcards.hugo_name}/structure_selection/pdbminer/
        cd {wildcards.hugo_name}/structure_selection/pdbminer/
        cp ../../../$readme .
        PDBminer -g {wildcards.hugo_name} -u {wildcards.uniprot_ac} -f csv
        '''

rule pdbminer_complexes:
    input:
        "{hugo_name}/structure_selection/pdbminer/results/{uniprot_ac}/{uniprot_ac}_all.csv"
    output:
        "{hugo_name}/structure_selection/pdbminer_complexes/{uniprot_ac}_filtered.csv",
    shell:
        '''
        mkdir -p "{wildcards.hugo_name}/structure_selection/pdbminer_complexes/"
        cd {wildcards.hugo_name}/structure_selection/pdbminer_complexes/
        readme={modules[structure_selection][pdbminer_complexes][readme]}
        script={modules[structure_selection][pdbminer_complexes][script]}
        cp  ../../../${{readme}} .
        cp  ../../../${{script}} .
        python find_PDBminer_complexes.py\
               -i ../pdbminer/results/{wildcards.uniprot_ac}/{wildcards.uniprot_ac}_all.csv\
               -o {wildcards.uniprot_ac}_filtered.csv\
               --binding_interface -d 10
        if ls *.pdb 1> /dev/null 2>&1; then
            mkdir -p {wildcards.uniprot_ac}_pdb_complexes/
            mv *.pdb {wildcards.uniprot_ac}_pdb_complexes/
        fi
        '''

rule procheck:
    input:
        trimmed_pdb_dir="{hugo_name}/structure_selection/trimmed_model/"  # The directory containing the PDB files
    output:
        directory("{hugo_name}/structure_selection/procheck/")  # The directory for the output .sum files
    params:
        chain="A",
        resolution="2.0",
    run:
        # Get path to the directory containing PDB files
        trimmed_pdb_dir = os.path.abspath(input.trimmed_pdb_dir)
        procheck_env=modules['structure_selection']['procheck']
        # List all PDB files in the trimmed_model directory
        pdb_files = [f for f in os.listdir(trimmed_pdb_dir) if f.endswith(".pdb")]

        # Ensure output directory exists
        output_dir = os.path.abspath(output[0])
        os.makedirs(output_dir, exist_ok=True)
        # Iterate over each PDB file
        for pdb_file in pdb_files:
            pdb_input_path = os.path.join(trimmed_pdb_dir, pdb_file)
            pdb_filename = os.path.splitext(pdb_file)[0]  # Remove the .pdb extension to construct the output file name

            # Output .sum file path
            sum_output_path = os.path.join(output_dir, f"{pdb_filename}.sum")

            # Run PROCHECK for each PDB file
            shell(
        """
                set +eu; {config[modules][structure_selection][procheck][procheck_env]}; set -eu
                readme={modules[structure_selection][procheck][readme]}
                cd {output_dir}
                cp ../../../$readme .
                procheck.scr {pdb_input_path} {params.chain} {params.resolution}
            """)

            # Ensure the .sum file is generated (assuming the procheck.scr command generates the .sum file in the current directory)
            if not os.path.exists(sum_output_path):
                raise FileNotFoundError(f"PROCHECK failed to generate {sum_output_path}.")

############################### Interactome #################################

rule mentha2pdb:
    output:
        "{hugo_name}/interactome/mentha2pdb/{uniprot_ac}.csv"
    shell:
       '''
       mkdir -p {wildcards.hugo_name}/interactome/mentha2pdb/
       readme={modules[interactome][mentha2pdb][readme]}
       script={modules[interactome][mentha2pdb][script]}
       bash_script={modules[interactome][mentha2pdb][bash]}
       Huri={modules[interactome][mentha2pdb][AF_Huri_path]}
       Humap={modules[interactome][mentha2pdb][AF_HuMAP_path]}
       Huri_Humap={modules[interactome][mentha2pdb][AF_Huri_HuMAP_path]}
       cd {wildcards.hugo_name}/interactome/mentha2pdb/
       cp  ../../../${{readme}} .
       cp  ../../../${{script}} .
       cp  ../../../${{bash_script}} .
       ln -snf {modules[interactome][mentha2pdb][database_path]}/{modules[interactome][mentha2pdb][database_date]}
       echo {wildcards.uniprot_ac} > target_uniprot_ac.txt
       python mentha2pdb.py\
              -i {modules[interactome][mentha2pdb][database_date]}\
              -t target_uniprot_ac.txt\
              -s 0.2\
              -o {wildcards.uniprot_ac}.csv\
              -p\
              -a\
              -extra ${{Huri}} ${{Humap}}\
              -af ${{Huri_Humap}}\
              -ec 0.2
       '''

rule hpc_atlas:
    output:
        "{hugo_name}/interactome/hpc_atlas/{hugo_name}.out"
    shell:
        '''
        mkdir -p "{wildcards.hugo_name}/interactome/hpc_atlas"
        cd "{wildcards.hugo_name}/interactome/hpc_atlas"
        cp ../../../{modules[interactome][hpc_atlas][readme]} .
        ln -snf {modules[interactome][hpc_atlas][database_path]}
        grep -w "{wildcards.hugo_name}" HPC-Atlas_gene.txt > \
                "{wildcards.hugo_name}.out" ||
        echo "No entries in hpc_atlas database for "\
             "{wildcards.hugo_name} entry" > myoutput
        '''


###############################Efoldmine####################################

rule efoldmine:
    output:
        "{hugo_name}/efoldmine/{uniprot_ac}.tabular"
    shell:
        '''
        mkdir -p "{wildcards.hugo_name}/efoldmine/"
        readme={modules[efoldmine][readme]}
        cd "{wildcards.hugo_name}/efoldmine/"
        cp ../../${{readme}} .
        wget https://rest.uniprot.org/uniprotkb/{wildcards.uniprot_ac}.fasta
        {modules[efoldmine][environment]} &&
        b2bTools -i {wildcards.uniprot_ac}.fasta -t $(basename {output}) -o $(basename {output}) --efoldmine
        '''


################## mutations retrieval and aggregation ######################

rule clinvar:
    output:
        "{hugo_name}/clinvar_gene/genes_output.csv"
    run:
        # create the input for the clinvar.py script
        clinvar_input = df[['protein', 'ref_seq']]
        clinvar_input.rename(columns={'protein' : 'gene',
                                      'ref_seq' : 'isoform'}, inplace=True)
        clinvar_input = clinvar_input[(clinvar_input['gene'] == \
                                       wildcards.hugo_name)]
        clinvar_input.to_csv(os.path.join(wildcards.hugo_name,
                                          "clinvar_gene",
                                          "gene.csv"),
                                           sep=";",
                                           index=False)
        script=modules['ClinVar_database']['clinvar_gene']['script']
        readme=modules['ClinVar_database']['clinvar_gene']['readme']
        bash=modules['ClinVar_database']['clinvar_gene']['bash']
        shell(
            f"cd {os.path.join(wildcards.hugo_name, 'clinvar_gene')} && "
            f"cp ../../{readme} . && "
            f"cp ../../{script} . && "
            f"cp ../../{bash} . && "
            f"bash run.sh")

rule saturation_list:
    output:
        "{hugo_name}/saturation_mutlist/saturation_mutlist.txt"
    params:
        ranges=lambda wcs: df.loc[df['protein'] == wcs.hugo_name, 'trimmed'].tolist(),
        uniprot_ac=lambda wcs: df.loc[df['protein'] == wcs.hugo_name, 'uniprot_ac'].tolist()[0]
    run:
        module_dir = os.path.join(wildcards.hugo_name, 'saturation_mutlist')

        shutil.copy(modules['saturation_mutlist_generation']['py'], module_dir)
        shutil.copy(modules['saturation_mutlist_generation']['sh'], module_dir)
        shutil.copy(modules['saturation_mutlist_generation']['readme'], module_dir)
        ranges = params.ranges[0]
        for r in ranges:
            shell(f"cd {module_dir} && "
                  f"bash do.sh {params.uniprot_ac} {r}")
        shell(f"cd {module_dir} && "
              f"cat saturation_mutlist_* > saturation_mutlist.txt")

rule cancermuts:
    input:
        clinvar_output="{hugo_name}/clinvar_gene/genes_output.csv",
        saturation_mutlist=lambda wcs: f"{wcs.hugo_name}/saturation_mutlist/saturation_mutlist.txt" if modules['mutations_aggregation']['cancermuts']['saturation'] else [],
    output:
        f"{modules['mutations_aggregation']['cancermuts']['folder_name']}"+
        "{hugo_name}/metatable_pancancer_{hugo_name}.csv"
    resources:
        elm_connections=modules['mutations_aggregation']\
                            ['cancermuts']\
                            ['ELM_connections_per_run']
    retries:
        3
    run:
        input_files=[]

    # check the presence of external mutation lists and
    # copy them inside the protein folder

        for f in os.listdir('.'):
            filepath = os.path.join(".", f)
            if os.path.isfile(filepath):
                if f.startswith(wildcards.hugo_name) or\
                   f.startswith(wildcards.hugo_name.lower()) and\
                   f.endswith(".txt"):
                    shell("mkdir -p {wildcards.hugo_name}/"\
                                    "external_mutation_lists/")
                    shutil.copy(f,f'{wildcards.hugo_name}/'\
                                  f'external_mutation_lists/{f}')
                    input_files.append(f)

        path=modules['mutations_aggregation']['cancermuts']['folder_name']+\
                    str(wildcards.hugo_name)

        if input.saturation_mutlist:
            shutil.copy(input.saturation_mutlist, path)


        ##################### mutation list from clinvar ####################

        # read the mutation list in output from the clinvar.py script
        df_clinvar=pd.read_csv(input.clinvar_output,sep=";")

        # if there are mutations reported in clinvar create the corresponding
        # input file for the cancermuts run

        if not df_clinvar.empty:
            df_clinvar['name'] = ''

            df_clinvar['site'] = df_clinvar.apply(variant_name_to_HGVSp, axis=1)
            df_clinvar = df_clinvar[ df_clinvar['site'] != 'unexpected' ]

            df_clinvar['type'] = 'mutation'
            df_clinvar['function'] = ''
            df_clinvar['reference'] = ''

            df_clinvar['genomic_mutations'] = df_clinvar.apply(HGVSg_to_cancermuts, axis=1)

            df_clinvar.to_csv(f'{path}/clinvar.csv',
                                sep=";",
                                index=False)

        else:
            print("No mutations reported in Clinvar, Cancermuts will be run without the clinvar.csv")

         ################ mutation list from external sources ###############

        # write the cancermuts input file in which the mutations
        # have been converted in one letter code
        # check for external mutation lists
        if input_files:
            for i in input_files:
                df_external = pd.read_csv(i, header=None)
                df_external['mutations'] = df_external[0].apply(lambda x: \
                f"p.{IUPACData.protein_letters_1to3.get(x[0])}"\
                f"{x[1:-1]}{IUPACData.protein_letters_1to3.get(x[-1])}")
                df_external=df_external.rename(columns={"mutations":"site",
                                                        0:"name"})
                df_external["type"]="mutation"
                df_external["function"]=""
                df_external["reference"]=""
                output_file=i.split(".")[0]
                df_external.to_csv(f'{path}/{output_file}_input.csv',
                                   sep=";",
                                   index=False)

        ###################  mutation lists from COSMIC #################

         # Get the mutation related to a specific protein from COSMIC
         # database
        cosmic_database_path = modules["mutations_aggregation"]\
                                      ["cancermuts"]\
                                      ["Cosmic_database_path"]
        shell("cd {path}/ &&"\
                " head -n 1 {cosmic_database_path} > header.txt &&"\
                " (grep -w {wildcards.hugo_name} {cosmic_database_path} > content.txt || true) &&"\
                " cat header.txt content.txt > \
                 COSMIC_ini.csv && rm header.txt content.txt")

        # copy the metatable plotting script in the folder

        shell("cd {path}/ &&"\
               " cp ../mavisp_templates_final/plot.py .")

        #### run cancermuts depending on the input files availability ###
        env = modules["mutations_aggregation"]["cancermuts"]["source"]
        uniprot_id = df.loc[df['protein'] == wildcards.hugo_name,\
                            'uniprot_id'].iloc[0]
        uniprot_ac = df.loc[df['protein'] == wildcards.hugo_name,\
                            'uniprot_ac'].iloc[0]
        script = os.path.abspath(modules["mutations_aggregation"]\
                                        ["cancermuts"]\
                                        ["script_pancancer"])

        cancermuts_readme = modules["mutations_aggregation"]\
                                 ["cancermuts"]\
                                 ["readme"]

        cancermuts_input_script = modules["mutations_aggregation"]\
                                         ["cancermuts"]\
                                         ["script_inputs"]

        shell("cp {script}\
                  {cancermuts_readme}\
                  {cancermuts_input_script} {path}")

        if input.saturation_mutlist:
            saturation_csv = os.path.splitext(os.path.basename(input.saturation_mutlist))[0] + ".csv"

            if df_clinvar.empty:
                clinvar_option = ''
            else:
                clinvar_option = '-c clinvar.csv'

            shell("cd {path}/ &&"\
                  "python input_csv.py $(basename {input.saturation_mutlist}) && "\
                  "set +eu && . {env} && set -eu && "\
                  "mv input.csv {saturation_csv} && "\
                  "python {script} -p {wildcards.hugo_name}\
                                               -i {uniprot_id} \
                                               -a {uniprot_ac} \
                                               {clinvar_option} \
                                               -e {saturation_csv}")

        elif input_files:
            external_mutation_list=[]
            for mutlist in input_files:
                final_mutlist=mutlist.split(".")[0]+"_input.csv"
                external_mutation_list.append(final_mutlist)

        # pancancer_clinvar_others
        elif input_files and not df_clinvar.empty:
            external_mutation_list=" ".join(external_mutation_list)
            shell("cd {path} &&"\
                  " set +eu && . {env} &&"\
                  " set -eu && python {script} -p {wildcards.hugo_name}\
                                               -i {uniprot_id} \
                                               -a {uniprot_ac} \
                                               -c clinvar.csv \
                                               -e {external_mutation_list}")
        # pancancer_clinvar
        elif not input_files and not df_clinvar.empty:
            shell("cd {path} &&"\
                  " set +eu && . {env} &&"\
                  " set -eu && python {script} -p {wildcards.hugo_name} \
                                               -i {uniprot_id} \
                                               -a {uniprot_ac} \
                                               -c clinvar.csv")
        # pancancer_others
        elif input_files and df_clinvar.empty:
            external_mutation_list=" ".join(external_mutation_list)
            shell("cd {path} &&"\
                  " set +eu && . {env} &&"\
                  " set -eu && python {script} -p {wildcards.hugo_name} \
                                               -i {uniprot_id} \
                                               -a {uniprot_ac} \
                                               -e {external_mutation_list}")
        # pancancer
        elif not input_files and df_clinvar.empty:
            shell("cd {path} &&"\
                  " set +eu && . {env} &&"\
                  " set -eu && python {script} -p {wildcards.hugo_name} \
                                               -i {uniprot_id} \
                                               -a {uniprot_ac}")

################ Mutlists generation and protein annotations ################

rule mutlist:
    input:
        "{hugo_name}/structure_selection/trimmed_model/",
        f"{modules['mutations_aggregation']['cancermuts']['folder_name']}"+
        "{hugo_name}"+"/metatable_pancancer_{hugo_name}.csv"

    output:
        directory("{hugo_name}/cancermuts")
    run:
        # list with the ranges of the trimmed models to
        # obtain from the AF model
        resrange = df.loc[df['protein'] == wildcards.hugo_name,\
                                            'trimmed'].iloc[0]
        uniprot_ac = df.loc[df['protein'] == wildcards.hugo_name,\
                                        'uniprot_ac'].iloc[0]

        # list with the ranges expressed like start:end
        col_resrange=[ re.sub('-', ':', i) for i in resrange ]

        # create a list with as many uniprot_ac as
        # many trimmed models we have
        uniprot_list=[]
        for i in col_resrange:
            uniprot_list.append(uniprot_ac)

        # create a list with the trimmed models expressed as
        # uniprot_ac_start-end.pdb
        protein_list = []
        for i,l in zip(uniprot_list,resrange):
            protein_list.append(i+"_"+l+".pdb")

        pdbs = " ".join(protein_list) # put all the models in one
                                        # string separeted with a space

        ren  = " ".join(col_resrange) # put all the ranges in one
                                        # string separeted with a space

        # path with the trimmed models
        trimmed_pdb_path = f"../structure_selection/"\
                            f"trimmed_model/"
        script = modules['mutlist_generation']['script']
        readme = modules['mutlist_generation']['readme']
        shell("set +u && "\
                "{config[modules][mutations_aggregation][mutlist][source]} && "\
                "set -u && "\
                "mkdir -p {output} &&"\
                " cd {output} &&"\
                " cp ../../{input[0]}/*.pdb . &&"\
                " cp ../../{script} . &&"\
                " cp ../../{readme} . &&"\
                " python get_mutlists.py -m {input[1]}"\
                                " -d {ren}"\
                                " -M"\
                                " -R"\
                                " -H"\
                                " -p {pdbs}")

        ######### move the cancermuts file in the correct path ###########

        path = modules['mutations_aggregation']['cancermuts']['folder_name']

        mutlist_clinvar = f"{path}/{wildcards.hugo_name}/clinvar.csv"

        external_mutlist = f"{wildcards.hugo_name}/"\
                            f"external_mutation_lists"

        research_field  = df.loc[df['protein'] == wildcards.hugo_name,\
                                                'research_field'].iloc[0]
        ti_c = os.path.getmtime(input[1])
        s_ti = time.ctime(ti_c)
        t_obj = time.strptime(s_ti)
        date=time.strftime("%d%m%Y", t_obj)

        if os.path.exists(mutlist_clinvar) and\
            os.path.exists(external_mutlist):
            cancermode="pancancer_clinvar_others"

        if not os.path.exists(mutlist_clinvar) and\
            os.path.exists(external_mutlist):
            cancermode="pancancer_others"

        if not os.path.exists(external_mutlist):
            cancermode="pancancer_clinvar_saturation"

        final_path = f"{path}/{research_field}/{wildcards.hugo_name.lower()}/"\
                        f"{cancermode}/{date}"

        shell("mkdir -p {final_path} &&"\
                " mv {path}/{wildcards.hugo_name}/* {final_path} &&"\
                " rm -r {path}/{wildcards.hugo_name} &&"\
                " cp "\
                "{final_path}/metatable_pancancer_{wildcards.hugo_name}.csv\
                {wildcards.hugo_name}/cancermuts")


rule domains:
    input:
        directory("{hugo_name}/cancermuts/")
    output:
        "{hugo_name}/structure_selection/domain_annotations/"\
        "domains_mutlist.csv",
	"{hugo_name}/structure_selection/domain_annotations/"\
        "summary.csv"
    run:
        uniprot_ac = df.loc[df['protein'] == wildcards.hugo_name,\
                                                'uniprot_ac'].iloc[0]

        # save in mutlist string the mutation_list with the cancermuts
        # date obtained in the mutlist rule
        mutlist = ""
        pattern = "mutlist_\d{8}\.txt"
        script  = modules['domain_annotations']['script']
        readme  = modules['domain_annotations']['readme']
        for filename in os.listdir(str(input)):
            if re.match(pattern,filename):
                mutlist = filename

        # run the script for the domain module

        shell("mkdir -p {wildcards.hugo_name}/structure_selection/"\
                        "domain_annotations/ &&"\
                " cd {wildcards.hugo_name}/structure_selection/"\
                "domain_annotations/ && "\
                " cp ../../../{script} . &&"\
                " cp ../../../{readme} . &&"\
                " ln -snf ../../cancermuts/{mutlist} mutlist.txt &&"\
                " python ../../../{script} -u {uniprot_ac} -m mutlist.txt")

rule netphos:
     output:
         "{hugo_name}/netphos/netphos.out"
     run:
        uniprot_ac = df.loc[df['protein'] == wildcards.hugo_name,\
                                             'uniprot_ac'].iloc[0]
        path=f"{wildcards.hugo_name}/netphos"
        readme=modules['netphos']
        shell("mkdir -p {path} && cd {path} &&"\
              " cp ../../{readme} . &&"\
              " wget https://rest.uniprot.org/uniprotkb/"\
              "{uniprot_ac}.fasta &&"\
              " netphos {uniprot_ac}.fasta > netphos.out")

rule denovo_phospho:
    input:
        cancermuts_dir="{hugo_name}/cancermuts/"
    output:
        "{hugo_name}/denovo_phospho/results/aggregated_filtered_output.csv",
    threads:
        workflow.cores
    params:
        uniprot_ac = lambda wcs: df.loc[df['protein'] == wcs.hugo_name,
                                           'uniprot_ac'].iloc[0]

    run:

        # save in the mutlist string the mutation_list with the cancermuts
        # date
        pattern = "mutlist_\d{8}\.txt"
        outdir = os.path.dirname(os.path.dirname(output[0]))

        fnames = [f for f in os.listdir(input.cancermuts_dir) if re.match(pattern, f)]
        assert len(fnames) == 1
        mutation_fname = fnames[0]

        shutil.copy(f"{input.cancermuts_dir}/{mutation_fname}", outdir)
        shutil.copy(modules['denovo_phospho']['snakefile'], outdir)
        shutil.copy(modules['denovo_phospho']['readme'], outdir)
        urllib.request.urlretrieve(f"https://rest.uniprot.org/uniprotkb/{params.uniprot_ac}.fasta",
                                   f"{outdir}/{params.uniprot_ac}.fasta")

        dnp_config = {'mutlist'    : os.path.basename(mutation_fname),
                     'fasta_file' : f"{params.uniprot_ac}.fasta",
                     'output_directory' : 'results'}

        with open(f"{outdir}/config.yaml", 'w') as fh:
            yaml.dump(dnp_config, fh, default_flow_style=False)

        shell("cd {outdir} && snakemake -c {threads} --rerun-incomplete")

rule ptm_stability:
    input:
        data=lambda wcs: f"{modules['mutations_aggregation']['mutatex']['repository']}/" +\
                         f"{df.loc[df['protein'] == wcs.hugo_name, 'research_field'].iloc[0]}/{wcs.hugo_name.lower()}/" +\
                         f"free/stability/mutatex_runs/{wcs.structure_source}_" +\
                         f"{wcs.resrange}/model_" +\
                         f"{df.loc[df['protein'] == wcs.hugo_name, 'model'].iloc[0]}/saturation/" +\
                         f"{df.loc[df['protein'] == wcs.hugo_name, 'uniprot_ac'].iloc[0]}_scan/results/mutation_ddgs/final_averages/",
        pdb=lambda wcs: f"{modules['mutations_aggregation']['mutatex']['repository']}/" +\
                        f"{df.loc[df['protein'] == wcs.hugo_name, 'research_field'].iloc[0]}/{wcs.hugo_name.lower()}/" +\
                        f"free/stability/mutatex_runs/{wcs.structure_source}_" +\
                        f"{wcs.resrange}/model_" +\
                        f"{df.loc[df['protein'] == wcs.hugo_name, 'model'].iloc[0]}/saturation/" +\
                        f"{df.loc[df['protein'] == wcs.hugo_name, 'uniprot_ac'].iloc[0]}_scan/" +\
                        f"{df.loc[df['protein'] == wcs.hugo_name, 'uniprot_ac'].iloc[0]}_trimmed_model0_checked.pdb",
        mutlist_dir="{hugo_name}/cancermuts/"
    output:
        summary="{hugo_name}/ptm/{structure_source}_{resrange}/mutatex/summary_stability.txt",
    run:
        mutlist = ""
        pattern = "mutlist_mutatex_P_\d{8}\.txt"
        for filename in os.listdir(input.mutlist_dir):
            if re.match(pattern,filename):
                mutlist = filename

        outdir = os.path.dirname(output.summary)
        shutil.copy(f"{input.mutlist_dir}/{mutlist}", f"{outdir}")
        shutil.copy(modules['mutations_aggregation']['ptm']['mutatex']['readme'], outdir)
        shutil.copy(modules['mutations_aggregation']['ptm']['mutatex']['script'], outdir)
        shutil.copy(modules['mutations_aggregation']['ptm']['mutatex']['mutlist'], outdir)
        if os.path.islink(f"{outdir}/final_averages") or os.path.exists(f"{outdir}/final_averages"):
            os.remove(f"{outdir}/final_averages")
        os.symlink(input.data, f"{outdir}/final_averages")
        if os.path.islink(f"{outdir}/{os.path.basename(input.pdb)}") or os.path.exists(f"{outdir}/{os.path.basename(input.pdb)}"):
            os.remove(f"{outdir}/{os.path.basename(input.pdb)}")
        os.symlink(input.pdb,  f"{outdir}/{os.path.basename(input.pdb)}")

        shell(f"""cd {outdir} &&\
                  bash run_ddgs.sh {os.path.basename(input.pdb)} {mutlist}
               """)

rule ptm_sas:
    output:
        rsa="{hugo_name}/ptm/{structure_source}_{resrange}/naccess/{uniprot_ac}_trimmed_model0_checked.rsa"
    input:
        pdb=lambda wcs: f"{modules['mutations_aggregation']['mutatex']['repository']}/" +\
                        f"{df.loc[df['protein'] == wcs.hugo_name, 'research_field'].iloc[0]}/{wcs.hugo_name.lower()}/" +\
                        f"free/stability/mutatex_runs/{wcs.structure_source}_{wcs.resrange}/model_" +\
                        f"{df.loc[df['protein'] == wcs.hugo_name, 'model'].iloc[0]}/saturation/" +\
                        f"{df.loc[df['protein'] == wcs.hugo_name, 'uniprot_ac'].iloc[0]}_scan/" +\
                        f"{df.loc[df['protein'] == wcs.hugo_name, 'uniprot_ac'].iloc[0]}_trimmed_model0_checked.pdb"
    run:
        outdir = os.path.dirname(output.rsa)
        pdb_fname = os.path.basename(input.pdb)
        shutil.copy(modules['mutations_aggregation']['ptm']['naccess']['readme'], outdir)
        shutil.copy(input.pdb, outdir)
        shell(f"cd {outdir} && naccess {pdb_fname}")



######################### Mutations classifiers #############################

rule demask_config:
    output:
        "{hugo_name}/demask/config.ini"
    run:
        # Demask: Create the cofig.ini file
        uniprot_ac = df.loc[df['protein'] == wildcards.hugo_name,\
                                             'uniprot_ac'].iloc[0]
        readme = modules['mutations_classifier']['demask']['readme']
        script = modules['mutations_classifier']['demask']['script']
        config = configparser.ConfigParser()
        config.add_section('demask')
        config.set('demask', 'blastp', modules['mutations_classifier']\
                                              ['demask']\
                                              ['blastp_path'])
        config.set('demask', 'db', modules['mutations_classifier']\
                                          ['demask']\
                                          ['database_path'])
        if not os.path.exists(f"{wildcards.hugo_name}/demask"):
            os.mkdir(f"{wildcards.hugo_name}/demask")
        shell("cp {readme} {wildcards.hugo_name}/demask && "\
              "cp {script} {wildcards.hugo_name}/demask")

        with open(os.path.join(wildcards.hugo_name, "demask",
                                              "config.ini"), "w") as f:
            config.write(f)

rule demask_homologs:
    input:
        "{hugo_name}/demask/config.ini"
    output:
        "{hugo_name}/demask/myquery_homologs.a2m"
    params:
        n_threads = 4,
        uniprot_ac = lambda wcs: df.loc[df['protein'] == wcs.hugo_name,
                                           'uniprot_ac'].iloc[0]

    shell:
        """
        cd {wildcards.hugo_name}/demask
        wget https://rest.uniprot.org/uniprotkb/{params.uniprot_ac}.fasta
        set +eu
        source {config[modules][mutations_classifier][demask][source]}
        python3 -m demask.homologs \
        -s {params.uniprot_ac}.fasta \
        -o myquery_homologs.a2m \
        -c config.ini \
        -t {params.n_threads}\
        """

rule demask_prediction:
    input:
        "{hugo_name}/demask/myquery_homologs.a2m"
    output:
        "{hugo_name}/demask/myquery_predictions.txt"
    params:
        n_threads = 4,
        uniprot_ac = lambda wcs: df.loc[df['protein'] == wcs.hugo_name,
                                           'uniprot_ac'].iloc[0]

    shell:
        """
        cd {wildcards.hugo_name}/demask
        wget https://rest.uniprot.org/uniprotkb/{params.uniprot_ac}.fasta
        set +eu
        source {config[modules][mutations_classifier][demask][source]}
        python3 -m demask.predict \
        -i myquery_homologs.a2m \
        -o myquery_predictions.txt
        """

rule alphamissense:
    output:
        "{hugo_name}/alphamissense/am.tsv.gz"
    params:
        uniprot_ac = lambda wcs: df.loc[df['protein'] == wcs.hugo_name,
                                           'uniprot_ac'].iloc[0]
    shell:
        """
        cd $(dirname {output})
        cp ../../{config[modules][mutations_classifier][alphamissense][readme]} .
        cp ../../{config[modules][mutations_classifier][alphamissense][script]} .
        bash do.sh {params.uniprot_ac}
        """

############################## Calculations #################################
rule rasp_workflow:
    input:
        lambda wcs: f"{wcs.hugo_name.upper()}/structure_selection/trimmed_model/",
    output:
        directory("{path}{research_field}/{hugo_name}/free/{structure_source}_{resrange}/{model}_model/")
    shell:
        """
        mkdir -p {output}
        cp {config[modules][rasp][readme]} {output}

        pdb_name=$(ls {input}/*_{wildcards.resrange}.pdb | xargs basename)
        processed_pdb=$(echo "$pdb_name" | sed 's/\.pdb$/_processed.pdb/')

        cp {input}/"$pdb_name" {output}/"$pdb_name"

        cd {output}

        # Step 1: Add chain if missing
        add_chain() {{
            awk -v chain="A" '{{
                if ($0 ~ /^ATOM/ || $0 ~ /^HETATM/) {{
                    if (substr($0, 22, 1) == " ") {{
                        print substr($0, 1, 21) chain substr($0, 23);
                    }} else {{
                        print $0;
                    }}
                }} else {{
                    print $0;
                }}
            }}' "$pdb_name" > tmp_chain.pdb
        }}

        # Step 2: Remove solvent and hydrogen
        remove_solvent_and_hydrogen() {{
            pdb_element tmp_chain.pdb | pdb_delelem -H > tmp_noH.pdb
        }}

        # Step 3: Fix amino acid names
        fix_amino_acid_names() {{
            sed -E 's/HIE/HIS/g; s/HID/HIS/g; s/HIP/HIS/g; s/LYN/LYS/g; s/ASH/ASP/g; s/GLH/GLU/g; s/CYX/CYS/g' tmp_noH.pdb > tmp_fixed_names.pdb
        }}

        # Step 4: Remove ending lines
        remove_end_lines() {{
            awk '/^(ATOM|HETATM)/' tmp_fixed_names.pdb > tmp_no_end_lines.pdb
        }}

        # Step 5: Pad the file
        pad_file() {{
            sed -E 's/.{{1,79}}$/&                                                                            /' tmp_no_end_lines.pdb | cut -c1-80 > "$processed_pdb"
        }}

        # Execute workflow
        add_chain "$pdb_name"
        remove_solvent_and_hydrogen
        fix_amino_acid_names
        remove_end_lines
        pad_file "$processed_pdb"

        # Clean up
        rm tmp_chain.pdb tmp_noH.pdb tmp_fixed_names.pdb tmp_no_end_lines.pdb

        # Activate conda environment
        set +u; source {config[modules][rasp][conda_activation]}
        conda activate {config[modules][rasp][rasp_conda_env]}; set -u

        # Run RaSP workflow
        RaSP_workflow -i $processed_pdb \
                      -r cpu \
                      -p /usr/local/envs/RaSP_workflow/RaSP_workflow/src/ \
                      -o . \
                      -n 4 \
                      -c A
        RaSP_postprocess -i output/predictions/cavity_pred_*.csv
        """
rule rosetta_relax:
    input:
        "{hugo_name}/structure_selection/trimmed_model/"
    output:
        "{path}/{research_field}/{hugo_name}/free/{structure_source}_{resrange}/{model}_model/"\
        "ref2015_cartesian2020/relax/relax_{uniprot_ac}_{resrange}_0001.pdb"
    params:
        rosetta_module = modules['rosetta_relax']['rosetta_module'],
        rosetta_folder = modules['rosetta_relax']['rosetta_folder'],
        mpi = modules['rosetta_relax']['rosettampi_yaml'],
        readme = modules['rosetta_relax']['readme'],
        yaml = modules['rosetta_relax']['relax_yaml']
    shell:
        """
        set +u; source {config[modules][rosetta_relax][rosetta_env]}; set -u && \
        mkdir -p {config[modules][rosetta_relax][rosetta_folder]}/{wildcards.research_field}/{wildcards.hugo_name}/free/{wildcards.structure_source}_{wildcards.resrange}/{wildcards.model}_model/ref2015_cartesian2020/ && \
        cp {input}/{wildcards.uniprot_ac}_{wildcards.resrange}.pdb \
           {config[modules][rosetta_relax][rosetta_folder]}/{wildcards.research_field}/{wildcards.hugo_name}/free/{wildcards.structure_source}_{wildcards.resrange}/{wildcards.model}_model/ref2015_cartesian2020/ && \
        cd {config[modules][rosetta_relax][rosetta_folder]}/{wildcards.research_field}/{wildcards.hugo_name}/free/{wildcards.structure_source}_{wildcards.resrange}/{wildcards.model}_model/ref2015_cartesian2020/ && \
        cp {params.mpi} . && \
        cp {params.readme} . && \
        cp {params.yaml} . && \
        rosetta_ddg_run -p {wildcards.uniprot_ac}_{wildcards.resrange}.pdb \
                        -cr {params.yaml} \
                        -n 1 \
                        -r {params.rosetta_module}\
                        -cs {params.mpi}
        """

'''
        expand("{hugo_name}/long_range/"\
               "allosigma2/"\
               "{uniprot_ac}_{resrange}/"\
               "4.allosigma_filtering/"\
               "cmpsn_all.csv",
               zip, hugo_name = df_exploded['protein'],
                                            uniprot_ac = \
                                            df_exploded['uniprot_ac'],
                                            resrange = \
                                            df_exploded['trimmed']),

rule allosigma_1:
    output:
        directory("{hugo_name}/long_range/allosigma2/"\
                  "{uniprot_ac}_{resrange}/1.allosteric_signalling_map/")
    run:
        # 1.allosigma_signalling_map

        allosigma_file=df.loc[df['protein']==wildcards.hugo_name,\
                                             'allosigma_file'].iloc[0]

        # prepare the allosigma runs

        readme = modules["allosigma"]["allosigma1"]["readme"]

        shell("mkdir -p {output}/ && \
              cp {readme} {output} && \
              cp {allosigma_file} {output}")

rule allosigma_2_3:
    input:
        directory("{hugo_name}/cancermuts"),
        directory("{hugo_name}/long_range/allosigma2/"\
                  "{uniprot_ac}_{resrange}/1.allosteric_signalling_map/"),


    output:
        "{hugo_name}/long_range/allosigma2/{uniprot_ac}_{resrange}/"\
        "2.allosigma_classify/allosigma_mut.txt",
        "{hugo_name}/long_range/allosigma2/{uniprot_ac}_{resrange}/"\
        "3.allosigma_heatmap/all/up_mutations.tsv",
        "{hugo_name}/long_range/allosigma2/{uniprot_ac}_{resrange}/"\
        "3.allosigma_heatmap/all/down_mutations.tsv",
    run:

        output_dest="/".join(input[1].split("/")[:-1])+\
                    "/2.allosigma_classify"

        readme_2 = modules["allosigma"]["allosigma2"]["readme"]
        readme_3 = modules["allosigma"]["allosigma3"]["readme"]

        output_dest2="/".join(input[1].split("/")[:-1])+\
                     "/3.allosigma_heatmap/all"
        allosigma_file=df.loc[df['protein']==wildcards.hugo_name,\
                                 'allosigma_file'].iloc[0]

        date=""
        directory=os.getcwd()
        for filename in os.listdir(input[0]):
            match=re.search(r'(_\d{8}\.)',filename)
            if match!= None:
                date=match.group()[1:-1]

        classify=modules['allosigma']['allosigma2']['script']
        aminoacids=modules['allosigma']['allosigma1']['aminoacids']
        heatmap=modules['allosigma']['allosigma3']['script']

        shell("mkdir -p {output_dest} && \
               mkdir -p {output_dest2} && \
               cp {heatmap} {output_dest2} && \
               cp {readme_2} {output_dest} && \
               cp {readme_3} {output_dest2} && \
               cp {input[0]}/mutlist_{date}.txt {output_dest}/muts.dat && \
               cp {classify} {output_dest} && cp {aminoacids} {output_dest}")

        shell("cd {output_dest} &&  \
              ./allosigma-classify \
              -o {directory}/{output[0]} \
              -c 5 muts.dat aminoacids.dat &&"\
              " cd ../../../../../{output_dest2} && "\
              "./allosigma-heatmap -x 10 \
                                   -y 40  \
                                   -f 6  \
                                   ../../1.allosteric_signalling_map/"\
                                 "{allosigma_file}  {directory}/{output[0]}")
        # to add a checking step for the presence of the heatmaps-->
        # no heatmaps if the ziop file is thre wrong one

rule allosigma4:
    input:
        directory("{hugo_name}/structure_selection/trimmed_model/"),
        "{hugo_name}/long_range/allosigma2/{uniprot_ac}_{resrange}/"\
        "3.allosigma_heatmap/all/up_mutations.tsv",
        "{hugo_name}/long_range/allosigma2/{uniprot_ac}_{resrange}/"\
        "3.allosigma_heatmap/all/down_mutations.tsv",

    output:
        "{hugo_name}/long_range/allosigma2/{uniprot_ac}_{resrange}/"\
        "4.allosigma_filtering/cmpsn_all.csv"
    run:
        output_dest4 = "/".join(str(output).split("/")[:-1])
        directory = os.getcwd()
        filtering = modules['allosigma']['allosigma4']['script']
        readme = modules['allosigma']['allosigma4']['readme']

        # 4.allosigma.filtering
        shell("mkdir -p {output_dest4} && \
               cp {filtering} {output_dest4} &&\
               cp {readme} {output_dest4} && \
               cp {input[0]}/{wildcards.uniprot_ac}_{wildcards.resrange}.pdb\
                {output_dest4} && cd {output_dest4} && "\
              "cp {directory}/{input[1]} {directory}/{input[2]} . && "\
              "touch {directory}/{output} && "\
              "./allosigma-filtering -s "\
              "{wildcards.uniprot_ac}_{wildcards.resrange}.pdb \
              -i up_mutations.tsv -t 2 -d 8 -a 20 && "\
              "./allosigma-filtering -s "\
              "{wildcards.uniprot_ac}_{wildcards.resrange}.pdb \
              -i down_mutations.tsv -t 2 -d 8 -a 20")
'''

rule metadata:
    output:
        metadata  = "{hugo_name}/metadata/metadata.yaml",
        importing = "{hugo_name}/metadata/importing.yaml"
    run:
        import getpass
        from datetime import datetime

        gene    = wildcards.hugo_name
        meta_dir = os.path.dirname(output.metadata)
        os.makedirs(meta_dir, exist_ok=True)

        # ---- write metadata.yaml ----
        curator_map = {
            'pablosanchezb': {
                'full_name': 'Pablo Sanchez-Izquierdo',
                'affiliation': ['DTU, Denmark', 'DCI, Denmark']
            },
        }
        user  = getpass.getuser()
        entry = curator_map.get(user, {
            'full_name': user,
            'affiliation': ['<Your Affiliation>']
        })

        uniprot_ac       = df.loc[df['protein'] == gene, 'uniprot_ac'].iloc[0]
        refseq_id        = df.loc[df['protein'] == gene, 'ref_seq'].iloc[0]
        structure_source = df.loc[df['protein'] == gene, 'structure_source'].iloc[0]
        pdb_id           = df.loc[df['protein'] == gene, 'input_pdb'].fillna('').iloc[0]

        with open(output.metadata, 'w') as out:
            out.write("curators:\n")
            out.write(f"  {entry['full_name']}:\n")
            out.write("    affiliation:\n")
            for aff in entry['affiliation']:
                out.write(f"      - {aff}\n")
            out.write(f"uniprot_ac: {uniprot_ac}\n")
            out.write(f"refseq_id: {refseq_id}\n")
            out.write("allosigma_distance_cutoff: [ 15 ]\n")
            out.write("review_status: 0\n")
            out.write(f"structure_source: {structure_source}\n")
            out.write("linker_design: False\n")
            out.write(f"pdb_id: {pdb_id}\n")

        # ---- write importing.yaml ----
        project  = df.loc[df['protein'] == gene, 'research_field'].iloc[0]
        trimmed_list  = df.loc[df['protein'] == gene, 'trimmed'].iloc[0]
        starting_list = [f"{structure_source}_{part}" for part in trimmed_list]

        base_dir  = (
            "/data/raw_data/computational_data/cancermuts_data/"
            f"{project}/{gene.lower()}/pancancer_clinvar_saturation"
        )
        candidates = []
        for name in os.listdir(base_dir):
            path = os.path.join(base_dir, name)
            if os.path.isdir(path):
                try:
                    dt = datetime.strptime(name, "%d%m%Y")
                    candidates.append((dt, name))
                except ValueError:
                    pass
        if not candidates:
            raise ValueError(f"No valid date folders in {base_dir}")
        date = max(candidates, key=lambda x: x[0])[1]

        with open(output.importing, 'w') as out:
            out.write(f"gene: {gene}\n")
            out.write(f"project: {project}\n")
            out.write(f"date: {date}\n")
            out.write("cancermuts: pancancer\n")
            for idx, st in enumerate(starting_list, 1):
                out.write(f"starting_structure{idx}: {st}\n")

rule collect_outputs:
    input:
        clinvar_genes=lambda wcs: f"{wcs.hugo_name}/clinvar_gene/genes_output.csv",
        ptm_stability=lambda wcs: f"{wcs.hugo_name}/ptm/{wcs.structure_source}_{wcs.resrange}/mutatex/summary_stability.txt",
        ptm_sas=lambda wcs: f"{wcs.hugo_name}/ptm/{wcs.structure_source}_{wcs.resrange}/naccess/{wcs.uniprot_ac}_trimmed_model0_checked.rsa",
        demask=lambda wcs: f"{wcs.hugo_name}/demask/myquery_predictions.txt",
        alphamissense=lambda wcs: f"{wcs.hugo_name}/alphamissense/am.tsv.gz",
        cancermuts=lambda wcs: f"{wcs.hugo_name}/cancermuts/",
        efoldmine=lambda wcs: f"{wcs.hugo_name}/efoldmine/{wcs.uniprot_ac}.tabular",
        structure_rasp=lambda wcs: f"/data/raw_data/computational_data/rasp_data/{wcs.research_field}/{wcs.hugo_name.lower()}/free/{wcs.structure_source}_{wcs.resrange}/{wcs.model}_model",
        structure_foldx5=lambda wcs: f"/data/raw_data/computational_data/mutatex_data/{wcs.research_field}/{wcs.hugo_name.lower()}/free/stability/mutatex_runs/{wcs.structure_source}_{wcs.resrange}/model_{wcs.model}/saturation/{wcs.uniprot_ac}_table/energies.csv",
	pfam=lambda wcs: f"{wcs.hugo_name}/structure_selection/domain_annotations/summary.csv",
	alphafold=lambda wcs: f"{wcs.hugo_name}/structure_selection/original_model/{wcs.hugo_name.lower()}/{wcs.uniprot_ac}.csv"
    output:
        temp("{hugo_name}/simple_mode/collection_{research_field}_{structure_source}_{resrange}_{uniprot_ac}_{model}.done")
    shell:
        """
        mkdir -p {wildcards.hugo_name}/simple_mode/clinvar
        cp {input.clinvar_genes} {wildcards.hugo_name}/simple_mode/clinvar/variants_output.csv

        mkdir -p {wildcards.hugo_name}/simple_mode/ptm
        cp {input.ptm_stability} {wildcards.hugo_name}/simple_mode/ptm/summary_stability.txt
        cp {input.ptm_sas} {wildcards.hugo_name}/simple_mode/ptm/sasa.rsa
	cp {input.cancermuts}/metatable_pancancer_{wildcards.hugo_name}.csv {wildcards.hugo_name}/simple_mode/ptm/metatable.csv

        mkdir -p {wildcards.hugo_name}/simple_mode/sas
        cp {input.ptm_sas} {wildcards.hugo_name}/simple_mode/sas/sasa.rsa

        mkdir -p {wildcards.hugo_name}/simple_mode/demask
        cp {input.demask} {wildcards.hugo_name}/simple_mode/demask/myquery_predictions.txt

        mkdir -p {wildcards.hugo_name}/simple_mode/alphamissense
        cp {input.alphamissense} {wildcards.hugo_name}/simple_mode/alphamissense/am.tsv.gz

        mkdir -p {wildcards.hugo_name}/simple_mode/cancermuts
        cp {input.cancermuts}/metatable_pancancer_{wildcards.hugo_name}.csv {wildcards.hugo_name}/simple_mode/cancermuts/

        mkdir -p {wildcards.hugo_name}/simple_mode/efoldmine
        cp {input.efoldmine}    {wildcards.hugo_name}/simple_mode/efoldmine/

        mkdir -p {wildcards.hugo_name}/simple_mode/stability/{wildcards.structure_source}_{wildcards.resrange}/{wildcards.structure_source}/model_{wildcards.model}/rasp
        cp {input.structure_rasp}/output/predictions/post_processed_*.csv \
        {wildcards.hugo_name}/simple_mode/stability/{wildcards.structure_source}_{wildcards.resrange}/{wildcards.structure_source}/model_{wildcards.model}/rasp/

        mkdir -p {wildcards.hugo_name}/simple_mode/stability/{wildcards.structure_source}_{wildcards.resrange}/{wildcards.structure_source}/model_{wildcards.model}/foldx5
        cp {input.structure_foldx5} \
        {wildcards.hugo_name}/simple_mode/stability/{wildcards.structure_source}_{wildcards.resrange}/{wildcards.structure_source}/model_{wildcards.model}/foldx5/energies.csv
	
        cp {wildcards.hugo_name}/metadata/metadata.yaml    {wildcards.hugo_name}/simple_mode/metadata.yaml

	mkdir -p {wildcards.hugo_name}/simple_mode/pfam
	cp {input.pfam} {wildcards.hugo_name}/simple_mode/pfam/summary.csv
	
	mkdir -p {wildcards.hugo_name}/simple_mode/alphafold
	cp {input.alphafold} {wildcards.hugo_name}/simple_mode/alphafold/{wildcards.uniprot_ac}.csv
	
    	mkdir -p {wildcards.hugo_name}/simple_mode/mutation_list
    	f=$(ls {wildcards.hugo_name}/cancermuts/mutlist_[0-9][0-9]*.txt \
         	| grep -E '^.*/mutlist_[0-9]{{8}}\.txt$')
    	b=$(basename "$f")
    	awk 'BEGIN{{print "mutation\tPMID"}} {{print $0 "\thttps://doi.org/10.1101/2022.10.22.513328"}}' \
        	"$f" \
    	> {wildcards.hugo_name}/simple_mode/mutation_list/${{b/mutlist_/mutations_pmid_}}
	
	touch {output}
        """

