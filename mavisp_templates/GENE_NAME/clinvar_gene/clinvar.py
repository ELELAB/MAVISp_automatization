import re
import requests
import sys
import pandas as pd
import numpy as np
import xmltodict
import pprint
import time
import traceback
import argparse
from tabulate import tabulate

d = {'Cys': 'C', 'Asp': 'D', 'Ser': 'S', 'Gln': 'Q', 'Lys': 'K',
     'Ile': 'I', 'Pro': 'P', 'Thr': 'T', 'Phe': 'F', 'Asn': 'N',
     'Gly': 'G', 'His': 'H', 'Leu': 'L', 'Arg': 'R', 'Trp': 'W',
     'Ala': 'A', 'Val':'V', 'Glu': 'E', 'Tyr': 'Y', 'Met': 'M'}

stars = {'practice guideline':'4', 
         'reviewed by expert panel':'3', 
         'criteria provided, multiple submitters, no conflicts':'2',
         'criteria provided, multiple submitters': '1.2',
         'criteria provided, conflicting interpretations':'1',
         'criteria provided, conflicting classifications':'1',
         'criteria provided, single submitter':'1',
         'no assertion for the individual variant':'0', 
         'no interpretation for the single variant':'0',
         'no assertion criteria provided':'0', 
         'no assertion provided':'0',
         'no classification provided':'0'}



#----------------------------------------- utility functions -----------------------------------------------#

def mutations_convert(string,dic,print_message):

    '''Convert the mutation in three letter code in a one letter code mutation
    
    The function takes a string of a mutation expressed in the following format 
    p.[A-Z][a-z][a-z][0-9]+[A-Z][a-z][a-z] (Clinvar annotation i.e p.Ala128Ile)
    and returns the mutation in the following format: [A-Z][0-9][A-Z].

    Parameters
    ----------
    string: string
           mutation expressed in three letters code (i.e p.Ala128Ile)
    dic: dictionary
           dictionary containing as keys three letters annotated aminoacids
           and values one letters annotated aminoacids

    Returns
    -------
    string containing the mutation expresed as one letter code.

    '''

    #this matches strings like p.Ala128Ile and creates one group 
    #for each part of the string we need
    match = re.search("p\.([A-Z][a-z][a-z])([0-9]+)([A-Z][a-z][a-z])", string)
    match_output_df=re.search("([A-Z])([0-9]+)([A-Z])", string)
    if match is not None and len(match.groups()) == 3:
        #some variants are annotated as follow: p.Glu4Ter. 
        #So we need the if statement to match the second element 
        #after the number.
        if match.group(3) in d.keys():
            substring_list = list(match.groups())
            return f"{dic[substring_list[0]]}{substring_list[1]}{dic[substring_list[2]]}"
    
    else:
        if match_output_df is None and print_message==True :
            print(f"WARNING: {string} the function requires the variant expressed with"\
                " the three letters code, another annotation could have been used in the Clinvar database")

def from_321(string,dic):

    '''Convert the mutation in one letter code in a three letter code mutation
    
    The function takes a string of a mutation expressed in the following 
    format [A-Z][0-9]+[A-Z] (Clinvar annotation i.e A128I)
    and returns the mutation in the following format: 
    p.[A-Z][a-z][a-z][0-9]+[A-Z][a-z][a-z].

    Parameters
    ----------
    string: string
           mutation expressed in one letters code (A128I)
    dic: dictionary
           dictionary containing as keys three letters annotated aminoacids
           and values one letters annotated aminoacids

    Returns
    -------
    string containing the mutation expresed as three letter code.

    '''
    # this matches strings like A128I and creates one group 
    # for each part of the string we need
    match = re.search("([A-Z])([0-9]+)([A-Z])", string) 
    if match is not None and len(match.groups()) == 3:
        #some variants can have another symbol after the position number. 
        #So we need the if statement to match the second element 
        #after the number.
        if match.group(3) in d.values():
            substring_list = list(match.groups())
            wt_residue     = {i for i in dic if dic[i] == substring_list[0]}
            mutant_residue = {i for i in dic if dic[i] == substring_list[2]}
            return "(p." + str(wt_residue)[2:-2] + substring_list[1] + str(mutant_residue)[2:-2] + ")"

def melting_dictionary(variants_annotation):
    '''
    Count the number of keys in the second element of variants annotation values
    dictionary, create as many lists as the number of counted keys, inserting in each list
    the corresponding value of the key and populate the output list.

    The function takes as input a dictionary containing the ClinVar ID as keys and a list
    with the variant information from ClinVar as value organized in dictionaries: mutation, 
    gene, classification, condition, review status, and methods. The function counts how 
    many classifications (keys) are reported in the classification, condition, and review 
    status dictionaries and after checking that the keys per position are correct create as 
    many lists as the number of classifications (keys) in which each list contains the corresponding 
    element of the key. All the generated lists are grouped in the output_list.

    Parameters
    ----------
    variants_annotation: dict
        Input dictionary containing the ClinVar ID as key and a list with the following elements 
        organized in dictionaries as value: mutation, gene, classification, condition, 
        review status, methods.

    Returns
    -------
    output_lists: list of lists
        List of lists from the variants_annotation dictionary.
    '''
    output_list = []
    for clinvar_id in variants_annotation.keys():
        for index in range(len(variants_annotation[clinvar_id][2][0].keys())):
            classification_check = []
            classification_check.append(list(variants_annotation[clinvar_id][2][0].keys())[index])
            classification_check.append(list(variants_annotation[clinvar_id][3][0].keys())[index])
            classification_check.append(list(variants_annotation[clinvar_id][4][0].keys())[index])
            if len(set(classification_check)) == 1:
                classification = list(variants_annotation[clinvar_id][2][0].values())[index]
                condition = list(variants_annotation[clinvar_id][3][0].values())[index]
                review_status = list(variants_annotation[clinvar_id][4][0].values())[index]
                features = [clinvar_id,
                            variants_annotation[clinvar_id][0],
                            variants_annotation[clinvar_id][1],
                            [classification],
                            condition,
                            [review_status]]
                if args.methods:
                    features.append(variants_annotation[clinvar_id][5][0])
                output_list.append(features)

    return output_list


#----------------------------------------------- URL query functions -----------------------------------------------------#


def URL_response_check(URL, error_message, function):

    ''' check the response from the queried URL

    The function tries to access to the given URL. if it succeeds, 
    it will return the response, 
    otherwise it'll try again for 200 times before giving error.

    Parameters
    ----------
    URL: string
         contains the URL to be checked

    Returns
    -------
    response 
        output from the function requests.get()
    '''
    max_retries = 200
    attempts = 0
    delay = 3
    req_succeeded = False
    while not req_succeeded and attempts < max_retries:
        if attempts== max_retries:
            print("ERROR: request failed after 155 attempts; exiting...")
            exit(1)
        response = requests.get(URL)
        if response.status_code != 200:
            print(f"WARNING: {function} failed to access to {error_message};"\
                   f" will try again in {delay} seconds (attempt {attempts+1}/{max_retries})")
            req_succeeded = False
            attempts +=1
            time.sleep(delay)
        else:
            req_succeeded = True
    if not req_succeeded:
        print("ERROR: request failed after 155 attempts; exiting...")
        exit(1)
    return response


def VCV_summary_retriever(clinvar_id):
    '''
    Retrieve the XML summary file associated with a specific ClinVar ID.

    The function takes a clinvar_id as input, retrieves the corresponding VCV
    code in order to retrieve the summary XML file associated with that 
    clinvar_id, returning it as output.

    Parameters
    ----------
    clinvar_id: str
        ClinVar ID code to access the VCV code from which to retrieve the 
        summary XML associated with a specific ClinVar ID.

    Returns
    -------
    parse_VCV: OrderedDict
        XML summary output associated with the VCV code provided as input.
    '''
    URL_summary="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=clinvar&id="+clinvar_id
    error_message=f"the variant name in the Clinvar summary for the variant_id {clinvar_id}"
    function = "VCV_summary_retriever"
    parse_summary=xmltodict.parse(URL_response_check(URL_summary,error_message,function).content)
    VCV=(parse_summary['eSummaryResult']\
                      ["DocumentSummarySet"]\
                      ['DocumentSummary']['accession'])
    URL_VCV="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=clinvar&rettype=vcv&id="+VCV
    error_message=f"the XML for the following VCV accession: {VCV}"
    parse_VCV=xmltodict.parse(URL_response_check(URL_VCV,error_message,function).content)

    return parse_VCV


def filtered_variants_extractor(URL_filter,gene,isoform,mutation_type):

    '''
    Retrieve all the variant IDs belonging to certain ClinVar classes (e.g., missense,
    pathogenic, likely pathogenic, likely benign, or benign classifications) from an XML file
    associated with a specific query on the ClinVar database.

    The function takes as input a URL filter with which to perform the specific query,
    the gene name in the Hugo name format, and one of the RefSeq codes associated 
    with that gene. It extracts all the ClinVar IDs whose mutations are classified 
    according to the filter provided as input and returns a dictionary with the classifications
    as keys and a list of variant IDs belonging to missense mutations as values.

    Parameters
    ----------
    gene: str
        Input gene name in the Hugo name format.
    isoform: str
        RefSeq belonging to that gene.
    url_filter: str
        URL filter used to perform the specific query on ClinVar.

    Returns
    -------
    classified_missense_ids: dict 
        Dictionary in which the keys are the classifications from ClinVar database and
        the values are lists of variant IDs belonging to missense mutations.
    '''

    function="filtered_variants_extractor" 

    ##################################################################################################################################################
    #                                                                                                                                                #
    #                           Extraction of ids with a ClinVar classification belonging to a mutations category                                    #
    #                                                                                                                                                #
    ##################################################################################################################################################

    
    URL="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=clinvar&term="+URL_filter
    error_message=f"the XML file with the total number of {mutation_type} mutations for {gene} gene"
    if "clinvar id" in mutation_type:
        error_message=f"the XML file with the total number of missense mutation in {gene} gene"

    first_search=xmltodict.parse(URL_response_check(URL,error_message,function).content)

    #---------------------------- extract the ClinVar ids associated to a mutations with a ClinVar classification -------------------------------#

    # return clinvar_id if there are variants associated to the query

    if first_search["eSearchResult"]["Count"] != "0":
        if "clinvar_id" in mutation_type:
            clinvar_ids = variant_ids_extractors(first_search)
        else:
            tot_variant=first_search["eSearchResult"]["Count"]
            error_message=f"the number of missense variants annotated in Clinvar for {gene} gene"
            filtered_URL="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=clinvar&term="+URL_filter+"&retmax="+str(tot_variant)
            parse_search=xmltodict.parse(URL_response_check(filtered_URL,error_message,function).content)
            clinvar_ids = variant_ids_extractors(parse_search)

            print("Gene "+gene+" has " +str(tot_variant)+" "+mutation_type+" variants.")

        return clinvar_ids,""
                    
    # if no variants are annotated in Clinvar Database the script return a WARNING message
    else:
        if "missense" in URL_filter and "clinvar_id" not in mutation_type:
            print("No missense variants are annotated in Clinvar Database for "+gene+" "+\
                  "gene, the variants provided as input will be annotated in entry_not_found.csv file.")

        return [],gene



#------------------------------------- functions to parse the ClinVar Xml file --------------------------------------------#



def coding_region_variants_extractor(clinvar_VCV_xml,clinvar_code):
    '''
    Retrieve the HGVS annotations for a given VCV XML file associated with a 
    ClinVar variant ID.

    The function takes as input a VCV XML file from a ClinVar variant ID and 
    returns the list of HGVS annotations belonging to coding and protein 
    expression types.

    Parameters
    ----------
    clinvar_VCV_xml: file-like object
        XML file from a VCV ClinVar entry.
    clinvar_id: str
        ClinVar variant ID associated with the clinvar_VCV_xml variable.

    Returns
    -------
    hgvss: list 
        List of HGVS annotations belonging to coding and protein expression 
        types.
    '''

    try:
        simple_allele_accession = clinvar_VCV_xml['ClinVarResult-Set']\
                                                 ['VariationArchive']\
                                                 ['ClassifiedRecord']\
                                                 ['SimpleAllele']\
        
        if "HGVSlist" in simple_allele_accession.keys():
            try:
                HGVS_accession = clinvar_VCV_xml['ClinVarResult-Set']\
                                                ['VariationArchive']\
                                                ['ClassifiedRecord']\
                                                ['SimpleAllele']\
                                                ["HGVSlist"]\
                                                ["HGVS"]
            except:
                print(simple_allele_accession["HGVSlist"].keys())

            try:
                hgvss = [ x for x in HGVS_accession if x['@Type'] == 'coding' and "ProteinExpression" in x.keys()]
                    
                # there are cases in which the field "InterpretedRecord" is "Included Record", so epeat the previous step with the 
                # new key

            except Exception:
                print("The clinvar_id "+clinvar_code+ " does not have any annotation referring to variant in an expressed protein coding region")
                hgvss=[]
        else:
            hgvss=[]
    except:
        print("The clinvar_id "+clinvar_code+ " has a different annotation structure. It will be annotated in variants_to_check.csv")
        hgvss=[]


    return hgvss

def missense_variants_extractor(hgvss, gene, clinvar_id, isoform_to_check):
    '''
    Given an HGVS list from a VCV XML file associated with a ClinVar variant ID, 
    returns the coding mutations in HGVS format.

    The function takes as input an HGVS list obtained from a VCV XML file, the 
    corresponding gene, the ClinVar variant ID associated with the VCV XML file, 
    and an isoform associated with the input gene. It extracts the coding mutations 
    associated with the given gene in the HGVS format.

    Parameters
    ----------
    hgvss: list
        HGVS list belonging to coding and protein expression type.
    gene: str
        Gene to which the mutations belong.
    clinvar_id: str
        ClinVar variant ID associated with the VCV XML file from which the HGVS 
        has been extracted.
    isoform_to_check: str
        RefSeq associated with the gene.

    Returns
    -------
    correct_variant: str
        Missense coding variant expressed in the HGVS format.
    '''
    correct_variant = ""
    for mut_info in hgvss:
        #look for isoform identifier in the hgvss list of dictionaries
        if "@sequenceAccession" in mut_info["ProteinExpression"]:
            
            identifier = mut_info["ProteinExpression"]["@sequenceAccession"]
            variant = mut_info["ProteinExpression"]["@change"]            
            # use the information associated to the right isoform to build the correct entry
            if identifier == isoform_to_check:
                index = mut_info['NucleotideExpression']['Expression'].rfind(":")
                correct_variant = correct_variant + mut_info['NucleotideExpression']['Expression'][:index]+ f"({gene})"+\
                                  mut_info['NucleotideExpression']['Expression'][index:]+ f" ({variant})"
                return correct_variant

            else:
                print(f"{identifier} associated to the variant_id {clinvar_id}"\
                        f" is an alternative isoform of {gene} gene, only the mutations"\
                        f" belonging to the {isoform_to_check} isoform will be considered")

def variant_ids_extractors(clinvar_VCV_xml):
    '''
    Retrieve all the variant IDs from an XML file associated with a specific query on ClinVar.

    The function takes as input a VCV XML file from a ClinVar variant ID and 
    returns the list of ClinVar IDs contained in the XML file and associated
    with a specific query on the ClinVar database.

    Parameters
    ----------
    clinvar_VCV_xml: file-like object
        XML file from a VCV ClinVar entry.

    Returns
    -------
    ids: list 
        List with the ClinVar IDs of interest.
    '''
    
    ids = clinvar_VCV_xml["eSearchResult"]\
                         ["IdList"]\
                         ["Id"]
    if isinstance(ids, list):
        return ids
    else:
        return [ids]

def classification_methods_extractor(clinvar_VCV_xml):
    '''
    Retrieve all the methods used for the variants classification
    reported in the XML file associated with a specific ClinVar ID
    from the ClinVar database.

    The function takes as input a VCV XML file from a ClinVar variant ID and 
    returns the list of methods used for the classification of the variant 
    associated with that specific variant ID.

    Parameters
    ----------
    clinvar_VCV_xml: xml file
        XML file from a VCV ClinVar entry 

    Returns
    -------
    methods: list 
        List of methods associated with a specific variant ID
    '''
    method = clinvar_VCV_xml['ClinVarResult-Set']\
                            ['VariationArchive']\
                            ['ClassifiedRecord']\
                            ['ClinicalAssertionList']\
                            ['ClinicalAssertion']
    try:
        if isinstance(method, list):
            methods = []
            for i in method:
                methods.append(i['AttributeSet']['Attribute']['#text'])
            return list(set(methods))
        else:
            return [method['AttributeSet']['Attribute']['#text']]

    except:
        return ["not provided"],clinvar_VCV_xml


def conditions_extractor(clinvar_VCV_xml):

    ''' retrive all the condition reported in the xml file associated to a
    specific clinvar id from ClinVar database 

    The function takes as input a VCV xml file from a ClinVar variant id and 
    return the dictionary of conditions associated to that specific variant id

    Parameters
    ----------
    clinvar_VCV_xml: xml file
            xml file from a VCV ClinVar entry 

    Returns
    -------
    conditions: dictionary 
       dictionary of conditions associated to a specific variant id
    '''
    classification_types = ['GermlineClassification',
                            'SomaticClinicalImpact',
                            'OncogenicityClassification']
    condition_out = {}
    condition_value =[]

    for classification_type in classification_types:
        try:               
            condition = clinvar_VCV_xml['ClinVarResult-Set']\
                                       ['VariationArchive']\
                                       ['ClassifiedRecord']\
                                       ['Classifications']\
                                       [classification_type]

            traits = condition['ConditionList']["TraitSet"]
            if not isinstance(traits, list):
                traits = [traits]
            for trait in traits:
                if isinstance(trait['Trait'],list):
                    trait_entity = trait['Trait'][0]['Name']
                else:
                    trait_entity = trait['Trait']['Name']

                if isinstance(trait_entity,list):
                    for element in trait_entity:
                        if element['ElementValue']["@Type"] == "Preferred":
                            condition_value.append(element['ElementValue']['#text'])
                else:
                    if trait_entity['ElementValue']["@Type"] == "Preferred":
                            condition_value.append(trait_entity['ElementValue']['#text'])

            condition_out[classification_type] = condition_value

        except:
            print(f"No condition associated to {classification_type} classification for variant {clinvar_id}")

    return condition_out


def classifications_extractor(clinvar_VCV_xml):
    '''
    Retrieve all the ClinVar classifications reported in the XML file 
    associated with a specific ClinVar ID from the ClinVar database.

    The function takes as input a VCV XML file from a ClinVar variant ID and
    returns the list of classifications associated with that specific variant ID.

    Parameters
    ----------
    clinvar_VCV_xml: file-like object
        XML file from a VCV ClinVar entry.

    Returns
    -------
    classifications: dictionary 
       dictionary of classifications associated with a specific ClinVar variant ID.
    '''

    ########################### Add description, review status and condition to the dictionary #####################
    # Description adding 
    classification_types = ['GermlineClassification',
                            'SomaticClinicalImpact',
                            'OncogenicityClassification']
    classifications = {}

    for classification_type in classification_types:
        try:               
            classification = clinvar_VCV_xml['ClinVarResult-Set']\
                                            ['VariationArchive']\
                                            ['ClassifiedRecord']\
                                            ['Classifications']\
                                            [classification_type]['Description']

            classifications[classification_type] = classification
        except:
           print(f"No {classification_type} classification for variant {clinvar_id}")
    return classifications

# add Review status information

def review_status_extractor(clinvar_VCV_xml):
    ''' retrive the ClinVar review_status in the input xml file
    associated to a specific clinvar id from ClinVar database

    The function takes as input a VCV xml file from a ClinVar variant id and
    return the list of conditions associated to that specific variant id

    Parameters
    ----------
    clinvar_VCV_xml: xml file
            xml file from a VCV ClinVar entry 

    Returns
    -------
    review_status: dictionary
       dictionary containing the review_status associated to a specific ClinVar variant id
    '''
    review_status = {}
    
    clinical_assertions = clinvar_VCV_xml['ClinVarResult-Set']\
                                         ['VariationArchive']\
                                         ['ClassifiedRecord']\
                                         ['Classifications']


    clinvar_id = clinvar_VCV_xml['ClinVarResult-Set']\
                                ['VariationArchive']\
                                ['ClassifiedRecord']\
                                ['SimpleAllele']\
                                ['@VariationID']
                                         

    for classification_type in clinical_assertions:
        try:
            review_status[classification_type] = clinical_assertions[classification_type]['ReviewStatus']
        except:
            print(f"no review status associated to {classification_type} classification for variant {clinvar_id}")



    return review_status

###############################################################################################################
#                                                                                                             #
#                                              MAIN SCRIPT                                                    #
#                                                                                                             #
###############################################################################################################

parser=argparse.ArgumentParser(description = 'ClinVar.py: \
                                              Retrieve information from ClinVar database\
                                              about condition, classification, review status \
                                              and methods for all missense variants for a specific\
                                              gene or specific missense variants of interest\
                                              depending on the provided input.')

group_input = parser.add_mutually_exclusive_group()
group_search = parser.add_mutually_exclusive_group()

group_input.add_argument('-g', "--gene_file",
                         dest = "gene_file", 
                         required = False, 
                         help = "csv \";\" separated input file " \
                                " containing the gene hugo name"\
                                " and isoform code (ref_seq)")

group_input.add_argument('-v',  "--variants_file",
                         dest = "variants_file", 
                         required = False,
                         help = "csv \";\" separated input file " \
                                " containing the missense_variants, "\
                                " the corresponding gene hugo name and "\
                                " its isoform code (ref_seq)")

group_search.add_argument('-c', "--cross_check",
                          action='store_true', 
                          dest = "cross_check", 
                          default = False,
                          help = "check if there is correspondence between variants annotated"\
                                " as missense variants in ClinVar database"\
                                " and missense variants annotated as Pathogenic,"\
                                " Benign, Likely Pathogenic, Likely Benign,"\
                                " Uncertain Significance and Conflicting interpretation "
                                "of Pathogenicity. This check is useful only when an input"\
                                "with only the gene and the isoform is provided.")

parser.add_argument("-m","--methods_inclusion",
                       dest = "methods", 
                       default = False, 
                       action = "store_true",
                        help = "choose to generate an output file including "\
                               " the method with which each variant"\
                               "has been classified")

parser.add_argument("-o","--output_csv",
                       dest = "out_csv", 
                       default = False, 
                       help = "output name")


args=parser.parse_args()

if __name__ == '__main__':
    ############################################################################################################################################
    #                                             ERROR PARSING + OBJECTS INITIALIZATION                                                       #
    ############################################################################################################################################

    if not args.variants_file and not args.gene_file: 
        print("The command line must contain the following entries: python3 [script_name.py]"\
              " [list_of_variant.csv] [output.csv].")

    if args.cross_check and args.variants_file:
        print("The cross_check mode of the variants in ClinVar database do not have sense when"\
              " a file with the mutations is provided. Please, change the command line to avoid"\
              " waiste of time and make the script faster. Exiting... ")
        exit(1)

    #------------------------------ Input Reading --------------------------#

  
    if args.gene_file:
        input_file = args.gene_file
    if args.variants_file:
        input_file = args.variants_file

    df = pd.read_csv(input_file, sep = ";")    # input csv file

    if len(df.columns.tolist()) != 3 and len(df.columns.tolist()) != 2:
        print(f"The input file must contain the following columns ';' separated depending on the analysis:\n"\
              f" if you're looking for all mutations information for a given gene the input must have the following columns\n\n"\
              f" gene;iso \n"\
              f" (gene_name);(isoform)\n\n"\
              f" if you're looking for information about specific mutations the input must have the following columns: \n\n"\
              f" gene;protein_var;iso\n"\
              f" (gene_name);(mutation);(isoform)\n\n")

        exit(1)

    #----------------------- Initial object assignment ---------------------#

    # objects common to gene and variants searches

    input_rows = df.to_numpy().tolist()   # convert input file into a numpy object with input rows as list
    not_found_var=[]                      # list of variants not found
    not_found_gene=[]                     # list of genes not found
    strange_variant_annotation={}         # dictionary of variants not canonically annotated, in which the key is 
                                          # the clinavr id and the value a list in which in the first position there is 
                                          # the mutation and in second position the gene name
    ids=[]                                # list with the Clinvar id associated to the mutations for a given gene
    
    # objects for gene search

    gene_name=[]                          # list with the genes name provided in the input file
    selected_isoform=[]                   # list with the isoform provided in the input file
    variants = []                         # list with missense mutations in coding region expressed in hgvs format
    filter_ids = {}                       # list with the Clinvar id associated to a specifica classification in ClinVar 
    inconsistent_annotations = {}         # list with the mutations that should be annotated in the missense class in ClinVar but actually not
    key_to_remove=[]                      # list of keys that will be removed form the output dictionary    
    missense_variants={}                  # dictionary to return in ehcih the key is the clinvar id and as value the
                                          # a list containing mutation gene name, classification, condition and review status

    # objects for variants search                             
    var_right_iso=0                       # count the misssense mutations in the isoform provided as input
    var_other_iso=0                       # count the missense mutations in other isoforms
    uncanonical_annotation=0              # count the missense mutations in the isoform provided as input with
                                          # a non-canonical annotation
    gene_search = False
    variant_search = False

    # distinguish which kind of research to perform

    if len(df.columns.tolist()) == 3:
        variant_search = True
        
    if len(df.columns.tolist()) == 2:
        gene_search = True

    ###########################################################################################################################################
    #                                     EXTRACTION OF CLINVAR IDS FROM MISSENSE MUTATIONS                                                   #
    ###########################################################################################################################################

    string = 50*" "+"Clinvar id extraction"+50*" "
    header = [string]
    table =[""]
    print(tabulate(table, header, tablefmt="fancy_grid"))
    
    #get access to the xml file in dictionary format for each gene
    for input_row in input_rows:
        gene = input_row[0]          # gene name from input file
        ref_seq = input_row[1]       # ref seq from input file
        if variant_search:
            variant = input_row[1]   # variant name in case a variant input file is provided
            ref_seq = input_row[2]   # ref_seq in case of a variant input file is provided
        i = 0                        # while loop stopping condition
        attempt = 0                  # starting attempt number
        max_attempts = 200           # ending attempt number
        delay = 3                    # delay between two consecutive queries

        while i == 0 and attempt < max_attempts:
            if attempt == max_attempts:
                print("ERROR: request failed after 200 attempts; exiting...")
                exit(1)  

            ########################## Get the clinvar ids associated to the missense mutations of interest ########################
            try:
                if gene_search:
                    term = "\"missense variant\""
                    filter_missense_variants="("+str(gene)+"%5Bgene%5D%20AND%20((%22missense%20variant%22%5Bmolecular%20consequence%5D%20OR%20%22SO%200001583%22%5Bmolecular%20consequence%5D)))"
                    clinvar_ids,out_gene = filtered_variants_extractor(filter_missense_variants,gene,ref_seq,"missense")
                    if clinvar_ids:
                        selected_isoform.append(ref_seq)
                        gene_name.append(gene)
                        ids.append(clinvar_ids)
                    else:
                        not_found_gene.append(out_gene)

                if variant_search:
                    URL_filter=str(gene)+"[gene]+AND+("+str(variant)+")[all]"
                    clinvar_ids, out_gene = filtered_variants_extractor(URL_filter, gene, ref_seq, "clinvar_id" )
                    if clinvar_ids:
                        selected_isoform.append(ref_seq)
                        gene_name.append(gene)
                        ids.append(clinvar_ids)
                        variants.append(variant)
                    else:
                        not_found_var.append(variant)
                        not_found_gene.append(out_gene)


                ################################### Get clinvar ids associated to the classified mutations #########################
                # URLs for selecting the Pathogenic, Likely Pathognic, Likely Benign and Benign ClinVar ids associated to the input gene
                if args.cross_check:
                    filter_pathogenic = "("+str(gene)+"%5Bgene%5D%20AND%20((%22clinsig%20pathogenic%22%5BProperties%5D%20or%20%22clinsig%20pathogenic%20low%20penetrance%22%5BProperties%5D%20or%20%22clinsig%20established%20risk%20allele%22%5BProperties%5D)"
                    filter_benign = "("+str(gene)+"%5Bgene%5D%20AND%20(%22clinsig%20benign%22%5BProperties%5D))"
                    filter_likely_pathogenic = "("+str(gene)+"%5Bgene%5D%20AND%20((%22clinsig%20likely%20pathogenic%22%5BProperties%5D%20or%20%22clinsig%20likely%20pathogenic%20low%20penetrance%22%5BProperties%5D%20or%20%22clinsig%20likely%20risk%20allele%22%5BProperties%5D)"
                    filter_likely_benign = "("+str(gene)+"%5Bgene%5D%20AND%20(%22clinsig%20benign%22%5BProperties%5D))%20AND%20(%22clinsig%20likely%20benign%22%5BProperties%5D))"
                    filter_conflicting = "("+str(gene)+"%5Bgene%5D%20AND%20(%22clinsig%20has%20conflicts%22%5BProperties%5D))"
                    filter_uncertain = "("+str(gene)+"%5Bgene%5D%20AND%20((%22clinsig%20vus%22%5BProperties%5D%20or%20%22clinsig%20uncertain%20risk%20allele%22%5BProperties%5D)))"

                    filters_on_classification = [filter_pathogenic, 
                                                 filter_benign, 
                                                 filter_likely_pathogenic, 
                                                 filter_likely_benign, 
                                                 filter_uncertain, 
                                                 filter_conflicting]

                    classifications = ["Pathogenic",     
                                       "Benign", 
                                       "Likely Pathogenic", 
                                       "Likely Benign", 
                                       "vus", 
                                       "Conflicting"]

                    for classification,URL in zip(classifications,filters_on_classification):
                        clinvar_ids,out_gene = filtered_variants_extractor(URL,gene,ref_seq,classification)
                        if clinvar_ids: 
                            filter_ids[classification] = [clinvar_ids, [gene], [ref_seq]]
                i=1
            except Exception:
                traceback.print_exc()
                print(f"WARNING: Request failed; will try again in {delay} seconds (attempt {attempt+1}/{max_attempts}")
                attempt+=1
                time.sleep(delay)

    ##########################    Dictionary  with clinvar_id as key and ref_seq and gene as values     ############################

    for index in range(len(ids)):
        for clinvar_id in range(len(ids[index])):
            missense_variants.update({ids[index][clinvar_id]: [selected_isoform[index],[gene_name[index]]]})
            if variant_search:
                 missense_variants.update({ids[index][clinvar_id]: [selected_isoform[index],[gene_name[index]],variants[index]]})


###############################################################################################################################################
#                                       ANNOTATION OF MUTATION INFORMATION FOR EACH CLINVAR_ID                                                #
###############################################################################################################################################

    ###### Parse each clinvar id to check if the mutation are associated to the main isoform and the mutations are annotated canonically ######

    string = 40*" "+"Missense variants information annotation"+41*" "
    header = [string]
    table =[""]
    print(tabulate(table, header, tablefmt="fancy_grid"))

    variants_to_check = missense_variants
    counter=0
    for clinvar_id in variants_to_check:
        gene = variants_to_check[clinvar_id][1][0]            # gene name from the entry in the dictionary
        isoform_to_check = variants_to_check[clinvar_id][0]   # isoform to check 
        variants_to_check[clinvar_id].pop(0)                  # remove the isoform information to provide a dictionary in the right format
        if variant_search:
            variant = variants_to_check[clinvar_id][1]
            variants_to_check[clinvar_id].pop(1) 
        correct_variant = ""                                
        counter+=1                  
        variation=[]  
        i=0
        attempt=0
        max_attempts=200
        delay=3
        while i==0 and attempt<max_attempts:
            if attempt==max_attempts:
                print("ERROR: request failed after 200 attempts; exiting...")
                exit(1)

            ##################################### get access to specific summary VCV xml file access ######################################### 

            try:
                VCV = ""
                print("Processing information for "+clinvar_id+" clinvar id from Clinvar. Progress --> "+str(counter)+"/"+str(len(variants_to_check)))

                ##################################### get access to specific summary VCV xml file access #####################################
                parse_VCV = VCV_summary_retriever(clinvar_id)

                ############################## get access to the hgvss list associated to the clinvar_id #####################################

                hgvss = coding_region_variants_extractor(parse_VCV,clinvar_id)

                if not hgvss:
                    if gene_search:
                        strange_variant_annotation.update({clinvar_id:["to check",variants_to_check[clinvar_id][0][0]]})
                    if variant_search:
                        strange_variant_annotation.update({clinvar_id:[variant,variants_to_check[clinvar_id][0][0]]})
                    key_to_remove.append(clinvar_id)

                else:

                    correct_variant = missense_variants_extractor(hgvss, gene, clinvar_id, isoform_to_check)

                    if correct_variant == None:
                        print("The clinvar_id "+clinvar_id+ " has a different annotation structure. It will be annotated in variants_to_check.csv")
                        strange_variant_annotation.update({clinvar_id:["to check",gene]})
                        key_to_remove.append(clinvar_id)


                    ################################################# retrieve classification ####################################################

                    classifications = classifications_extractor(parse_VCV)

                    ################################################### retrieve conditions ######################################################
                    
                    conditions = conditions_extractor(parse_VCV)

                    ############################################### retrieve assessment method ###################################################

                    methods = classification_methods_extractor(parse_VCV)

                    if isinstance(methods,tuple):
                        methods = methods[0]

                    ############################################### retrieve review status #######################################################
                    
                    review_status = review_status_extractor(parse_VCV)
                      
                    #the clinvar_id for which there are not mutations in the isoform provided as input will be added to the entry_not_found file
                    if gene_search:
                        if len(str(correct_variant))==0:
                            print(f"The {clinvar_id} clinavr_id for {gene} gene does not contain any mutation in the isoform provided as input.",\
                                  f" The entry will be added to entry_not_found.csv")
                            var_other_iso=var_other_iso+1
                            not_found_var.append(clinvar_id)
                            not_found_gene.append(gene)
                            key_to_remove.append(clinvar_id)
                        else:
                            var_right_iso=var_right_iso+1
                            if re.search("p.[A-Z][a-z][a-z][0-9]+[A-Z][a-z][a-z]",str(correct_variant)) or re.search("p.\w+delins\w+",str(correct_variant)):
                                variation.append(correct_variant)
                                variants_to_check[clinvar_id].insert(0,variation)
                                for feature in [classifications,conditions,review_status,methods]:
                                    variants_to_check[clinvar_id].append([feature])
                            else:
                                if not correct_variant == None:
                                    uncanonical_annotation=uncanonical_annotation+1
                                    strange_variant_annotation.update({clinvar_id:[correct_variant,gene]})
                                    key_to_remove.append(clinvar_id)
                                
                    if variant_search:
                        if correct_variant:
                            reference_variant = mutations_convert(correct_variant,d,True)
                            # compile the dictionary with the mutation associated to the right isoform
                            if reference_variant == str(variant):
                                variation.append(reference_variant)
                                variants_to_check[clinvar_id].insert(0,variation)
                                for feature in [classifications,conditions,review_status,methods]:
                                    variants_to_check[clinvar_id].append([feature])

                            # if the isoform is correct but the input mutation is not in the hgvss dicitonary put the mutation 
                            # and the correpsonding gene in not found lists.
                            else:
                                strange_variant_annotation.update({clinvar_id:[variant,gene]})
                                key_to_remove.append(clinvar_id)
                        else:
                            key_to_remove.append(clinvar_id)


                i=1

            except Exception:
                traceback.print_exc()
                print(f"WARNING: Request failed; will try again in {delay} seconds (attempt {attempt+1}/{max_attempts}")
                attempt+=1
                time.sleep(delay)

    if gene_search:

        print(f"{var_right_iso}/{len(variants_to_check)} missense mutations belong to the input isoform of which {var_right_iso-uncanonical_annotation}"\
              f" with the canonical annotation and {uncanonical_annotation} with a different"\
              f" annotation, {var_other_iso}/{len(variants_to_check)} missense mutations do not belong to the input isoform")

    # remove entries in which the clinvar id is associated with a non-proper mutation (i.e p.Gly123Ter) 
    for i in key_to_remove:
        del variants_to_check[i]

    ###############################################################################################################################################
    #                                                      CLINVAR ANNOTATION CONSISTENCY CHECKING                                                #
    ###############################################################################################################################################

    if args.cross_check:
        string = 38*" "+"Consistency check between ClinVar annotations"+38*" "
        header = [string]
        table =[""]
        print(tabulate(table, header, tablefmt="fancy_grid"))
        misannotated = {}
        for classification,ids in filter_ids.items():
            counter+=1                  
            variation=[] 
            i=0
            attempt=0
            max_attempts=200
            delay=3
            while i==0 and attempt<max_attempts:
                try:
                    if attempt==max_attempts:
                        print("ERROR: request failed after 200 attempts; exiting...")
                        exit(1)
                    clinvar_ids = ids[0]
                    variants = []
                    for clinvar_id in clinvar_ids:
                        variants_set = []
                        parse_VCV = VCV_summary_retriever(clinvar_id)
                        hgvss = coding_region_variants_extractor(parse_VCV,clinvar_id)
                        if hgvss:
                            correct_variant = missense_variants_extractor(hgvss, ids[1][0], clinvar_id, ids[2][0])
                            if re.search("p.[A-Z][a-z][a-z][0-9]+[A-Z][a-z][a-z]",str(correct_variant)) and "Ter" not in correct_variant:
                                if not clinvar_id in variants_to_check.keys():
                                    variants.append([clinvar_id,ids[1][0]])
                                    classifications = classifications_extractor(parse_VCV)
                                    conditions = conditions_extractor(parse_VCV)
                                    methods = classification_methods_extractor(parse_VCV)
                                    if isinstance(methods,tuple):
                                        methods = methods[0]
                                    review_status = review_status_extractor(parse_VCV)
                                    value = [[correct_variant],[ids[1][0]]]
                                    for feature in [classifications,conditions,review_status,methods]:
                                            value.append([feature])
                                    inconsistent_annotations[clinvar_id] = value

                    misannotated[classification] = variants

                    i=1

                except Exception:
                    traceback.print_exc()
                    print(f"WARNING: Request failed; will try again in {delay} seconds (attempt {attempt+1}/{max_attempts}")
                    attempt+=1
                    time.sleep(delay)

        for classification,variants in misannotated.items():
            if variants:
                for clinvar_id in variants:
                    print(f"annotation inconsistency detected for {clinvar_id[0]} clinvar_id in {clinvar_id[1]} gene")
            else:
                print(f"{classification} mutations passed the consistency check")
 


###############################################################################################################################################
#                                                              OUTPUT COMPILING                                                               #
###############################################################################################################################################

if not variants_to_check:
    print("WARNING: none of the variants given as input are included in"\
          " Clinvar database or an isoform code "\
          "not corresponding to the main isoform has" \
          " been provided: an empty output file will be generated.")

    df = pd.DataFrame.from_dict(variants_to_check, 
                                orient='index',
                                columns=['variant_id',
                                         'variant_name', 
                                         'gene_name', 
                                         'interpretation', 
                                         'condition',
                                         'review_status'])
    df.to_csv(args.out_csv, sep=";")

    # Process the input dataframe:
    
else:
    dfs = {}
    inconsistency_output_list = []
    output_list = melting_dictionary(variants_to_check)
    if args.cross_check:
        if inconsistent_annotations:
            inconsistency_output_list = melting_dictionary(inconsistent_annotations)
    
    columns = ['variant_id',
               'variant_name', 
               'gene_name', 
               'interpretation', 
               'condition',
               'review_status']
    if args.methods:
        columns.append('methods')
   
# assign the the values to the corresponding columns
    df = pd.DataFrame(output_list, columns=columns)
    dfs[args.out_csv] = df
    if inconsistency_output_list:
        df2 = pd.DataFrame(inconsistency_output_list, columns=columns)
        dfs["inconsitency_annotations.csv"] = df2

    for file_name,df in dfs.items():
        df['variant_name']      = df['variant_name'].str.get(0)
        df['gene_name']         = df['gene_name'].str.get(0)
        df['interpretation']    = df['interpretation'].str.get(0)
        df['review_status']     = df['review_status'].str.get(0)
        df['number_of_stars']   = [stars[x] for x in df['review_status']]
        df["condition"]         = (df["condition"].apply(lambda x: ",".join(x) if isinstance(x, list) else x))
        df["mutation"]          = df.apply(lambda x: mutations_convert(x["variant_name"],d,False),axis=1)
        if args.methods:
            df['methods']       = df['methods'].str.get(0)
        

        if file_name == args.out_csv:
            mutation_dict           = dict(zip(df.gene_name,df.mutation))

            # create a file with only the mutation (one letter code)
            if pd.notnull(df["mutation"]).any():
                for key in set(mutation_dict.keys()):
                    df_mutations=df.loc[df["gene_name"] == key]
                    with open(f"{key}_mutation_list.txt","a") as f:
                        for mutations in df_mutations["mutation"].tolist():
                            if mutations is not None:
                                f.write(mutations+"\n")
                    f.close()

        # remove useless columns
        df=df.drop(["mutation"],axis=1)
        df=df.drop(["review_status"],axis=1)
        if args.methods:
            out_csv_no_method = file_name.split(".")[0]+"_methods"+".csv"
            df.to_csv(out_csv_no_method, sep=";", index = False) 
            df = df.drop(["methods"],axis=1)
        df.to_csv(file_name,sep=";", index = False)

# write the file with the mutations to check manually
if strange_variant_annotation:
    data_strange=strange_variant_annotation
    df_strange=pd.DataFrame.from_dict(data_strange, 
                                          orient="index",
                                          columns=['variant_name','gene_name'])

    df_strange.index.name = 'clinvar_id'
    df_strange.to_csv(f"variants_to_check.csv",sep=";")


# write the file with the the entry not found (only genes if the list with the not found variants is empty)

if len(not_found_gene) != 0:
    if len(not_found_var) != 0:
        pd.DataFrame(list(zip(not_found_gene,not_found_var)), 
                    columns=['gene_name',"variant_name"]).to_csv(f"entry_not_found.csv", 
                    sep=";", 
                    index=False)
   
    else:
        pd.DataFrame(list(not_found_gene), 
                     columns=['gene_name']).to_csv(f"entry_not_found.csv", 
                     sep=";", 
                     index=False)