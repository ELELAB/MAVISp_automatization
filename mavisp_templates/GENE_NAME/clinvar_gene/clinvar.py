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


three_one_letter_annotation = {'Cys': 'C', 'Asp': 'D', 'Ser': 'S', 'Gln': 'Q', 'Lys': 'K',
                               'Ile': 'I', 'Pro': 'P', 'Thr': 'T', 'Phe': 'F', 'Asn': 'N',
                               'Gly': 'G', 'His': 'H', 'Leu': 'L', 'Arg': 'R', 'Trp': 'W',
                               'Ala': 'A', 'Val':'V', 'Glu': 'E', 'Tyr': 'Y', 'Met': 'M'}

stars = {'practice guideline':'4', 
         'reviewed by expert panel':'3', 
         'criteria provided, multiple submitters, no conflicts':'2',
         'criteria provided, multiple submitters': 'NA',
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
        if match.group(3) in dic.keys():
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

def melting_dictionary(variants_annotation, add_method = False):
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
    classifications = ['GermlineClassification',
                       'SomaticClinicalImpact',
                       'OncogenicityClassification']

    for clinvar_id in variants_annotation.keys():
        variant_information = [clinvar_id,variants_annotation[clinvar_id]["variant"],variants_annotation[clinvar_id]["gene"]]
        for classification in classifications:
            classification_check = []
            if classification in variants_annotation[clinvar_id]["classifications"].keys():
                classification_check.append(variants_annotation[clinvar_id]["classifications"][classification])
            if classification in variants_annotation[clinvar_id]["conditions"].keys():
                classification_check.append(variants_annotation[clinvar_id]["conditions"][classification])
            if classification in variants_annotation[clinvar_id]["review_status"].keys():
                classification_check.append(variants_annotation[clinvar_id]["review_status"][classification])
            if add_method:
                classification_check.append(variants_annotation[clinvar_id]["methods"])
            if len(classification_check) == 3 or len(classification_check) == 4:
                variant_features = variant_information+classification_check
                output_list.append(variant_features)
            else:
                print(f"No information for {classification} associated to the clinvar id {clinvar_id}")

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
    while attempts < max_retries:
        try:
            response = requests.get(URL)
        except:
            print(f"WARNING: An error occurred during the ClinVar databse query."
                  f" Will try again in {delay} seconds (attempt {attempts+1}/{max_retries})")
            attempts +=1
            time.sleep(delay)

        if response and response.status_code == 200:
            return response
        elif response:
            print(f"WARNING: The query to the ClinVar database to access {error_message}"\
                  f" returned {response.status_code} as response."\
                  f" Will try again in {delay} seconds (attempt {attempts+1}/{max_retries})")
            attempts +=1
            time.sleep(delay)
        else:
            print("WARNING: No response received. Retrying...")
            attempts += 1
            time.sleep(delay)
        
    raise RuntimeError("ERROR: request failed after 200 attempts; exiting...")
    


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
    error_message=f"summary xml file for the variant_id {clinvar_id}"
    function = "VCV_summary_retriever"
    parse_summary=xmltodict.parse(URL_response_check(URL_summary,error_message,function).content)
    try:
        VCV=(parse_summary['eSummaryResult']\
                          ["DocumentSummarySet"]\
                          ['DocumentSummary']['accession'])
    except KeyError:
        raise KeyError("Error in parsing the summary XML file. Check the ClinVar summary XML structure. Exiting...")
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

    try:
        first_search["eSearchResult"]["Count"]
    except KeyError:
        raise KeyError("Error in parsing the XML file at the Count key level. Check the ClinVar XML structure")

    if first_search["eSearchResult"]["Count"] != "0":
        if "clinvar_id" in mutation_type:
            clinvar_ids = variant_ids_extractors(first_search)
        else:
            tot_variant=first_search["eSearchResult"]["Count"]
            error_message=f"the number of missense variants for {gene} gene"
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

        return None,gene



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
                                                 ['SimpleAllele']
    except KeyError:
        print("Error with 'SimpleAllele' key: the clinvar_id "+clinvar_code+ " could have a different annotation structure. It will be annotated in variants_to_check.csv")
        simple_allele_accession={}
        
    if "HGVSlist" in simple_allele_accession.keys():
        try:
            HGVS_accession = clinvar_VCV_xml['ClinVarResult-Set']\
                                            ['VariationArchive']\
                                            ['ClassifiedRecord']\
                                            ['SimpleAllele']\
                                            ["HGVSlist"]\
                                            ["HGVS"]
        except KeyError:
            print("Error with 'HGVS' key: the clinvar_id "+clinvar_code+ " could have a different annotation structure for the mutations. It will be annotated in variants_to_check.csv")
            return list()

        try:
            hgvss = [ x for x in HGVS_accession if x['@Type'] == 'coding' and "ProteinExpression" in x.keys()]
                
            # there are cases in which the field "InterpretedRecord" is "Included Record", so epeat the previous step with the 
            # new key

        except Exception:
            print("The clinvar_id "+clinvar_code+ " does not have any annotation referring to variant in an expressed protein coding region")
            return list()
    else:
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
    return "wrong isoform"

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
    
    try:
        ids = clinvar_VCV_xml["eSearchResult"]\
                             ["IdList"]\
                             ["Id"]
    except:
        raise KeyError("Error in parsing the XML file at the IdList key level. Check the ClinVar XML structure")

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

group_input = parser.add_mutually_exclusive_group(required = True)


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

parser.add_argument('-c', "--cross_check",
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
                       default = "output.csv", 
                       help = "output name")


args=parser.parse_args()

if __name__ == '__main__':
    ############################################################################################################################################
    #                                             ERROR PARSING + OBJECTS INITIALIZATION                                                       #
    ############################################################################################################################################

    if args.cross_check and args.variants_file:
        print("The cross_check mode to check the consistency between variants "
              "annotated as missense variants and variants classified as Pathogenic,"\
              " Benign, Likely Pathogenic, Likely Benign,"\
              " Uncertain Significance and Conflicting interpretation makes sense only when"\
              " a file with the gene is provided. Exiting... ")
        exit(1)

    #------------------------------ Input Reading --------------------------#

    gene_search = False 
    variant_search = False

    if args.gene_file:
        gene_search = True  
        input_file = args.gene_file
        required_cols = set(['gene', 'isoform'])
        columns_name = {'gene':'gene','isoform':'ref_seq'}

    if args.variants_file:
        variant_search = True
        input_file = args.variants_file
        required_cols = set(['gene','protein_var','iso'])
        columns_name = {'gene':'gene','protein_var':'variant','iso':'ref_seq'}

    df = pd.read_csv(input_file, sep = ";")   # input csv file
   
    df = df.astype(str)

    if not required_cols.issubset(set(df.columns.tolist())):
        print(f"The input file must contain the following columns ';' separated:\n"\
              f" {';'.join(columns_name)}\n\n")

        exit(1)    

    df = df.rename(columns = columns_name)

    #----------------------- Initial object assignment ---------------------#

    # objects common to gene and variants searches

    not_found_var = []                      # List to store variants that were not found.
    not_found_gene = []                     # List to store genes that were not found.
    strange_variant_annotation = {}         # Dictionary to store variants not canonically annotated. Key is the ClinVar ID, 
                                            # and the value is a list containing the mutation in the first position and 
                                            # the gene name in the second position.
    ids = []                                # List of ClinVar IDs associated with mutations for a given gene.
    
    # objects for gene search

    gene_name = []                          # List of gene names provided in the input file.
    selected_isoform = []                   # List of isoforms provided in the input file.
    variants = []                           # List of missense mutations in coding region expressed in HGVS format.
    filter_ids = {}                         # Dictionary with ClinVar IDs associated with a specific classification in ClinVar.
    inconsistent_annotations = {}           # Dictionary with mutations that should be annotated in the missense class in ClinVar but are not.
    key_to_remove=[]                        # List of keys that will be removed from the output dictionary.
    missense_variants={}                    # dictionary to return in ehcih the key is the clinvar id and as value the
                                            # Dictionary to return, where the key is the ClinVar ID and the value is a list 
                                            # containing mutation, gene name, classification, condition, and review status.

    ###########################################################################################################################################
    #                                     EXTRACTION OF CLINVAR IDS FROM MISSENSE MUTATIONS                                                   #
    ###########################################################################################################################################

    print(f"|{'-' * 81}|")
    print(f"|{'Clinvar id extraction' : ^81}|")
    print(f"|{'-' * 81}|")
    
    #get access to the xml file in dictionary format for each gene
    for _, row in df.iterrows():                                     

    ########################## Get the clinvar ids associated to the missense mutations of interest ########################
        if gene_search:
            filter_missense_variants = "("+row["gene"]+"%5Bgene%5D%20AND%20((%22missense%20variant%22%5Bmolecular%20consequence%5D%20OR%20%22SO%200001583%22%5Bmolecular%20consequence%5D)))"
            mutation_type = "missense"
        if variant_search:
            filter_missense_variants = row["gene"]+"[gene]+AND+("+row["variant"]+")[all]"
            mutation_type = "clinvar_id"

        try:
            clinvar_ids,out_gene = filtered_variants_extractor(filter_missense_variants,row["gene"],row["ref_seq"],mutation_type)
        except RuntimeError as e:
            print(e)
            exit(1)

        if clinvar_ids:
            selected_isoform.append(row["ref_seq"])
            gene_name.append(row["gene"])
            ids.append(clinvar_ids)
            if variant_search:
                variants.append(row["variant"])
        else:
            not_found_gene.append(out_gene)
            if variant_search:
                not_found_var.append(row["variant"])

        ################################### Get clinvar ids associated to the classified mutations #########################
        # URLs for selecting the Pathogenic, Likely Pathognic, Likely Benign and Benign ClinVar ids associated to the input gene
        if args.cross_check:
            filter_pathogenic = "("+row["gene"]+"%5Bgene%5D%20AND%20((%22clinsig%20pathogenic%22%5BProperties%5D%20or%20%22clinsig%20pathogenic%20low%20penetrance%22%5BProperties%5D%20or%20%22clinsig%20established%20risk%20allele%22%5BProperties%5D)"
            filter_benign = "("+row["gene"]+"%5Bgene%5D%20AND%20(%22clinsig%20benign%22%5BProperties%5D))"
            filter_likely_pathogenic = "("+row["gene"]+"%5Bgene%5D%20AND%20((%22clinsig%20likely%20pathogenic%22%5BProperties%5D%20or%20%22clinsig%20likely%20pathogenic%20low%20penetrance%22%5BProperties%5D%20or%20%22clinsig%20likely%20risk%20allele%22%5BProperties%5D)"
            filter_likely_benign = "("+row["gene"]+"%5Bgene%5D%20AND%20(%22clinsig%20benign%22%5BProperties%5D))%20AND%20(%22clinsig%20likely%20benign%22%5BProperties%5D))"
            filter_conflicting = "("+row["gene"]+"%5Bgene%5D%20AND%20(%22clinsig%20has%20conflicts%22%5BProperties%5D))"
            filter_uncertain = "("+row["gene"]+"%5Bgene%5D%20AND%20((%22clinsig%20vus%22%5BProperties%5D%20or%20%22clinsig%20uncertain%20risk%20allele%22%5BProperties%5D)))"

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
                
                try:
                    clinvar_ids,out_gene = filtered_variants_extractor(URL,row["gene"],row["ref_seq"],classification)
                except RuntimeError as e:
                    print(e)
                    exit(1)

                if clinvar_ids: 
                    filter_ids[classification] = [clinvar_ids, [row["gene"]], [row["ref_seq"]]]

    ##########################    Dictionary  with clinvar_id as key and ref_seq and gene as values     ############################

    for index in range(len(ids)):
        for clinvar_id in range(len(ids[index])):
            variant_information = {"gene":gene_name[index],
                                   "isoform":selected_isoform[index],
                                   }
            if variant_search:
                variant_information["variant"] = variants[index]
            missense_variants[ids[index][clinvar_id]] = variant_information


###############################################################################################################################################
#                                       ANNOTATION OF MUTATION INFORMATION FOR EACH CLINVAR_ID                                                #
###############################################################################################################################################

    ###### Parse each clinvar id to check if the mutation are associated to the main isoform and the mutations are annotated canonically ######

    print(f"|{'-' * 81}|")
    print(f"|{'Missense variants information annotation' : ^81}|")
    print(f"|{'-' * 81}|")

    counter=0

    # objects for variants search                             
    var_right_iso=0                       # count the misssense mutations in the isoform provided as input
    var_other_iso=0                       # count the missense mutations in other isoforms
    uncanonical_annotation=0              # count the missense mutations in the isoform provided as input with
                                          # a non-canonical annotation
    for clinvar_id in missense_variants:
        gene = missense_variants[clinvar_id]['gene']                     # remove the isoform information gene name from the entry in the dictionary
        isoform_to_check = missense_variants[clinvar_id].pop('isoform')  # isoform to check                          
        if variant_search:
            variant = missense_variants[clinvar_id].pop('variant')
        correct_variant = ""                                
        counter += 1

        ##################################### get access to specific summary VCV xml file access ######################################### 

        print("Processing information for "+clinvar_id+" clinvar id from Clinvar. Progress --> "+str(counter)+"/"+str(len(missense_variants)))

        ##################################### get access to specific summary VCV xml file access #####################################
        
        try:
            parse_VCV = VCV_summary_retriever(clinvar_id)
        except (RuntimeError,KeyError) as e:
            print(e)
            exit(1)

        ############################## get access to the hgvss list associated to the clinvar_id #####################################

        hgvss = coding_region_variants_extractor(parse_VCV,clinvar_id)

        if not hgvss:
            clinvar_id_features = {}
            clinvar_id_features["gene"] = gene
            if gene_search:
                clinvar_id_features["variant"] = "to check"
            if variant_search:
                clinvar_id_features["variant"] = variant
            strange_variant_annotation[clinvar_id] = clinvar_id_features

            key_to_remove.append(clinvar_id)
            continue

        correct_variant = missense_variants_extractor(hgvss, gene, clinvar_id, isoform_to_check)

        if correct_variant is None:
            print("The clinvar_id "+clinvar_id+ " has a different annotation structure. It will be annotated in variants_to_check.csv")
            uncanonical_annotation += 1
            clinvar_id_features = {}
            clinvar_id_features["variant"] = "to check"
            clinvar_id_features["gene"] = gene
            strange_variant_annotation[clinvar_id] = clinvar_id_features
            key_to_remove.append(clinvar_id)
            continue


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
            if correct_variant=="wrong isoform":
                print(f"The {clinvar_id} clinvar id for {gene} gene does not contain any mutation in the isoform provided as input.",\
                      f" The entry will be added to entry_not_found.csv")
                var_other_iso+=1
                not_found_var.append(clinvar_id)
                not_found_gene.append(gene)
                key_to_remove.append(clinvar_id)
                continue
            var_right_iso+=1
            if re.search("p.[A-Z][a-z][a-z][0-9]+[A-Z][a-z][a-z]",str(correct_variant)) or re.search("p.\w+delins\w+",str(correct_variant)):
                missense_variants[clinvar_id]['variant'] = correct_variant
                variant_features = {}
                for feature,key_name in zip([classifications,conditions,review_status,methods],['classifications','conditions','review_status','methods']):
                    variant_features[key_name] = feature
                missense_variants[clinvar_id].update(variant_features)

            else:
                clinvar_id_features = {}
                clinvar_id_features["variant"] = correct_variant
                clinvar_id_features["gene"] = gene
                uncanonical_annotation=uncanonical_annotation+1
                strange_variant_annotation[clinvar_id] = clinvar_id_features
                key_to_remove.append(clinvar_id)
                    
        if variant_search:
            if not correct_variant:
                key_to_remove.append(clinvar_id)
                continue 

            reference_variant = mutations_convert(correct_variant,three_one_letter_annotation,True)
            # compile the dictionary with the mutation associated to the right isoform
            if reference_variant == variant:
                missense_variants[clinvar_id]['variant'] = variant
                variant_features = {}
                for feature,key_name in zip([classifications,conditions,review_status,methods],['classifications','conditions','review_status','methods']):
                    variant_features[key_name] = feature
                missense_variants[clinvar_id].update(variant_features)
            else:
                # if the isoform is correct but the input mutation is not in the hgvss dicitonary put the mutation 
                # and the corresponding gene in not found lists.
                clinvar_id_features = {}
                clinvar_id_features["variant"] = variant
                clinvar_id_features["gene"] = gene
                strange_variant_annotation[clinvar_id] = clinvar_id_features
                key_to_remove.append(clinvar_id)

            
                
    if gene_search:

        print(f"{var_right_iso}/{len(missense_variants)} missense mutations belong to the input isoform of which {var_right_iso-uncanonical_annotation}"\
              f" with the canonical annotation and {uncanonical_annotation} with a different"\
              f" annotation, {var_other_iso}/{len(missense_variants)} missense mutations do not belong to the input isoform")

    # remove entries in which the clinvar id is associated with a non-proper mutation (i.e p.Gly123Ter) 
    for i in key_to_remove:
        missense_variants.pop(i)

    ###############################################################################################################################################
    #                                                      CLINVAR ANNOTATION CONSISTENCY CHECKING                                                #
    ###############################################################################################################################################

    if args.cross_check:
        print(f"|{'-' * 81}|")
        print(f"|{'Consistency check between ClinVar annotations' : ^81}|")
        print(f"|{'-' * 81}|")
        misannotated = {}

        for classification,ids in filter_ids.items():
            clinvar_ids = ids[0]
            variants = []
            for clinvar_id in clinvar_ids:
                variants_set = []
                try:
                    parse_VCV = VCV_summary_retriever(clinvar_id)
                except (RuntimeError,KeyError) as e:
                    print(e)
                    exit(1)

                hgvss = coding_region_variants_extractor(parse_VCV,clinvar_id)
                if hgvss and re.search("p.[A-Z][a-z][a-z][0-9]+[A-Z][a-z][a-z]",str(correct_variant)) and "Ter" not in correct_variant:
                    correct_variant = missense_variants_extractor(hgvss, ids[1][0], clinvar_id, ids[2][0])
                    if not clinvar_id in missense_variants.keys():
                        variants.append([clinvar_id,ids[1][0]])
                        classifications = classifications_extractor(parse_VCV)
                        conditions = conditions_extractor(parse_VCV)
                        methods = classification_methods_extractor(parse_VCV)
                        if isinstance(methods,tuple):
                            methods = methods[0]
                        review_status = review_status_extractor(parse_VCV)
                        value = {"variant":correct_variant,"gene":ids[1][0]}
                        for feature,key_name in zip([classifications,conditions,review_status,methods],["classifications","conditions","review_status","methods"]):
                                value[key_name] = feature
                        inconsistent_annotations[clinvar_id] = value

            misannotated[classification] = variants

        for classification,mut_annotations in misannotated.items():
            if mut_annotations:
                for clinvar_id in mut_annotations:
                    print(f"annotation inconsistency detected for {clinvar_id[0]} clinvar id in {clinvar_id[1]} gene")
            else:
                print(f"{classification} mutations passed the consistency check")
 


###############################################################################################################################################
#                                                              OUTPUT COMPILING                                                               #
###############################################################################################################################################

columns = ['variant_id',
           'variant_name', 
           'gene_name', 
           'interpretation', 
           'condition',
           'review_status']
if args.methods:
    columns.append('methods')

if not missense_variants:
    print("WARNING: none of the variants given as input are included in"\
          " Clinvar database or an isoform code "\
          "not corresponding to the main isoform has" \
          " been provided: an empty output file will be generated.")

    df = pd.DataFrame.from_dict(missense_variants, 
                                orient='index',
                                columns=columns)
    df.to_csv(args.out_csv, sep=";")

    # Process the input dataframe:
    
else:
    dfs = {}
    inconsistency_output_list = []
    output_list = melting_dictionary(missense_variants,args.methods)
    if args.cross_check and inconsistent_annotations:
        inconsistency_output_list = melting_dictionary(inconsistent_annotations,args.methods)
   
# assign the values to the corresponding columns
    df = pd.DataFrame(output_list, columns=columns)
    dfs[args.out_csv] = df
    if inconsistency_output_list:
        df2 = pd.DataFrame(inconsistency_output_list, columns=columns)
        dfs["inconsitency_annotations.csv"] = df2

    df["condition"] = df["condition"].apply(lambda x: sorted(x))
    if args.methods:
            df["methods"] = df["methods"].apply(lambda x: sorted(x))

    for file_name,df in dfs.items():
        df['variant_name']      = df['variant_name']
        df['gene_name']         = df['gene_name']
        df['interpretation']    = df['interpretation']
        df['number_of_stars']   = [stars[x] for x in df['review_status']]
        df["condition"]         = (df["condition"].apply(lambda x: ",".join(x) if isinstance(x, list) else x))
        df["mutation"]          = df.apply(lambda x: mutations_convert(x["variant_name"],three_one_letter_annotation,False),axis=1)
        if args.methods:
            df['methods']       = (df["methods"].apply(lambda x: ",".join(x) if isinstance(x, list) else x))
        

        if file_name == args.out_csv:
            mutation_dict           = dict(zip(df.gene_name,df.mutation))

            # create a file with only the mutation (one letter code)
            if pd.notnull(df["mutation"]).any():
                for key in set(mutation_dict.keys()):
                    df_mutations=df.loc[df["gene_name"] == key]
                    with open(f"{key}_mutation_list.txt","w") as f:
                        for mutations in df_mutations["mutation"].tolist():
                            if mutations is not None:
                                f.write(mutations+"\n")

        # remove useless columns
        df = df.drop(["mutation"],axis=1)
        df = df.drop(["review_status"],axis=1)
        if args.methods:
            out_csv_no_method = file_name.split(".")[0]+"_methods"+".csv"
            df.to_csv(out_csv_no_method, sep=";", index = False) 
            df = df.drop(["methods"],axis=1)
        df.to_csv(file_name,sep=";", index = False)

# write the file with the mutations to check manually

if strange_variant_annotation:
    data_strange={}
    for clinvar_id,information in strange_variant_annotation.items():
        data_strange[clinvar_id] = [information["variant"],information["gene"]]
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


