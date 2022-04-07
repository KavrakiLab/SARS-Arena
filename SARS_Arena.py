"""
SARS-Arena:  A Pipeline for Selection and Structural HLA Modeling of Conserved Peptides of SARS-related

## Installation

SARS-Arena is made available through [Docker Hub](https://hub.docker.com/r/kavrakilab/hla-arena), under the tag "sars-arena".

Installation instructions are also provided in our [github page](TODO).

## Modeller license key

SARS-Arena Workflow 2 rely on Modeller to perform the homology modeling of a given HLA receptor. This modeling task is integrated into a specific HLA-Arena function (more details below). However, using Modeller requires you to register and obtain your own license key, if you do not already have one. First, follow instructions on the [Modeller registration page](https://salilab.org/modeller/registration.html). 

Once you have the key, you can permanently update the SARS-Arena container with your key. For that, you should execute the commands below, replacing `MODELLER_KEY` with the correct key. 

    docker run -it kavrakilab/hla-arena:sars-arena
    sed -i "s/XXXX/MODELLER_KEY/g" /conda/envs/apegen/lib/modeller-9.20/modlib/modeller/config.py
    exit
    docker commit $(docker ps -a | sed -n 2p | awk '{ print $1 }') kavrakilab/hla-arena:sars-arena
    docker container rm $(docker ps -a | sed -n 2p | awk '{ print $1 }')

Note that this modification is permanent, in the sense that will not be lost when you close the container. However, it will be required again when you update the container (e.g., docker pull kavrakilab/hla-arena).

There is also the option of adding the key temporarily, when running a specific workflow. For that, just add the content below as one of the first cells to be executed in the workflow. Remember to replace `MODELLER_KEY` with the correct key.

    from subprocess import call
    call(["sed -i \"s/XXXX/" + MODELLER_KEY + "/g\" /conda/envs/apegen/lib/modeller-9.20/modlib/modeller/config.py"], shell=True)


## Using Jupyter Notebook

Each file with a '.ipynb' extension in this folder is a Jupyter notebook allowing you to run one of SARS-Arena workflows. Note that Jupyter Notebook is already installed in this docker image of SARS-Arena. If you are new to Jupyter Notebook, you can check this [tutorial on how to interact with its interface](https://www.youtube.com/watch?v=HW29067qVWk&feature=youtu.be&t=274). Numerous other resources are available online.


## Available workflows

[Workflow_1A.ipynb](http://127.0.0.1:8888/notebooks/ProjectDevelopment/Workflows/Peptide_Extraction_Workflow_1A.ipynb) (open link in a new tab)

Workflow 1A will allow that users run the multiple sequence alignment of SARS-CoV-2 proteins in loco. To avoid processing crashes, we recommend using this workflow to run no more than 50,000 protein sequences. Workflow 1A consists of five steps: (i) fetch dataset from NCBI, (ii) extract and filter the sequence file, (iii) Multiple Sequence Alignment, (iv) computing conservation score, and (v) computing conserved peptides.

[Workflow_1B.ipynb](http://127.0.0.1:8888/notebooks/ProjectDevelopment/Workflows/Peptide_Extraction_Workflow_1B.ipynb) (open link in a new tab)

Workflow 1B, differently from Workflow 1A, allows users to recover information from a pre-computed multiple sequence alignment based on a specific date. We recommend the use of this workflow for cases where there is a need to analyze a large number of protein sequences (e.g. more than 50.000 sequences). This workflow consists of three steps: (i) fetch Pre-computed MSA dataset, (ii) computing conservation score, (iii) computing conserved peptides.

[Workflow_1C.ipynb](http://127.0.0.1:8888/notebooks/ProjectDevelopment/Workflows/Peptide_Extraction_Workflow_1C.ipynb) (open link in a new tab)

The Workflow 1C, contrary to Workflows 1A and 1B, has a different purpose. Here, instead of analyzing only protein sequences from SARS-CoV-2, the user can also retrieve information from multiple sequence alignment with sequences from SARS-related coronaviruses. This workflow consists of four steps: (i) fetch Pre-computed Multiple Sequence Alignment (MSA) from SARS-CoV-2, (ii) multiple Sequence Alignment, (iii) computing conservation score, (iv) computing conserved peptides.

[Workflow_2.ipynb](http://127.0.0.1:8888/notebooks/ProjectDevelopment/Workflows/Peptide-HLA_Binding_Prediction_Workflow_2.ipynb) (open link in a new tab)

Workflow 2 provides a way to model the three-dimensional structure of selected peptides in the context of different HLAs. This workflow consists of five steps: (i) obtain peptide and HLA sequences for prediction, (ii) filter peptides using a sequence-based affinity prediction tool, (iii) model HLAs, (iv) model peptide-HLA complexes with APE-Gen, and (v) structural scoring functions.

"""

# -------------------------------------------------------------------------------- #
## Imports

# For utility functions used within the code
from utils.helper_funcs import *

# Subprocess module
from subprocess import check_output, STDOUT, PIPE, run, call, Popen
import multiprocessing
import time

#For printing purposes
import re
from itertools import groupby, zip_longest
from dateutil import parser
from datetime import datetime

#For the loading bar
from tqdm.notebook import tqdm

#For visualization and interaction purposes
import matplotlib.pyplot as plt
import seaborn as sns
from ipywidgets import *

#For data processing
import numpy as np
import pandas as pd
import math

#For consensus sequence
from Bio import AlignIO, SeqIO
from Bio.Align import AlignInfo

#For modelling HLA sequences
import HLA_Arena as arena

# For p3HLA scoring function
import pyrosetta
import pickle as pk


# -------------------------------------------------------------------------------- #
## WORKFLOW 1 functions

def call_ncbi_datasets(proteins, refseq_only, annotated_only, complete_only, host, released_since):
    
    """
    **Function**: Fetches protein sequence data from the NCBI datasets tool (parameter explanations taken from the [NCBI OpenAPI 3.0 REST API Docs](https://www.ncbi.nlm.nih.gov/datasets/docs/reference-docs/rest-api/))

    **Parameters**: <br />
    &ensp;&ensp;&ensp; - *proteins (str or list[str])*: Which proteins to retrieve in the data package. <br />
    &ensp;&ensp;&ensp; - *refseq_only (bool)*: If true, limit results to RefSeq genomes. <br />
    &ensp;&ensp;&ensp; - *annotated_only (bool)*: If true, limit results to annotated genomes. <br />
    &ensp;&ensp;&ensp; - *complete_only (bool)*: Only include complete genomes. <br />
    &ensp;&ensp;&ensp; - *host (str)*: If set, limit results to genomes extracted from this host (Taxonomy ID or name). <br />
    &ensp;&ensp;&ensp; - *released_since (str)*: If set, limit results to viral genomes that have been released after a specified date and time. April 1, 2020 midnight UTC should be formatted as follows: 2020-04-01T00:00:00.000Z. <br />

    **Returns**: <br />
    &ensp;&ensp;&ensp; - *protein_sequence_name (str)*: The name of the *.faa* protein sequence file downloaded from the NCBI datasets tool.

    """

    query_string = "https://api.ncbi.nlm.nih.gov/datasets/v1alpha/virus/taxon/sars2/protein/" + proteins \
                    + "/download?refseq_only=" + str(refseq_only) \
                    + "&annotated_only=" + str(annotated_only) \
                    + "&released_since=" + str(released_since) \
                    + "&host=" + host \
                    + "&complete_only=" + str(complete_only) \
                    + "&include_annotation_type=PROT_FASTA" \
    
    protein_sequence_name = fetch_file_and_unzip(query_string)

    return protein_sequence_name

def create_tab(workflow_dir):
    
    """
    **Function**: Creates a UI tab for selecting arguments, in order to fetch data from the NCBI Virus database

    **Parameters**: <br />
    &ensp;&ensp;&ensp; - *workflow_dir*: The workflow to download infor about pangolin lineage. <br />

    **Returns**: <br />
    &ensp;&ensp;&ensp; - *tab (widgets.Tab)*: An ipywidgets tab that apparears on the Workflow screen
    """
    
    Virus_type = widgets.Dropdown(options=['Sars-CoV', 'Sars-CoV-2', 'Both'],
                              value='Sars-CoV-2',
                              description='Virus:')

    Protein = widgets.Dropdown(options=['ORF1ab polyprotein', 'ORF1a polyprotein', 'leader protein', 'nsp2', 'nsp3',
                                        'nsp4', '3C-like proteinase', 'nsp6', 'nsp7', 'nsp8', 'nsp9', 'nsp10', 
                                        'RNA-dependent RNA polymerase', 'helicase', "3'-to-5' exonuclease",
                                        'endoRNAse', "2'-o-ribose methyltransferase", 'nsp11', 'surface glycoprotein',
                                        'ORF3a', 'envelope protein', 'membrane glycoprotein', 'ORF6', 'ORF7a',
                                        'ORF7b', 'ORF8', 'nucleocapsid phosphoprotein', 'ORF10'],
                               value='nucleocapsid phosphoprotein',
                               description='Protein:')

    Completeness = widgets.Dropdown(options=['Complete', 'Partial', 'Both'],
                                    value='Complete',
                                    description='Completeness:',
                                    style={'description_width': 'initial'})

    Host = widgets.Dropdown(options=['Human', 'All'],
                            value='Human',
                            description='Host:')
    
    RefSeq = widgets.Dropdown(options=['RefSeq', 'GenBank', 'Both'],
                              value='RefSeq',
                              description='Sequence Type:',
                              style={'description_width': 'initial'})

    general_accordion = widgets.Accordion(children=[Virus_type, Protein, Completeness, Host, RefSeq])
    general_accordion_titles = ['Virus', 'Protein', 'Completeness', 'Host', 'Sequence Type']
    for i, title in enumerate(general_accordion_titles):
        general_accordion.set_title(i, title)

    Isolation_source = widgets.SelectMultiple(options=['blood', 'feces', 'lung', 'lung, oronasopharynx',
                                                       'oronasopharynx', 'oronasopharynx, oronasopharynx',
                                                       'placenta', 'saliva, oronasopharynx', 'swab',
                                                       'urine'],
                                              value=[],
                                              description='Isolation source:',
                                              style={'description_width': 'initial'})

    Release_Date_From = widgets.DatePicker(
        description='From',
        value = datetime.today().replace(day=1).date(),
        disabled=False
    )

    Release_Date_To = widgets.DatePicker(
        description='To',
        value = datetime.today().date(),
        disabled=False
    )

    date_accordion = widgets.Accordion(children=[Release_Date_From, Release_Date_To])
    date_accordion.set_title(0, 'From')
    date_accordion.set_title(1, 'To')

    Geography = widgets.Dropdown(options=['Continent', 'Country', 'USA State'],
                                 value='Continent',
                                 description='Selection of geographic type:',
                                 style={'description_width': 'initial'})

    Continent = widgets.SelectMultiple(options=['Africa', 'Antartica', 'Asia', 'Europe', 'North America', 'Oceania', 
                                                'Oceans and Seas', 'South America'],
                                       value=[],
                                       description='Selection of continent:',
                                       style={'description_width': 'initial'})

    geography_accordion = widgets.VBox([Geography, widgets.HBox([Continent])])
    
    pangolin_storage = workflow_dir + '/lineage_notes.txt'
    fetch_pangolin_lineage(pangolin_storage)
    lineage_data = pd.read_csv(pangolin_storage, sep='\t', header=0)
    lineage_cleaning_query = (lineage_data["Lineage"].str.startswith('*')) | (lineage_data["Lineage"].str.startswith('X'))
    lineage_data = lineage_data[~lineage_cleaning_query]['Lineage']
    pangolin_lineage = widgets.SelectMultiple(options=lineage_data.values.tolist(),
                                              value=[],
                                              description='Pangolin Lineage:',
                                              style={'description_width': 'initial'})
    
    tab = widgets.Tab()
    tab_titles = ['General Information', 'Geographic Region', 'Isolation Source', 'Pangolin Lineage', 'Release Date']
    tab.children = [general_accordion, geography_accordion, Isolation_source, pangolin_lineage, date_accordion]
    for i, title in enumerate(tab_titles):
        tab.set_title(i, title)
    return tab
    
def dataset_selection(tab):
    
    """
    **Function**: Displays the UI ipywidgets tab for selecting arguments, in order to fetch data from the NCBI Virus database

    **Parameters**: <br />
    &ensp;&ensp;&ensp; - *tab*: The UI ipywidgets tab. <br />
    """
    
    display(tab)
    tab.children[1].children[0].observe(lambda x: update_country_tab(x, tab), names='value')

def call_ncbi_virus(Virus_Type, Protein, Completeness, Host, Refseq, Geographic_region, 
                         Isolation_source, Pangolin_lineage, Released_Dates):
    
    """
    **Function**: Fetches protein sequence data from the [NCBI Virus tool](https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/virus?SeqType_s=Protein))

    **Parameters**: <br />
    &ensp;&ensp;&ensp; - *Virus_Type (str)*: Which SARS virus to work with. <br />
    &ensp;&ensp;&ensp; - *Protein (str)*: Which proteins to retrieve in the data package. <br />
    &ensp;&ensp;&ensp; - *Completeness (str)*: Include complete/partial genomes. <br />
    &ensp;&ensp;&ensp; - *Host (str)*: If set, limit results to genomes extracted from this host. <br />
    &ensp;&ensp;&ensp; - *Refseq (str)*: Fetch only RefSeq/GenBank or both. <br />
    &ensp;&ensp;&ensp; - *Geographic_Region (str tuple)*: Tuple of Geographic type (Continent/Country/USA State) and the value <br />
    &ensp;&ensp;&ensp; - *Isolation_source (str or list[str])*: Different types of Isolation source <br />
    &ensp;&ensp;&ensp; - *Pangolin_lineage (str or list[str])*: Different types of Pangolin Lineages (Taken from [here](https://github.com/cov-lineages/pango-designation) <br />
    &ensp;&ensp;&ensp; - *Released_Dates (str tuple)*: Limit results to sequences that have been released in a specified date frame. <br />

    **Returns**: <br />
    &ensp;&ensp;&ensp; - *protein_sequence_name (str)*: The name of the *.fasta* protein sequence file downloaded from the NCBI Virus database.

    """
    
    query_string = 'q=*:*&fq={!tag=SeqType_s}SeqType_s:("Protein")'
    
    # Virus Type part
    if Virus_Type == "Sars-CoV":
        query_string += '&fq=VirusLineageId_ss:(694009)'
    elif Virus_Type == "Sars-CoV-2":
        query_string += '&fq=VirusLineageId_ss:(2697049)'
    elif Virus_Type == "Both":
        query_string += '&fq=VirusLineageId_ss:(2697049 OR 694009)'
    else: 
        return "Invalid Virus type, please set one of [Sars-CoV, Sars-CoV-2, Both] correctly!"   
    
    # Ambiguous parts removal
    query_string += '&fq={!tag=QualNum_i}QualNum_i:([0 TO 0])'
    
    # Protein part
    query_string += '&fq={!tag=ProtNames_ss}ProtNames_ss:("' + Protein + '")'
    
    # Sequence Type part
    if Refseq == "RefSeq":
        query_string += '&fq={!tag=SourceDB_s}SourceDB_s:("RefSeq")'
    elif Refseq == "GenBank":
        query_string += '&fq={!tag=SourceDB_s}SourceDB_s:("GenBank")'
    elif Refseq == "Both":
        pass
    else: 
        return "Invalid Sequence type, please set one of [RefSeq, GenBank, Both] correctly!" 
    
    # Completeness part
    if Completeness == "Complete":
        query_string += '&fq={!tag=Completeness_s}Completeness_s:("complete")'
    elif Completeness == "Partial":    
        query_string += '&fq={!tag=Completeness_s}Completeness_s:("partial")'
    elif Completeness == "Both": 
        pass
    else: 
        return "Invalid Completeness type, please set one of [Complete, Partial, Both] correctly!" 
    
    # Host part
    if Host == "Human":
        query_string += '&fq=HostLineageId_ss:(9606)'
    if Host != "Human" and Host != "All":
        return "Invalid Host, please set one of [Human, All]!"
    
    # Isolation source part
    list_length = len(Isolation_source)
    if list_length != 0:
        i = 0
        query_string += '&fq={!tag=Isolation_csv}Isolation_csv:('
        while i < list_length:
            query_string += '"' + Isolation_source[i] + '"'
            if i < list_length - 1:
                query_string += " OR "
            else:
                query_string += ")"
            i += 1
            
    # Pango lineage part
    list_length = len(Pangolin_lineage)
    if list_length >= 1025:
        return "Selection of pangolin lineages is too large, and the algorithm will fail! If you wish to choose all pangolin lineages, just leave the selection empty!"
    if list_length != 0:
        i = 0
        query_string += '&fq={!tag=Lineage_s}Lineage_s:('
        while i < list_length:
            query_string += '"' + Pangolin_lineage[i] + '"'
            if i < list_length - 1:
                query_string += " OR "
            else:
                query_string += ")"
            i += 1
    
    # Geographic region part
    if Geographic_region[0] == 'Continent':
        list_length = len(Geographic_region[1])
        if list_length != 0:
            i = 0
            query_string += '&fq={!tag=Region_s}Region_s:('
            while i < list_length:
                query_string += '"' + Geographic_region[1][i] + '"'
                if i < list_length - 1:
                    query_string += " OR "
                else:
                    query_string += ")"
                i += 1
    
    if Geographic_region[0] == 'Country':
        list_length = len(Geographic_region[1])
        if list_length != 0:
            i = 0
            query_string += '&fq={!tag=Country_s}Country_s:('
            while i < list_length:
                query_string += '"' + Geographic_region[1][i] + '"'
                if i < list_length - 1:
                    query_string += " OR "
                else:
                    query_string += ")"
                i += 1
            
    if Geographic_region[0] == 'USA State':
        list_length = len(Geographic_region[1])
        if list_length != 0:
            i = 0
            query_string += '&fq={!tag=USAState_s}USAState_s:('
            while i < list_length:
                query_string += '"' + Geographic_region[1][i] + '"'
                if i < list_length - 1:
                    query_string += " OR "
                else:
                    query_string += ")"
                i += 1
                
    # Date part
    try:
        date1 = datetime.combine(Released_Dates[0], datetime.min.time()).isoformat()
        date2 = datetime.combine(Released_Dates[1], datetime.min.time()).isoformat()
    except ValueError:
        return "Invalid isoformat dates"
    if(date1 > date2):
        return "Dates are reversed! Make sure the first is older than the newer one!"
    query_string += '&fq={!tag=CreateDate_dt}CreateDate_dt:([' + date1 + '.00Z TO ' + date2 + '.00Z])'
    
    # Final part
    query_string += '&cmd=download&sort=SourceDB_s desc,CreateDate_dt desc,id asc&dlfmt=fasta&fl=AccVer_s,Definition_s,Protein_seq'
    
    #print(query_string)
    
    protein_sequence_name = fetch_fasta_file(query_string)

    return protein_sequence_name


def count_sequences(protein_sequence_file):
    
    """
    **Function**: Prints and returns the total number of sequences in the *.faa* protein sequence file file downloaded from the NCBI datasets tool.

    **Parameters**: <br />
    &ensp;&ensp;&ensp; - *protein_sequence_name (str)*: The name of the *.faa* protein sequence file. <br />

    **Returns**: <br />
    &ensp;&ensp;&ensp; - *number_of_sequences (int)*: Number of sequences in the *.faa* file. 

    """ 
    
    command = "grep -c '>' " + protein_sequence_file
    results = check_output(command, stderr=STDOUT, shell=True)
    number_of_sequences = int(re.search(r'\d+', str(results)).group(0))
    print("Total number of sequences:", number_of_sequences)
    return number_of_sequences


def read_faa(input_file, output_file, N):

    """
    **Function**: Parses the *.faa* and filters sequences with unwanted characters. Additionally, only the first **N** sequences will be kept. 

    **Parameters**: <br />
    &ensp;&ensp;&ensp; - *input_file (str)*: The name of the *.faa* protein sequence file. <br />
    &ensp;&ensp;&ensp; - *output_file (str)*: The preferred name of the **filtered** *.faa* protein sequence file. <br />
    &ensp;&ensp;&ensp; - *N (int)*: The first **N** number of sequences that will be kept

    """ 
    
    i = 1
    with open(output_file, 'w') as output:
        for header,group in groupby(input_file, isheader):
            if header:
                line = next(group)
                ensembl_id = line
                if i == N + 1:
                    break
                i = i + 1
            else:
                temp_list = []
                X_flag = True
                for line in group:
                    if 'X' not in line:
                        temp_list.append(line)
                    else:
                        X_flag = False
                if X_flag:
                    output.write(ensembl_id)
                    for line in temp_list:
                        output.write(line)
    number_of_sequences = count_sequences(output_file)
    return str(N - number_of_sequences) + " sequences had invalid characters and were discarded"


def run_msa(protein_sequence_file, nthread, threshold):

    """
    **Function**: Runs multiple sequence alignment on the input sequences using [MAFFT](https://mafft.cbrc.jp/alignment/software/).

    **Parameters**: <br />
    &ensp;&ensp;&ensp; - *protein_sequence_name (str)*: The name of the *.faa* protein sequence file. <br />
    &ensp;&ensp;&ensp; - *nthread (int)*: The number of cores MAFFT will use to perform the alignment <br />
    &ensp;&ensp;&ensp; - *threshold (float)*: Threshold for calculating the consensus sequence of the alignment (frequencies below this threshold will have an unknown amino acid) <br />

    **Returns**: <br />
    &ensp;&ensp;&ensp; - *consensus_sequence (str)*: The consensus sequence obtained after the MSA. 

    """ 
    
    if(threshold > 1 or threshold < 0):
        return "Define a threshold between 0 and 1!"
    if(nthread < 1 or nthread > multiprocessing.cpu_count()):
        return "Define a proper cpu core number!"

    # Run alignment using MAFFT
    os.system("rm -f aligned.faa")
    command = ['mafft', '--auto', '--thread', str(nthread), protein_sequence_file]
    with open('aligned.faa', 'w') as f:
        call(command, stdout=f)
    
    #Store sequences in csv
    sequence_list = []
    for i, re in enumerate(SeqIO.parse('aligned.faa', 'fasta')):
        sequence_list.append((i + 1, str(re.seq)))
    pd.DataFrame(data=sequence_list).to_csv('aligned.csv')

    # Compute consensus sequence
    alignment = AlignIO.read('aligned.faa', 'fasta')
    summary_align = AlignInfo.SummaryInfo(alignment)

    return str(summary_align.gap_consensus(threshold))

def conservation_analysis(scoring_method, scoring_matrix):

    """
    **Function**: Computes a conservation score based on each position of the aligned sequences (source code is adopted from [here](https://compbio.cs.princeton.edu/conservation/)).

    **Parameters**: <br />
    &ensp;&ensp;&ensp; - *scoring_method (str)*: Preffered scoring method. Possible values are: <br />
        &ensp;&ensp;&ensp;&ensp;&ensp;&ensp; 1. `js_divergence` (Method by [Capra and Singh](https://academic.oup.com/bioinformatics/article/23/15/1875/203579)) <br />
        &ensp;&ensp;&ensp;&ensp;&ensp;&ensp; 2. `shannon_entropy` <br />
        &ensp;&ensp;&ensp;&ensp;&ensp;&ensp; 3. `property_entropy` (Method by [Mirny and Shakhnovich](https://www.sciencedirect.com/science/article/pii/S002228369992911X)) <br />
        &ensp;&ensp;&ensp;&ensp;&ensp;&ensp; 4. `vn_entropy` (Method by [Caffrey et al.](https://onlinelibrary.wiley.com/doi/full/10.1110/ps.03323604)) <br />
    &ensp;&ensp;&ensp; - *scoring_matrix (str)*: Preffered scoring matrix. This only applies to methods that actually use a scoring matrix for calculating conservation, like `js_divergence`, else, it is ignored (e.g. `shannon_entropy`). Possible values are: <br />
        &ensp;&ensp;&ensp;&ensp;&ensp;&ensp; 1. `blosum62` <br />
        &ensp;&ensp;&ensp;&ensp;&ensp;&ensp; 2. `blosum35` <br />
        &ensp;&ensp;&ensp;&ensp;&ensp;&ensp; 3. `blosum40` <br />
        &ensp;&ensp;&ensp;&ensp;&ensp;&ensp; 4. `blosum50` <br />
        &ensp;&ensp;&ensp;&ensp;&ensp;&ensp; 5. `blosum80` <br />
        &ensp;&ensp;&ensp;&ensp;&ensp;&ensp; 6. `blosum100` <br />

    **Returns**: <br />
    &ensp;&ensp;&ensp; - conservation_file (str): The name of the *.csv* conservation file.

    """ 
    
    scoring_matrix_path = '../conservation_code/matrix/' + scoring_matrix + '.bla'
    conservation_file = 'conservation.csv'
    os.system("rm -f " + conservation_file)
    command = ['python', '../conservation_code/score_conservation.py', '-m', scoring_matrix_path, 
               '-s', scoring_method, '-o', conservation_file, 'aligned.faa']
    result = run(command, stdout=PIPE, stderr=PIPE, universal_newlines=True)
    print(result.stdout, result.stderr)

    return conservation_file

def extract_peptides(min_len, max_len, aligned_sequences_df):
    """
    **Function**: Extracts all peptides of a given length of all the aligned sequences

    **Parameters**: <br />
    &ensp;&ensp;&ensp; - *min_len (int)*: The minimum peptide length. <br />
    &ensp;&ensp;&ensp; - *max_len (int)*: The maximum peptide length. <br />
    &ensp;&ensp;&ensp; - *aligned_sequences_df (pandas.DataFrame)*: Dataframe containing all the aligned sequences <br />

    **Returns**: <br />
    &ensp;&ensp;&ensp; - *extracted_peptides (list[str]): List of all the peptides extracted from all the aligned sequences.  

    """ 
    
    # Extract peptides of length `max_len`
    print("Extracting all peptides from sequences")
    Region_sequences = aligned_sequences_df['Aligned_Sequences'].tolist()
    peptide_list = []
    for sequence in tqdm(Region_sequences):
        peptides = [(i, i + max_len, max_len, sequence[i:i + max_len]) for i in range(len(sequence)- max_len + 1)]
        peptide_list.append(peptides)

    # Post-processing for all peptide lengths until `min_len`
    print("Post-processing for all peptide lengths")
    peptide_list = [peptide for peptide_sublist in peptide_list for peptide in peptide_sublist]
    peptide_list = sorted(list(set(peptide_list)), key=lambda element: (element[0], element[1]))
    extracted_peptides = []
    extracted_peptides.append(peptide_list)
    for pep_length in tqdm(range(min_len, max_len)):
        temp_list = peptide_list.copy()
        peptide_list = []
        for (start, stop, length, peptide) in temp_list:
            peptide_list.append((start, stop - 1, length - 1, peptide[0:(length - 1)]))
            peptide_list.append((start + 1, stop, length - 1, peptide[1:length]))
        peptide_list = sorted(list(set(peptide_list)), key=lambda element: (element[1], element[2]))    
        extracted_peptides.append(peptide_list)
    extracted_peptides = [peptide for final_sublist in extracted_peptides for peptide in final_sublist
                  if('-' not in peptide[3]) and ('X' not in peptide[3]) and ('J' not in peptide[3]) and ('B' not in peptide[3]) and ('Z' not in peptide[3])]
    extracted_peptides = sorted(list(set(extracted_peptides)), key=lambda element: (element[2], element[0], element[1]))

    return extracted_peptides


def interactive_plot_selection(conservation_df, extracted_peptides, min_len, max_len):
    """
    **Function**: Interactive plot of the conservation scores in each position. The user can interact with the conservation threshold, the rolling median window length and the peptide lengths to filter the desired peptides.

    **Parameters**: <br />
    &ensp;&ensp;&ensp; - *conservation_df (pandas.DataFrame)*: Dataframe that contains the alignment and the conservation score per sequence position. <br />
    &ensp;&ensp;&ensp; - *extracted_peptides (list[str])*: List of all the peptides extracted from all the aligned sequences. <br />
    &ensp;&ensp;&ensp; - *min_len (int)*: The minimum peptide length. <br />
    &ensp;&ensp;&ensp; - *max_len (int)*: The maximum peptide length. <br />

    """ 
    CV_slider = FloatSlider(value=conservation_df['Score'].mean(), min=round(conservation_df['Score'].min()) - 1, 
                      step=0.1, max=min(round(conservation_df['Score'].max()) + 1, 100), continuous_update=False)
    CV_cutoff = conservation_df['Score'].mean()

    RMW_slider = IntSlider(value=10, min=1, step=1, max=100, continuous_update=False)
    RMW_cutoff = 10

    Peptide_length_slider = IntSlider(value=min_len+1, min=min_len, step=1, max=max_len, continuous_update=False)
    Pep_length = min_len+1
    
    x = interact(handle_interact, CV_cutoff=CV_slider, RMW_cutoff=RMW_slider, Pep_length=Peptide_length_slider,
                 conservation_df=fixed(conservation_df), extracted_peptides=fixed(extracted_peptides))

    
def fetch_precomputed_sequences(year, month):

    """
    **Function**: Fetches the prealigned sequences from our repository. 

    **Parameters**: <br />
    &ensp;&ensp;&ensp; - *year (int)*: Sequences will be fetched from this year onwards. <br />
    &ensp;&ensp;&ensp; - *month (int)*: Sequences will be fetched from this month onwards. <br />

    **Returns**: <br />
    &ensp;&ensp;&ensp; - *protein_sequence_name (str)*: The name of the *.faa* protein sequence file that contains all requested aligned sequences.

    """
   
    query_string = "https://sars-arena.rice.edu:8000/get_msa/" + year + "/" + month

    protein_sequence_name = fetch_file_and_unzip(query_string)

    return protein_sequence_name   

# -------------------------------------------------------------------------------- #
## WORKFLOW 2 functions

def fetch_hla_sequences(hlas):
    """
    **Function**: Fetches the HLA sequences from the [IPD-IMGT/HLA Database](https://www.ebi.ac.uk/ipd/imgt/hla/)

    Parameters: <br />
    &ensp;&ensp;&ensp; - *hlas (list[str])*: List of HLAs.  It is important that the HLA name follows the pattern GENE*ALLELE GROUP:HLA PROTEIN (e.g. A*02:01, B*57:01, C*11:07) <br />

    Returns: <br />
    &ensp;&ensp;&ensp; - *hla_sequences (dict)*: dictionary where keys are HLAs and values are their sequences. 

    """         
    os.system("rm -f hla_prot.fasta")
    os.system("wget ftp://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/hla_prot.fasta")
    hla_reformatted = ['HLA-' + hla.replace("*", "").replace(":", "") for hla in hlas]
    hla_sequences = {}
    with open("hla_prot.fasta") as in_handle:
        for title, seq in SeqIO.FastaIO.SimpleFastaParser(in_handle):
            for i,hla in enumerate(hlas):
                if hlas[i] in title and hla_reformatted[i] not in hla_sequences.keys():
                    hla_sequences[hla_reformatted[i]] = seq
    return hla_sequences


def hla_filtering(hla_sequences):

    """
    **Function**: Filters the HLA sequences So that they are compatible with [MHCflurry](https://github.com/openvax/mhcflurry)

    Parameters: <br />
    &ensp;&ensp;&ensp; - *hla_sequences (dict)*: dictionary where keys are HLAs and values are their sequences. <br />

    Returns: <br />
    &ensp;&ensp;&ensp; - *hla_filtered_sequences (dict)*: dictionary where keys are filtered, MHCflurry-compatible HLAs and values are their sequences. <br />

    """   

    mhcflurry_alleles = list(pd.read_csv("../utils/MHCflurry_supported_alleles.txt", header=None)[0])
    hla_filtered_sequences = {}
    for key, value in hla_sequences.items():
        if key in mhcflurry_alleles:
            hla_filtered_sequences[key] = value

    return hla_filtered_sequences


def mhcflurry_scoring():
    """
    **Function**: Peptide-HLA pairs from Workflow 2 are being scored using [MHCflurry](https://github.com/openvax/mhcflurry)

    """     
    command = "mhcflurry-predict mhcflurry_input.csv --out predictions.csv"
    print("Calling MHCFlurry...")
    s = Popen(command, shell=True)
    s.wait()
    print("MHCflurry finished, collecting results...")
    f=IntProgress(min=0, max=100, description='MHCflurry:', bar_style='info')
    display(f)

    count = 0
    f.value = 0
    while count <= 100:
        time.sleep(.1)
        p1 = Popen(["cat", "predictions.csv"], stdout=PIPE)
        p2 = Popen(["wc", "-l"], stdin=p1.stdout, stdout=PIPE)
        val = int(p2.communicate()[0])
        if val == 0 and count < 80: val = 0.5
        count += float(val*100/40)
        f.value = count # signal to increment the progress bar

def mhcflurry_plot_selection(df_predictions, binder_cutoff):
    """
    **Function**: Interactive swarm plot of the binding affinity scores of peptide-HLA pairs. The x-axis denotes the different HLAs. The user can interact with the binding affinity threshold in order to filter the desired peptide-HLA pairs.

    **Parameters**: <br />
    &ensp;&ensp;&ensp; - *df_predictions (pandas.DataFrame)*: Dataframe with binding affinity predictions for each peptide-HLA pair. <br />
    &ensp;&ensp;&ensp; - *binder_cutoff (int)*: Binding affinity cutoff in order to filter peptide-HLA pairs (Default : 500nM) <br />

    """       
    max_thres = max(df_predictions['mhcflurry_prediction']) + 1000

    slider = IntSlider(value=500, min=0, step=50, max=int(max_thres), continuous_update=False)
    binder_cutoff.value = "500"

    x = interact(mhcflurry_handle_interact, df_predictions=fixed(df_predictions), cutoff=slider, binder_cutoff=fixed(binder_cutoff))

    
def model_hlas_MODELLER(hla_sequences):
    """
    **Function**: Calling MODELLER to model the input HLA sequences. 

    **Parameters**: <br />
    &ensp;&ensp;&ensp; - *hla_sequences (dict)*: dictionary where keys are HLAs and values are their sequences.  <br />

    """            
    hla_alleles = hla_sequences.keys()
    for hla_allele in hla_alleles:
        if os.path.exists(hla_allele + ".pdb"):
            print("Already found " + hla_allele + ".pdb")
            continue
        os.makedirs("./Modeller-files/" + hla_allele + "-modeller-output", exist_ok=True)
        os.chdir("./Modeller-files/" + hla_allele + "-modeller-output")
        with open("alpha_chain.fasta", "w") as f:
            f.write(">" + hla_allele)
            f.write("\n" + hla_sequences[hla_allele])
        f.close()
        arena.model_hla(("alpha_chain.fasta", hla_allele), num_models=2)
        call(["cp best_model.pdb ../../" + hla_allele + ".pdb"], shell=True)
        os.chdir("../../")

def model_structures(Filtered_peptides):
    """
    **Function**: Performs docking of peptide-HLA pairs with [APE-GEN](https://github.com/KavrakiLab/APE-Gen))

    **Parameters**: <br />
    &ensp;&ensp;&ensp; - *Filtered_peptides (pandas.DataFrame)*: Dataframe of peptide-HLA pairs to be modelled.  <br />

    **Returns**: <br />
    &ensp;&ensp;&ensp; - *best_scoring_confs (dict)*: dictionary where keys are peptide-HLA pairs and values are paths to the *.pdb* files that correspond to their best-modelled conformations. 

    """   
    best_scoring_confs = {}
    selected_hla_peptides = []
    list(Filtered_peptides.to_records(index=False))
    selected_hla_peptides = list(Filtered_peptides.to_records(index=False))
    i = 0
    fail_counter = 0
    while i < len(selected_hla_peptides):
        allele, peptide, mhcfpred = selected_hla_peptides[i]
        print ('-'*80)
        print("Running APE-Gen on HLA:" + allele +" peptide:"+ peptide)
        print ('-'*80)

        comp =  allele.replace("*", "")+"-"+peptide
        print(comp)
        root_dir = os.getcwd()
        call(["mkdir -p " + comp], shell=True)
        call(["cp " + allele+".pdb ./"+comp+"/"+allele+".pdb"], shell=True)
        os.chdir(comp)
        try:
            best_scoring_conf = arena.dock(peptide, "./"+allele+".pdb")
            print(best_scoring_conf)
            best_scoring_confs[comp.replace(":", "")] = best_scoring_conf
            os.chdir("../")
            i+=1
        except Exception as e:
            os.chdir(root_dir)
            call(["rm -r ./"+comp+"/"], shell=True)
            fail_counter +=1
            if fail_counter > 5:
                print("ERROR from APE-Gen generating structure more then five times "+comp)
                print("Skipping structure "+comp)
                fail_counter=0
                i+=1
            else:
                print("ERROR from APE-Gen generating structure "+comp)
                print("Repeating reconstruction of the structure "+comp)
    return best_scoring_confs

def score_structures(best_scoring_confs):
    """
    **Function**: Performs scoring of peptide-HLA pairs conformations with different scoring functions. 

    **Parameters**: <br />
    &ensp;&ensp;&ensp; - *best_scoring_confs (dict)*: dictionary where keys are peptide-HLA pairs and values are paths to the *.pdb* files that correspond to their best-modelled conformations. <br />

    **Returns**: <br />
    &ensp;&ensp;&ensp; - *scoring_results (pandas.DataFrame)*: DataFrame that contains the peptide-HLA pairs and the socres for each scoring function 

    """  
    init_rosetta()
    energies = {"Modeled HLAs": [], "peptide": [], "vina": [], "vinardo": [], "AutoDock4": [], "3pHLA": []}
    for key, value in best_scoring_confs.items():
        allele = key[4:9]
        peptide = key[10:]
        energy_vina = arena.rescore_complex_simple_smina(value, "vina")
        energy_vinardo = arena.rescore_complex_simple_smina(value, "vinardo")
        energy_ad4 = arena.rescore_complex_simple_smina(value, "ad4_scoring")
        energy_ppp = pyrosetta_ppp(allele, peptide, value)
        energies["Modeled HLAs"].append(allele)
        energies["peptide"].append(peptide)
        energies["vina"].append(energy_vina)
        energies["vinardo"].append(energy_vinardo)
        energies["AutoDock4"].append(energy_ad4)
        energies["3pHLA"].append(energy_ppp)
    scoring_results = pd.DataFrame(energies)
    return scoring_results

def energy_plot_selection(scoring_results, scoring_function, energy_cutoff):
    """
    **Function**: Interactive plot of the energy scores for each peptide-HLA pair. The user can interact with the energy cutoff so that wanted pHLA complexes are stored.

    **Parameters**: <br />
    &ensp;&ensp;&ensp; - *scoring_results (pandas.DataFrame)*: DataFrame that contains the peptide-HLA pairs and the scores for each scoring function. <br />
    &ensp;&ensp;&ensp; - *scoring_function (str)*: Scoring function to be chosen for energy filtering. Possible values are: <br />
    &ensp;&ensp;&ensp;&ensp;&ensp;&ensp; 1. `vina` <br />
    &ensp;&ensp;&ensp;&ensp;&ensp;&ensp; 2. `vinardo` <br />
    &ensp;&ensp;&ensp;&ensp;&ensp;&ensp; 3. `AutoDock4` <br />
    &ensp;&ensp;&ensp;&ensp;&ensp;&ensp; 4. `3pHLA` (In-house scoring method (link to paper/github pending)) <br />
    &ensp;&ensp;&ensp; - *energy_cutoff (str)*: Energy cutoff chosen by user for filtering pHLA complexes. Chosen by user (Default : Mean of the scoring function of the `scoring_results`)

    """      
    min_energy = np.min(scoring_results[scoring_function])
    max_energy = np.max(scoring_results[scoring_function])
    mean_energy = np.mean(scoring_results[scoring_function])
    scoring_specific_df = scoring_results[['Modeled HLAs', 'peptide', scoring_function]]
    slider = FloatSlider(min=min_energy-1, max=max_energy + 1, step=0.1, value = mean_energy, continuous_update=False)
    energy_cutoff.value = str(mean_energy)

    x = interact(energy_handle_interact, scoring_specific_df=fixed(scoring_specific_df), cutoff=slider,  min_energy=fixed(min_energy), energy_cutoff=fixed(energy_cutoff), scoring_function=fixed(scoring_function))

    

def store_best_structures(best_scoring_confs, selected, structures_storage_location):
    """
    **Function**: Stores the structures selected for further processing

    **Parameters**: <br />
    &ensp;&ensp;&ensp; - *best_scoring_confs (dict)*: dictionary where keys are peptide-HLA pairs and values are paths to the *.pdb* files that correspond to their best-modelled conformations. <br />
    &ensp;&ensp;&ensp; - *selected (dict)*: DataFrame that contains the filtered by the energy cutoff peptide-HLA pairs. <br />
    &ensp;&ensp;&ensp; - *structures_storage_location (str)*: Directory for the selected structures to be stored. 

    """ 
    for index, row in selected.iterrows():
        key = "HLA-" + row["Modeled HLAs"]+"-"+row["peptide"]
        path = best_scoring_confs[key]
        print("Writing "+key+" to " + structures_storage_location)
        call(["cp " + best_scoring_confs[key] + " " + structures_storage_location + "/" + key + ".pdb"], shell=True)
