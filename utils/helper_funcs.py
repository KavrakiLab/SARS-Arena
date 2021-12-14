# For fetching data from ncbi
from urllib.request import urlopen
from shutil import copyfileobj

# For protein file extraction
import os
import sys
from zipfile import ZipFile, is_zipfile

# Subprocess module
from subprocess import check_output, STDOUT, PIPE, run, call, Popen
import multiprocessing
import time

#For printing purposes
import re
from itertools import groupby, zip_longest

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

def fetch_file_and_unzip(query_string):
    
    with urlopen(query_string) as in_stream, open('/tmp/ncbi_dataset.zip', 'wb') as out_file:
        copyfileobj(in_stream, out_file)

    # Create a ZipFile Object and load the zipped dataset file in it
    with ZipFile('/tmp/ncbi_dataset.zip', 'r') as zipObj:

        # Get a list of all archived files from the zip
        listOfFileNames = zipObj.infolist()

        # Iterate over the files
        for fileName in listOfFileNames:
        
            # Check if filenime is a protein sequence file
            if fileName.filename.endswith('.faa'):
            
                # Extract this file from zip and store it in the currect directory (subject to change)
                fileName.filename = os.path.basename(fileName.filename)
                protein_sequence_name = os.path.basename(fileName.filename)
                zipObj.extract(fileName, os.getcwd())
            
            # Check if filenime is a aligned sequences file
            if fileName.filename.endswith('.csv'):
            
                # Extract this file from zip and store it in the currect directory (subject to change)
                fileName.filename = os.path.basename(fileName.filename)
                zipObj.extract(fileName, os.getcwd())
            
            # Check if filenime is a consensus sequence file
            if fileName.filename.endswith('.txt'):
            
                # Extract this file from zip and store it in the currect directory (subject to change)
                fileName.filename = os.path.basename(fileName.filename)
                zipObj.extract(fileName, os.getcwd())
    
    return protein_sequence_name

def isheader(line):
    return line[0] == '>'

def extract_conserved_peptides(CV_cutoff, RMW_cutoff, Pep_length, conservation_df, extracted_peptides):
    d, head = {}, None
    key_list = []
    conserved_positions = np.where(conservation_df['Score'].rolling(RMW_cutoff, center=True).median() > CV_cutoff)[0]
    for x in conserved_positions:
        if head is None or x != d[head].stop:
            head = x
            key_list.append(head)
        d[head] = range(head, x+1)
    
    limits_list = []
    for head in key_list:
        area_length = d[head].stop - head
        if area_length >= Pep_length:
            limits_list.append((head, d[head].stop, area_length))
    if len(limits_list) > 0:
        print("Conserved regions based on thresholds:")
        display(pd.DataFrame(data=limits_list, columns=['Start_Position', 'Stop_Position', 'Sequence_Length']))

    peptide_list = []
    for (start, stop, length) in limits_list:
        peptide_list.append([peptide for (beg, end, pep_len, peptide) in extracted_peptides 
                             if (beg >= start and end <= stop and pep_len == Pep_length)])
    peptide_list = [peptide for peptide_sublist in peptide_list for peptide in peptide_sublist]
    
    with open("peptides.list", "w") as outfile:
        outfile.write("\n".join(peptide_list))
    
    print("Number of unique peptides extracted:", len(peptide_list))
    return peptide_list

def plot_cutoff(ax, CV_cutoff, RMW_cutoff, conservation_df):
    ax.axhline(y=CV_cutoff, xmin=0, xmax=conservation_df.shape[0], color=sns.color_palette()[3])
    ax.axhspan(CV_cutoff, 100, facecolor=sns.color_palette()[3], alpha=0.1)
    ax.plot(conservation_df['Score'], marker='.', linestyle='-', alpha=0.3, linewidth=0.5, label='By Position')
    ax.plot(conservation_df['Score'].rolling(RMW_cutoff, center=True).median(),
    marker='.', linestyle='-', label=str(RMW_cutoff) + '-pos Rolling Mean')
    ax.set_xlabel('Conservation by Position')
    ax.set_ylabel('Conservation Value')
    ax.set_xlim(0, conservation_df.shape[0])
    ax.set_ylim(round(conservation_df[conservation_df['Score'] > .05]['Score'].min()) - 1, 
                min(round(conservation_df['Score'].max()) + 1, 100))
    ax.legend()

def handle_interact(CV_cutoff, RMW_cutoff, Pep_length, conservation_df, extracted_peptides):
    ax = plt.gca()
    plot_cutoff(ax, CV_cutoff, RMW_cutoff, conservation_df)
    plt.show()
    extract_conserved_peptides(CV_cutoff, RMW_cutoff, Pep_length, conservation_df, extracted_peptides)

def run_msa(protein_sequence_file, nthread, threshold):

    if(threshold > 1 or threshold < 0):
        return "Define a threshold between 0 and 1!"
    if(nthread < 1 or nthread > multiprocessing.cpu_count()):
        return "Define a proper cpu core number!"

    # Run alignment using MAFFT
    os.system("rm -f aligned.faa")
    command = ['mafft', '--auto', protein_sequence_file]
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

def helper_grouper(iterable, n, fillvalue=None):
    args = [iter(iterable)] * n
    return zip_longest(fillvalue=fillvalue, *args)

#------------------------------------
# WORKFLOW 2 FUNCTIONS:

def mhcf_pred(df_predictions, cutoff):
    hla_count = len(pd.unique(df_predictions['allele']))
    temp_df = df_predictions[df_predictions['mhcflurry_prediction'] < cutoff]
    print("No of selected peptides based on the cutoff:", temp_df.shape[0])
    display(temp_df[['allele', 'peptide']].groupby(['allele']).agg(['count']))
    print ('-'*80)
    print("Approximated time needed for modelling HLAs with MODELLER:")
    print("2 minutes per HLA")
    print("{:}".format(hla_count*2)+" min")
    print("{:.2f}".format(hla_count/30)+" h")
    print ('-'*80)
    print("Approximated time needed for modelling peptides with APE-Gen (8 cores per allele):")
    min_per_complex_APE_gen = 2
    print(str(min_per_complex_APE_gen)+" minutes per complex")
    print("{:}".format(temp_df.shape[0]*min_per_complex_APE_gen)+" min")
    print("{:.2f}".format(temp_df.shape[0]*min_per_complex_APE_gen/60)+" h")
    
#from utils import plot_data
def mhcflurry_plot_cutoff(ax, df_predictions, cutoff):
    ax.axhline(y = cutoff, xmin=0, xmax=50000, color=sns.color_palette()[3])
    ax.axhspan(0, cutoff, facecolor=sns.color_palette()[3], alpha=0.3)
    sns.swarmplot(x=df_predictions['allele'], y = df_predictions['mhcflurry_prediction'], ax=ax, data=df_predictions)

def mhcflurry_handle_interact(df_predictions, cutoff, binder_cutoff):
    ax = plt.gca()
    mhcflurry_plot_cutoff(ax, df_predictions, cutoff)
    mhcf_pred(df_predictions, cutoff)
    plt.show()
    binder_cutoff.value = str(cutoff)

# util functions for trating the data
def init_rosetta():
    pyrosetta.init()

def extract_ppp_energies(file_name, pep_len):
    scorefxn=pyrosetta.get_fa_scorefxn()
    #load rosetta
    pose = pyrosetta.pose_from_pdb(file_name)
    scorefxn(pose)
    #get energies
    res_ene = pose.energies().residue_total_energies_array()
    peptide_ene = res_ene[-pep_len:]
    #extract only relevant positions
    ppp_ene = []
    for i in [0, 1, 2, -2, -1]:
        for j in range(20):
            ppp_ene.append(peptide_ene[i][j])    
    return ppp_ene

def convert_nM(pred_val):
    return math.exp((1-pred_val)*math.log(50000))

def pyrosetta_ppp(allele, peptide, value):
    
    available_alleles = ['A0101', 'A0201', 'A0203', 'A0206', 'A0301', 'A1101', 'A2301',
                         'A2402', 'A2601', 'A2902', 'A3101', 'A6801', 'A6802', 'B0702',
                         'B0801', 'B1501', 'B1801', 'B2705', 'B3501', 'B3901', 'B4001',
                         'B4002', 'B4403', 'B5101', 'B5701', 'C0304', 'C0501', 'C1601']
    pep_len = len(peptide)

    if allele not in available_alleles:
        return None
    
    model_root = "../REGRmodels/"
    with open(model_root + allele + "_ppp_fulldata.pkl", 'rb') as pickle_file:
        model = pk.load(pickle_file)
        
    ene = extract_ppp_energies(value, pep_len)
    pred = convert_nM(model.predict([ene]))
    return pred

def energy_pred(cutoff, scoring_specific_df, scoring_function):
    print("Energy threshold: "+ str(cutoff) + " nM")
    print("Selected peptides:")
    display(scoring_specific_df[scoring_specific_df[scoring_function] < cutoff])

def energy_plot_cutoff(ax, energy_cutoff, min_energy):
    ax.axhline(y=energy_cutoff, xmin=0, xmax=15000, color=sns.color_palette()[3])
    ax.axhspan(min_energy-50, energy_cutoff, facecolor=sns.color_palette()[3], alpha=0.3)

def energy_handle_interact(scoring_specific_df, cutoff, min_energy, energy_cutoff, scoring_function):
    ax = plt.gca()
    energy_plot_cutoff(ax, cutoff, min_energy)
    energy_pred(cutoff, scoring_specific_df, scoring_function)
    g = sns.swarmplot(x="Modeled HLAs", y=scoring_function, ax=ax, data=scoring_specific_df)
    energy_cutoff.value = str(cutoff)
 