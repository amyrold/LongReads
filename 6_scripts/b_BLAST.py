#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 27 18:23:37 2023

@author: aaronmyrold
"""

import os
import pandas as pd
from Bio import SeqIO
# from Bio.Blast import NCBIWWW

# determine the path of the directory this file is located in
# idea taken from here: https://www.pythonanywhere.com/forums/topic/13464/
my_env = os.path.join(os.path.dirname(__file__)) #comment out unless running as script
# my_env = '/Users/aaronmyrold/GitHub/LongReads/6_scripts' #aaron
# my_env = '' #niru
# my_env = '' #japani
# my_env = '' #asad
# set the current working directory to that folder so that remaining paths can function properly
os.chdir(my_env+ '/..')



# PART 0 ----
# Make directories and paths
folder_names = ('1_raw_data', '2_filtered_data', '3_test_data', '4_output', '5_blast', '6_scripts')
p_raw_data = folder_names[0]
p_filt_data = folder_names[1]
p_test_data = folder_names[2]
p_out = folder_names[3]
p_blast = folder_names[4]
p_scripts = folder_names[5]

# this is primarily due to github not saving empty directories
# need to be able to create folder structure from scratch if needed  
for i in folder_names:
    if not os.path.exists(i):
        os.makedirs(i)

# PART 1 ----
# BLAST wgs sequences against 16S BLAST database
# Create 16S BLAST database:
entrezQuery = "J01859" #determine query
os.system(f'esearch -db nucleotide -query "{entrezQuery}" | efetch -format fasta > {p_blast}/16s.fasta')

# Use reference sequences to create BLAST database
os.system(f'makeblastdb -in {p_blast}/16s.fasta -out {p_blast}/16S -title 16S -dbtype nucl')

# Define BLASTn query
input_file = f'{p_raw_data}/wgs.fasta'
output_file = f'{p_blast}/myresults.csv'
# using the formatting requested
formatting = '10 qacc pident qstart qend evalue'
# Call the BLASTn query
os.system(f'blastn -query {input_file} -db {p_blast}/16S -out {output_file} -outfmt "{formatting}"')


# PART 2 ----
# Trim wgs reads to just the regions determined via 16S BLAST
# Read in the BLASTn output as a pandas table
pos16S = pd.read_csv(f'{output_file}', names=['accession','pident','start','end','ev']) 
# read in the multi-fasta from data_download script
wgs_dict = SeqIO.to_dict(SeqIO.parse(f'{p_raw_data}/wgs.fasta', 'fasta')) 
# define dictionary to store trimmed sequences
trim_dict = {}

# Function to use the BLASTn output to trim and store in new dict 
def trim_fa(accession, start, stop):
    trim_dict[accession] = wgs_dict[accession][start:stop]
    return

# create function to send variables to trim_fa()
def store_16S(row):
    # take accession, start, stop and use to store the shortened sequence into new seq list
    trim_fa(row['accession'],row['start'], row['end'])
    return

# apply the store_16S function across all rows in out BLAST output pd table
pos16S.apply(store_16S, axis=1)

# PART 3 ----
# format and write to output
with open(f"{p_filt_data}/trim.fasta", "w") as output_handle:
    SeqIO.write(trim_dict.values(), output_handle, "fasta")




