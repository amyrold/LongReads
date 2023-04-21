#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 27 18:23:37 2023

@author: aaron
"""
#%%
# PART 0 ----
import os
import pandas as pd
from Bio import SeqIO

# Dynamically set cwd to ../LongReads
if os.getcwd()[-10:] == '/LongReads':
    my_env = os.getcwd()
elif os.getcwd()[-10:] == '/6_scripts':
    os.chdir('..')
    my_env = os.getcwd()
else:
    my_env = os.path.join(os.path.dirname(__file__))

# Make directories and paths
folder_names = ('1_raw_data', '2_filtered_data', '3_test_data', '4_output', '5_blast', '6_scripts')
p_raw_data = folder_names[0]
p_filt_data = folder_names[1]
p_test_data = folder_names[2]
p_out = folder_names[3]
p_blast = folder_names[4]
p_scripts = folder_names[5]

# Create any missing directories
for i in folder_names:
    if not os.path.exists(i):
        os.makedirs(i)
        
#%%
# PART 1 ----
# Define required functions

# Functions 1&2 are used to assign unique names/identifiers to each 16S copy within a single genome.
# The input is the BLASTn result (myresults.csv) which contains query accession, % identity, start/stop, and e-value
# It contains the data for each genome in our whole genome (wgs.fasta) multi-fasta file.
# The output of the two function is trim.csv -- myresults.csv with unique accessions for each 16S copy

# FUNCTION 1
# Define a function to find each unique accession and apply the rename_acc function to each accession
def split_acc(seq_df):
    # Find unique accessions
    unique = pd.DataFrame(seq_df['accession'].unique())
    # apply rename function to each accession
    unique.apply(rename_acc, axis=1)
    return

# FUNCTION 2
# Define function to add _# to each 16S copy within the genome
def rename_acc(row):
    #use the accession provided from the apply() (row[0]) to filter our dataframe
    # we also have to drop NAs and reset the index
    subset = pos16S.where(pos16S['accession'] == row[0]).dropna().reset_index(drop='True')
    # For each accession in the subset, replace the accession with an _#
    for i in range(1, len(subset['accession'])+1):
        subset['accession'][i-1] = f'{row[0]}_{i}'
    # append the data to our trim.csv without headers or index
    subset.to_csv(f'{p_filt_data}/trim.csv',mode='a',index=False,header=False)
    return 

# Functions 3&4 are used to to construct our 16S region (trim.fasta) multi-fasta file
# They take trim.csv and a list of SeqIO objects from reading in our wgs.fasta file
# Using the new accessions and start/stop positions, they store each 16S copy as a unique SeqIO object

# FUNCTION 3
# Function to use the BLASTn output to trim and store in new dict 
def trim_fa(accession, start, stop):
    #keep accession no, not the _#, because in the wgs fasta file, there aren't any copy numbers
    period_loc = accession.find('.')
    record = wgs_dict[accession[0:(period_loc+2)]][start:stop]
    #setting the record.id = to the full accession, which has the copy number appended to it
    record.id = accession
    # slice the unique accession to the wgs accession and store the correct seq to trim dict
    trim_dict[accession] = record
    return

# FUNCTION 4
# create function to send variables to trim_fa()
def store_16S(row):
    # take accession, start, stop and use to store the shortened sequence into new seq list
    trim_fa(row['accession'],row['start'], row['end'])
    return

#%%
# PART 2 ----
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
formatting = '10 qacc pident qstart qend length evalue'
# Call the BLASTn query
os.system(f'blastn -query {input_file} -db {p_blast}/16S -out {output_file} -outfmt "{formatting}"')

#%%
# PART 3 ----
# Append a unique identifier to the end of each accession to differentiate between copies of a genome
# Read in the BLASTn output as a pandas table
pos16S = pd.read_csv(f'{output_file}', names=['accession','pident','start','end', 'length','ev']) 
# read in the multi-fasta from data_download script
wgs_dict = SeqIO.to_dict(SeqIO.parse(f'{p_raw_data}/wgs.fasta', 'fasta')) 
# define dictionary to store trimmed sequences
trim_dict = {}
# Create trim.csv output file
trim_temp = pd.DataFrame(columns=['accession','pident','start','end','length', 'ev'])
trim_temp.to_csv(f'{p_filt_data}/trim.csv', index=False)

# Call the split function on our myresults.csv dataframe (pos16S)
split_acc(pos16S)

#%%
# PART 4 ---
# Using the newly named data, store our 16S copies into a new multifasta
# Read in the resulting dataframe with proper formatting
trim_16S = pd.read_csv(f'{p_filt_data}/trim.csv', dtype={'start':'Int32','end':'Int32'})
# remove any reads that are shorter than our cutoff 
trim_16S_f = trim_16S.mask(trim_16S['length'] < 1400).dropna()

# apply the store_16S function across all rows in out BLAST output pd table
trim_16S_f.apply(store_16S, axis=1)

#%%
# PART 5 ----
# Using the dictionary created in part 4, write the 16S reads with their unique identifier to our output (trim.fasta)
with open(f"{p_filt_data}/trim.fasta", "w") as output_handle:
    SeqIO.write(trim_dict.values(), output_handle, "fasta")

# END ----
#%%
