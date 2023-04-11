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




# records = list(SeqIO.parse(f'{p_raw_data}/wgs.fasta', format = "fasta"))

# trim_file = open(f"{p_filt_data}/trim.fasta", "w")

# df = pd.read_csv(f'{p_blast}/myresults.csv', names = ["qacc", "pident", "qstart", "qend", "evalue"])
# for k in range(1, len(df.qacc)): #for each entry in the results.csv
# 	qacc = df.qacc[k]
# 	qstart = df.qstart[k]
# 	qend = df.qend[k]
# 	for record in records:
# 		if qacc in record.id:
# 			trim_file.write(">" + str(record.description) + "\n")
# 			trim_file.write(str(record.seq[qstart:qend]) + "\n")



# PART 2 ----
# Trim wgs reads to just the regions determined via 16S BLAST
# Read in the BLASTn output as a pandas table
pos16S = pd.read_csv(f'{output_file}', names=['accession','pident','start','end','ev']) 
# read in the multi-fasta from data_download script
wgs_dict = SeqIO.to_dict(SeqIO.parse(f'{p_raw_data}/wgs.fasta', 'fasta')) 
# define dictionary to store trimmed sequences
trim_dict = {}
# Create trim.csv output file
trim_temp = pd.DataFrame(columns=['accession','pident','start','end','ev'])
trim_temp.to_csv('trim.csv', index=False)

# Define a function to find each unique accession and apply the rename_acc function to each accession
def split_acc(seq_df):
    # Find unique accessions
    unique = pd.DataFrame(seq_df['accession'].unique())
    # apply rename function to each accession
    unique.apply(rename_acc, axis=1)
    return

# Define function to add _# to each 16S copy within the genome
def rename_acc(row):
    #use the accession provided from the apply() (row[0]) to filter our dataframe
    # we also have to drop NAs and reset the index
    subset = pos16S.where(pos16S['accession'] == row[0]).dropna().reset_index(drop='True')
    # For each accession in the subset, replace the accession with an _#
    for i in range(1, len(subset['accession'])+1):
        subset['accession'][i-1] = f'{row[0]}_{i}'
    # append the data to our trim.csv without headers or index
    subset.to_csv('trim.csv',mode='a',index=False,header=False)
    return 

# Call the split function on our myresults.csv dataframe (pos16S)
split_acc(pos16S)

# Read in the resulting dataframe with proper formatting
trim_16S = pd.read_csv('trim.csv', dtype={'start':'Int32','end':'Int32'}) 

# Function to use the BLASTn output to trim and store in new dict 
def trim_fa(accession, start, stop):
    # slice the unique accession to the wgs accession and store the correct seq to trim dict
    trim_dict[accession] = wgs_dict[accession[:13]][start:stop]
    return

# create function to send variables to trim_fa()
def store_16S(row):
    # take accession, start, stop and use to store the shortened sequence into new seq list
    trim_fa(row['accession'],row['start'], row['end'])
    return

# apply the store_16S function across all rows in out BLAST output pd table
trim_16S.apply(store_16S, axis=1)


# PART 3 ----
# format and write to output
with open(f"{p_filt_data}/trim2.fasta", "w") as output_handle:
    SeqIO.write(trim_dict.values(), output_handle, "fasta")
