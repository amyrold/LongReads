#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 20 18:07:04 2023

@author: aaronmyrold
"""

''' SCRIPT OUTLINE  '''
# 0 - Establish Variables and Functions
#  ** 0.1 - Import Packages
#  ** 0.2 - Create Directories
#  ** 0.3 - Define BLAST Functions
#  ** 0.4 - Define Edit Distance Functions
# 1 - Data Download
#  ** 1.1 - Download Metadata
#  ** 1.2 - Filter by Sequencing Technology
#  ** 1.3 - Download Long Read Genomes
# 2 - BLAST
#  ** 2.1 - Create BLAST Database
#  ** 2.2 - Make BLASTn Query
#  ** 2.3 - Append Unique Identifier to 16S Accessions
#  ** 2.4 - Tie 16S Sequences to Unique Accessions
#  ** 2.5 - Write the Multi-Fasta File
# 3 - Edit Distance
#  ** 3.1 - Calculating Edit Distance
#  ** 3.2 - Use ED to Create Matrix
#  ** 3.3 - Store Matrix and Vectors

#%% 0 - Establish Variables and Functions
#%% ** 0.1 - Import packages
import os
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
import itertools
import numpy as np
import editdistance
import argparse

# Set working directory
# my_env = './LongReads/'

parser = argparse.ArgumentParser(description='Quantify variation of the 16S region within and between bacterial strains of a specified species.')
parser.add_argument('--species','-s', type=str, help='The full species name of bacteria to be analyzed (i.e. Escherichia Coli)')
parser.add_argument('--num','-n', type=int, default=None, help='Determine the number of distinct genomes that should be downloaded. Leave blank to run analysis on all complete genomes associated with the specified species.')
args = parser.parse_args()

#%% ** 0.2 - Create Directories
# Create paths to each directory
folder_names = ('1_raw_data', '2_filtered_data', '3_output', '4_blast')
p_raw_data = folder_names[0]
p_filt_data = folder_names[1]
p_out = folder_names[2]
p_blast = folder_names[3]

# Create any missing directories
for i in folder_names:
    if not os.path.exists(i):
        os.makedirs(i)
        
#%% ** 0.3 - Define Part 1 Functions
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
    
    #remove overlap
    # unique.apply(drop_ol, axis=1)
    
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
def trim_fa(accession, start, stop, direction):
    period_loc = accession.find('.')
    record = wgs_dict[accession[:(period_loc+2)]][start:stop]
    record.id = accession
    if direction == 'plus':
        record.seq = record.seq.reverse_complement()
    # slice the unique accession to the wgs accession and store the correct seq to trim dict
    trim_dict[accession] = record
    return

# FUNCTION 4
# create function to send variables to trim_fa()
def store_16S(row):
    # take accession, start, stop and use to store the shortened sequence into new seq list
    trim_fa(row['accession'],row['start'], row['end'], row['direction'])
    return

#%% ** 0.4 - Define Edit Distance Functions
# Define Functions

# FUNCTION 1
#pass in the comb_list to a function, which will populate the edit distance matrix
def ed(seq_pair):
    seq1_header = seq_pair[0]
    seq2_header = seq_pair[1]
    
    seq1 = records[seq1_header].seq
    seq2 = records[seq2_header].seq
    
    ED = editdistance.eval(seq1,seq2)
    
    if ED > 700:
        seq1 = seq1.reverse_complement()
        ed_df[seq1_header][seq2_header] = editdistance.eval(seq1,seq2)
    else:
        ed_df[seq1_header][seq2_header] = ED

# FUNCTION 2    
#populates the within_coords and between_coords lists, depending on if the pairs are within the same genome or different genomes
def within_between(seq_pair):
    seq1 = seq_pair[0]
    seq2 = seq_pair[1]
    
    if seq1[0:13] == seq2[0:13]:
        within_coords.append(seq_pair) #appending the within
    else:
        between_coords.append(seq_pair) #appending the between

# FUNCTION 3        
#get edit distance as iterating through the within coordinates
def get_ed_within(seq_pair):
    seq1 = seq_pair[0]
    seq2 = seq_pair[1]
    
    ed = ed_df[seq1][seq2] 
    
    if seq_pair in within_coords:
        within_ed.append(ed)
 
# FUNCTION 4
#get edit distance as iterating through the between coordinates
def get_ed_between(seq_pair):
    seq1 = seq_pair[0]
    seq2 = seq_pair[1]
    
    ed = ed_df[seq1][seq2] 
    
    if seq_pair in between_coords:
        between_ed.append(ed)



#%% 1 - Data Download
#%% ** 1.1 - Download Metadata
# Download metadata
print('Downloading metadata...')
# Download the .tsv into ______ folder
organism = args.species
fields = '--fields accession,assminfo-sequencing-tech'
query = f"./datasets summary genome taxon '{organism}' --assembly-level 'complete'  --as-json-lines | ./dataformat tsv genome {fields}" + f" > {p_raw_data}/metadata.tsv"
os.system(query)

#%% ** 1.2 - Filter by Sequencing Technology
# read in the metadata and filter by sequencing technology
# We only want genomes that have been sequenced using long read sequencing

#replace the path with where the tsv file is
tsv = pd.read_csv(f'{p_raw_data}/metadata.tsv', delimiter= "\t") 
init_count = len(tsv)
print(f'There are {init_count} {organism} accessions')
# drop rows that are not long read (pacbio, nanopore, ???)
f = open(f'{p_raw_data}/longreads.tsv', 'w')

for row_num in tsv.index:
    acc = str(tsv['Assembly Accession'][row_num])
    seq_tech = str(tsv['Assembly Sequencing Tech'][row_num])
    
    remove = [';', 'and', ',', '+', 'Illumina']
    
    if ('Nanopore' in seq_tech or 'PacBio' in seq_tech) and not any(x in seq_tech for x in remove):
        f.write(acc + "\t" + seq_tech + '\n')

f.close()
#%% ** 1.3 - Download Long Read Genomes
# Read in the filtered metadata and download genomic fasta files into one multi-fasta
# Store the accession column as a list
acc = pd.read_csv(f'{p_raw_data}/longreads.tsv', sep = "\t", names = ['Accession', "Sequencing Technology"])
acc_list = list(acc['Accession'].values)
fin_count = len(acc_list)
print(f'There are {fin_count} {organism} long-read accessions')
num = args.num
if (num == None):
    # Iterate through the accession list and store the genomes in a multi fasta file
    for i in acc_list:
        os.system(f'esearch -db nucleotide -query "{i}" | efetch -format fasta >> {p_raw_data}/wgs.fasta')
else:
    # use the specified download number
    for i in acc_list[0:num]:
        os.system(f'esearch -db nucleotide -query "{i}" | efetch -format fasta >> {p_raw_data}/wgs.fasta')
print('Data Download is complete')
# END DATA DOWNLOAD----

#%% 2 - BLAST
#%% ** 2.1 - Create BLAST database
# BLAST wgs sequences against 16S BLAST database
print('Building BLAST database...')
# Create 16S BLAST database:
entrezQuery = "J01859" #determine query
os.system(f'esearch -db nucleotide -query "{entrezQuery}" | efetch -format fasta > {p_blast}/16s.fasta')

# Use reference sequences to create BLAST database
os.system(f'makeblastdb -in {p_blast}/16s.fasta -out {p_blast}/16S -title 16S -dbtype nucl')

#%% ** 2.2 - Make BLASTn Query
# Define BLASTn query
input_file = f'{p_raw_data}/wgs.fasta'
output_file = f'{p_blast}/myresults.csv'
# using the formatting requested
formatting = '10 qacc pident qstart qend length sstrand evalue'
# Call the BLASTn query
print('Running the BLASTn query...')
os.system(f'blastn -query {input_file} -db {p_blast}/16S -out {output_file} -outfmt "{formatting}"')
print('BLAST complete')

#%% ** 2.3 - Append Unique Identifier to 16S Accessions
# Append a unique identifier to the end of each accession to differentiate between copies of a genome
# Read in the BLASTn output as a pandas table
pos16S = pd.read_csv(f'{output_file}', names=['accession','pident','start','end', 'length', 'direction','ev']) 
copies_16S = len(pos16S)
print(f'There are {copies_16S} total copies of the 16S region')
# read in the multi-fasta from data_download script
wgs_dict = SeqIO.to_dict(SeqIO.parse(f'{p_raw_data}/wgs.fasta', 'fasta')) 
# define dictionary to store trimmed sequences
trim_dict = {}
# Create trim.csv output file
trim_temp = pd.DataFrame(columns=['accession','pident','start','end','length','direction', 'ev'])
trim_temp.to_csv(f'{p_filt_data}/trim.csv', index=False)

# Call the split function on our myresults.csv dataframe (pos16S)
split_acc(pos16S)

#%% ** 2.4 - Tie 16S Sequences to Unique Accessions
# Using the newly named data, store our 16S copies into a new multifasta
# Read in the resulting dataframe with proper formatting
trim_16S = pd.read_csv(f'{p_filt_data}/trim.csv', dtype={'start':'Int32','end':'Int32'})
# remove any reads that are shorter than our cutoff 
trim_16S_f = trim_16S.mask(trim_16S['length'] < 1400).dropna()

# apply the store_16S function across all rows in out BLAST output pd table
trim_16S_f.apply(store_16S, axis=1)

#%% ** 2.5 - Write the Multi-Fasta file
# Using the dictionary created in part 4, write the 16S reads with their unique identifier to our output (trim.fasta)
with open(f"{p_filt_data}/trim.fasta", "w") as output_handle:
    SeqIO.write(trim_dict.values(), output_handle, "fasta")

print('trim.fasta has been created')
# END PART 2 ----

#%% 3 - Edit Distance
#%% ** 3.1 - Calculating Edit Distances
print('Calculating Edit Distances')
records = SeqIO.to_dict(SeqIO.parse(f"{p_filt_data}/trim.fasta", format = "fasta"))

headers = [] 
for key,value in records.items(): #append id to headers list
    headers.append(key)

#initializing matrix
ed_df = pd.DataFrame(np.zeros((len(headers),len(headers))),
              columns=headers,
              index=headers)

#comb list will be populated with all the different pairings of accession nums (order doesn't matter: AB is the same as BA, so it will only pick one)
comb_list = []
for k in itertools.combinations(headers,2):
    comb_list.append(k)
    
#%% ** 3.2 - Use ED to Create Matrix 
#passing the ed function to the list of combinations list
matrix_result = list(map(ed, comb_list))

#initializing a list of within and between coordinates
within_coords = []
between_coords = []

#passing the within_between function to the combinations list
within_between_result = list(map(within_between, comb_list))

#initializing a list of within and between edit distances. We will later do some stats on these
within_ed = []
between_ed = []

#passing the get_ed_within function to the list of within coordinates
get_ed_within_result = list(map(get_ed_within,within_coords)) 

#passing the get_ed_between function to the list of between coordinates
get_ed_between_result = list(map(get_ed_between,between_coords))    

#%% ** 3.3 - Store the Matrix and Vectors
#sending the matrix, between_ed list, and within_ed list to csvs. 
#these csvs will be used in R 
ed_df.to_csv(f"{p_out}/matrix.csv")
pd.Series(within_ed).to_csv(f"{p_out}/within.csv", header = None)
pd.Series(between_ed).to_csv(f"{p_out}/between.csv", header = None)

print('Edit Distance has been calculated')
# END Part 3 ----
