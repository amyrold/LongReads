#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 28 13:11:31 2023

@author: aaronmyrold
"""
import os
import pandas as pd
from Bio import SeqIO
import itertools
import numpy as np
import editdistance

folder_names = ('1_raw_data', '2_filtered_data', '3_output', '4_blast')
p_raw_data = folder_names[0]
p_filt_data = folder_names[1]
p_out = folder_names[2]
p_blast = folder_names[3]

# Create any missing directories
for i in folder_names:
    if not os.path.exists(i):
        os.makedirs(i)
        
#%% ** 0.4 - Define Edit Distance Functions
# Define Functions

# FUNCTION 1
#pass in the comb_list to a function, which will populate the edit distance matrix
def ed(seq_pair):
    seq1_header = seq_pair[0]
    seq2_header = seq_pair[1]
    
    seq1 = records[seq1_header].seq
    seq2 = records[seq2_header].seq
    
    ed_df[seq1_header][seq2_header] = editdistance.eval(seq1,seq2)

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

#%% 3 - Edit Distance
#%% ** 3.1 - Calculating Edit Distances
print('Calculating Edit Distances')
records = SeqIO.to_dict(SeqIO.parse(f"{p_filt_data}/trim.fasta", format = "fasta"))

headers = [] 
hs_combo = []
for key,value in records.items(): #append id to headers list
    headers.append(key)
    hs_combo.append([key,value.seq])
#initializing matrix
ed_df = pd.DataFrame(np.zeros((len(headers),len(headers))),
              columns=headers,
              index=headers)

#comb list will be populated with all the different pairings of accession nums
comb_list = []
test_combo = []
for k in itertools.combinations(headers,2):
    comb_list.append(k)

for k in itertools.combination(hs_combo, 2):
    test_combo.append(k)
    
# How do we turn test_combo into dataframe OR how do we add seq1 seq2 cols to comb_list
test_df = pd.dataframe(test_combo)

#%% ** 3.? - Test Edit Distance


test_df['ed'] = editdistance.eval(test_df.loc('seq1'),test_df.loc('seq2'))

























    
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
#these csvs will be used in R for hierarchical clustering
ed_df.to_csv(f"{p_out}/matrix.csv")
pd.Series(within_ed).to_csv(f"{p_out}/within.csv", header = None)
pd.Series(between_ed).to_csv(f"{p_out}/between.csv", header = None)

print('Edit Distance has been calculated')
# END Part 3 ----