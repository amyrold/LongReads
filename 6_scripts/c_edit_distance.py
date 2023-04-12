#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  7 13:50:49 2023

@author: aaronmyrold
"""

import os
import pandas as pd
from Bio import SeqIO
import itertools
import numpy as np
import editdistance
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
records = SeqIO.to_dict(SeqIO.parse("/Users/nirushanbhag/Downloads/LongReads/output.fasta", format = "fasta"))

headers = [] 
for key,value in records.items(): #append id to headers list
    headers.append(key)

#initializing matrix
df = pd.DataFrame(np.zeros((len(headers),len(headers))),
             columns=headers,
             index=headers)

comb_list = []
for k in itertools.combinations(headers,2):
    comb_list.append(k)
    
#pass in the comb_list to a function, which will populate the edit distance matrix
def ed(seq_pair):
    seq1 = seq_pair[0]
    seq2 = seq_pair[1]
    
    #print(editdistance.eval(seq1,seq2))
    
    df[seq1][seq2] = editdistance.eval(seq1,seq2)
    
def within_between(seq_pair):
    seq1 = seq_pair[0]
    seq2 = seq_pair[1]
    
    if seq1[0:13] == seq2[0:13]:
        within_coords.append(seq_pair)
    else:
        between_coords.append(seq_pair)
        
def get_ed_within(seq_pair):
    seq1 = seq_pair[0]
    seq2 = seq_pair[1]
    
    ed = df[seq1][seq2] 
    
    within_ed.append(ed)
    
def get_ed_between(seq_pair):
    seq1 = seq_pair[0]
    seq2 = seq_pair[1]
    
    ed = df[seq1][seq2] 
    
    between_ed.append(ed)
    
        
matrix_result = list(map(ed, comb_list))

within_coords = []
between_coords = []
within_between_result = list(map(within_between, comb_list))

within_ed = []
between_ed = []
get_ed_within_result = list(map(get_ed_within,within_coords))
get_ed_between_result = list(map(get_ed_between,between_coords))    



#%%
def ed(x, seq1, seq2):
    index = x.index.values
    seq1 = records[seq1].seq
    seq2 = records[seq2].seq
    
    return(editdistance.eval(seq1,seq2))

ed("NZ_CP007265.1_2", "NZ_CP007390.1_9")
    

    
df.apply(lambda x: x.index.values)