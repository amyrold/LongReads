import os
import pandas as pd
from Bio import SeqIO
import itertools
import numpy as np
import editdistance

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
ed_df = pd.DataFrame(np.zeros((len(headers),len(headers))),
             columns=headers,
             index=headers)

#comb list will be populated with all the different pairings of accession nums
comb_list = []
for k in itertools.combinations(headers,2):
    comb_list.append(k)
    
#pass in the comb_list to a function, which will populate the edit distance matrix
def ed(seq_pair):
    seq1 = seq_pair[0]
    seq2 = seq_pair[1]
    
    ed_df[seq1][seq2] = editdistance.eval(seq1,seq2)
    
#populates the within_coords and between_coords lists, depending on if the pairs are within the same genome or different genomes
def within_between(seq_pair):
    seq1 = seq_pair[0]
    seq2 = seq_pair[1]
    
    if seq1[0:13] == seq2[0:13]:
        within_coords.append(seq_pair) #appending the within
    else:
        between_coords.append(seq_pair) #appending the between
        
#get edit distance as iterating through the within coordinates
def get_ed_within(seq_pair):
    seq1 = seq_pair[0]
    seq2 = seq_pair[1]
    
    ed = ed_df[seq1][seq2] 
    
    if seq_pair in within_coords:
        within_ed.append(ed)
    
#get edit distance as iterating through the between coordinates
def get_ed_between(seq_pair):
    seq1 = seq_pair[0]
    seq2 = seq_pair[1]
    
    ed = ed_df[seq1][seq2] 
    
    if seq_pair in between_coords:
        between_ed.append(ed)
    
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

#sending the matrix, between_ed list, and within_ed list to csvs. 
#these csvs will be used in R for hierarchical clustering
ed_df.to_csv("/Users/nirushanbhag/Downloads/LongReads/matrix.csv")
pd.Series(within_ed).to_csv("/Users/nirushanbhag/Downloads/LongReads/within.csv", header = None)
pd.Series(between_ed).to_csv("/Users/nirushanbhag/Downloads/LongReads/between.csv", header = None)
