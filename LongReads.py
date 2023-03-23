#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 22 18:55:52 2023

@author: aaronmyrold
"""

import os
import subprocess
from Bio import Entrez
from Bio import SeqIO
# from Bio.Blast import NCBIWWW

# determine the path of the directory this file is located in
# idea taken from here: https://www.pythonanywhere.com/forums/topic/13464/
my_env = os.path.join(os.path.dirname(__file__)) #comment out unless running as script
# my_env = '/Users/aaronmyrold/Desktop/PipelineProject/LongReads' #aaron
# my_env = '' #niru
# my_env = '' #japani
# my_env = '' #asad
# set the current working directory to that folder so that remaining paths can function properly
os.chdir(my_env)


# PART 0 ----
# Make directories and paths
folder_names = ('1_data_raw', '2_data_clean', '3_output', '4_data_test', '5_blast')
p_data_raw = folder_names[0]
p_data_clean = folder_names[1]
p_out = folder_names[2]
p_index = folder_names[3]
p_test = folder_names[4]
p_blast = folder_names[5]

# this is primarily due to github not saving empty directories
# need to be able to create folder structure from scratch if needed  
for i in folder_names:
    if not os.path.exists(i):
        os.makedirs(i)
        
