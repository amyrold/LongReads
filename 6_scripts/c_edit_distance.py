#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  7 13:50:49 2023

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