#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 22 18:55:52 2023

@author: aaronmyrold
"""
#%%
# PART 0 - Import packages
import os

# Dynamically set cwd to ../LongReads
if os.getcwd()[-10:] == '/LongReads':
    my_env = os.getcwd()
elif os.getcwd()[-10:] == '/6_scripts':
    os.chdir('..')
    my_env = os.getcwd()
else:
    my_env = os.path.join(os.path.dirname(__file__))


# Create paths to each directory
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
# PART 1 - Global Variable / key arguments (if necessary)

# This section would be used to define any global variables
# Or to store any information that we pass in via program arguments
# For now, I do not think we need to worry about this. Let's get the wrapper script working first

# We will end up using this if we decide to make the pipeline function with any bacterial species

# In terms of arguments, it would be nice to figure out a way to have a "testing" mode
# general tutorial: https://linuxhint.com/add_command_line_arguments_to_a_python_script/
# Ideally, there would be a '-t' flag that would cause the pipeline to run using only a small subset of data
# We would need to create a global variable and pass that as an argument to data_download 
# so that we only download 10, 100, or 1000 reads
# i.e. the user could pass '-t' (default = 10 reads), '-t 100', '-t 1000'
# if the '-t' flag were not passed, then our default would be to download all reads
# basically, the argument from our wrapper could define a variable that would then be passed to the data download script as it's argument



#%%
# PART 2 - Data Download

# Really all we need to do is call the data download script and build in some confirmation that it completed successfully
# so, first, figure out how to call the script
# then have some kind of print statement at the end of data download that confirms it was successful
# i.e. print('data download successful')


#%%
#PART 3 - BLAST

# Really all we need to do is call the BLAST script and build in some confirmation that it completed successfully
# so, first, figure out how to call the script
# then have some kind of print statement at the end of BLAST that confirms it was successful
# i.e. print('BLAST successful')


#%%
# PART 4 - Edit Distance

# Really all we need to do is call the edit distance script and build in some confirmation that it completed successfully
# so, first, figure out how to call the script
# then have some kind of print statement at the end of edit distance that confirms it was successful
# i.e. print('edit distance successful')



#%%
# PART 5 - Statistics

# Really all we need to do is call the data stats script and build in some confirmation that it completed successfully
# so, first, figure out how to call the script
# then have some kind of print statement at the end of stats script that confirms it was successful
# i.e. print('completed successfully')


#%%
# PART 6 - Clustering

# Really all we need to do is call the data clustering script and build in some confirmation that it completed successfully
# so, first, figure out how to call the script
# then have some kind of print statement at the end of clustering script that confirms it was successful
# i.e. print('completed successfully')


# END ----
#%%