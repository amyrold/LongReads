#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 27 18:23:37 2023

@author: aaronmyrold
"""
import os
import pandas as pd
from Bio import Entrez


# PART 0 ----
# Make directories and paths
folder_names = ('1_raw_data', '2_filtered_data', '3_test_data', '4_output', '5_blast', '6_scripts')
p_raw_data = folder_names[0]
p_filt_data = folder_names[1]
p_test_data = folder_names[2]
p_out = folder_names[3]
p_blast = folder_names[4]
p_scripts = folder_names[5]

cwd = os.getcwd()

os.chdir("..")
print(os.getcwd())

# Download the .tsv into ______ folder
organism = 'Escherichia coli'
fields = '--fields accession,assminfo-sequencing-tech'
query = f"datasets summary genome taxon '{organism}' --assembly-level 'complete'  --as-json-lines | dataformat tsv genome {fields}" + f" > {p_raw_data}/ecoli.tsv"
os.system(query)

tsv = pd.read_csv(f'{p_raw_data}/ecoli.tsv', delimiter= "\t") #replace the path with where the tsv file is

# drop rows that are not long read (pacbio, nanopore, ???)
f = open(f'{p_raw_data}/longreads.tsv', 'w') #replace this with the path to the longreads.tsv

for row_num in tsv.index:
    acc = str(tsv['Assembly Accession'][row_num])
    seq_tech = str(tsv['Assembly Sequencing Tech'][row_num])
    
    remove = [';', 'and', ',', '+', 'Illumina']
    
    if ('Nanopore' in seq_tech or 'PacBio' in seq_tech) and not any(x in seq_tech for x in remove):
        f.write(acc + "\t" + seq_tech + '\n')
        
f.close()
# Store the accession column as a list
acc = pd.read_csv(f'{p_raw_data}/longreads.tsv', sep = "\t", names = ['Accession', "Sequencing Technology"])
acc_list = list(acc['Accession'].values)

# Iterate through the accession list and store the genomes in a multi fasta file

for i in acc_list[0:5]:
    os.system(f'esearch -db nucleotide -query "{i}" | efetch -format fasta >> {p_raw_data}/wgs.fasta')




