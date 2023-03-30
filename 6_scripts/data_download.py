#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 27 18:23:37 2023

@author: aaronmyrold
"""
import os
import pandas as pd
from Bio import Entrez


# Download the .tsv into ______ folder

organism = 'Escherichia coli'
fields = '--fields accession,assminfo-sequencing-tech'
query = f"datasets summary genome taxon '{organism}' --assembly-level ‘complete’  --as-json-lines | dataformat tsv genome {fields}"

# convert to pandas table
cwd = os.getcwd()
tsv = pd.read_csv(cwd + '/compbio.tsv', delimiter= "\t") #replace the path with where the tsv file is

# drop rows that are not long read (pacbio, nanopore, ???)
f = open(cwd + '/longreads.tsv', 'w') #replace this with the path to the longreads.tsv

for row_num in tsv.index:
    acc = str(tsv['Assembly Accession'][row_num])
    seq_tech = str(tsv['Assembly Sequencing Tech'][row_num])
    
    remove = [';', 'and', ',', '+', 'Illumina']
    
    if ('Nanopore' in seq_tech or 'PacBio' in seq_tech) and not any(x in seq_tech for x in remove):
        f.write(acc + "\t" + seq_tech + '\n')
        
f.close()
# Store the accession column as a list
acc = pd.read_csv(cwd + '/longreads.tsv', sep = "\t", names = ['Accession', "Sequencing Technology"])
acc_list = list(acc['Accession'].values)
acc_list

# iterate or query full list?

# either way, use accessions to download the appropriate genomes

# store genomes in data_raw folder
os.makedirs(cwd + '/data_raw')
os.chdir(cwd+ '/data_raw')

for k in acc_list:
    os.system("datasets download gene accession " + k)
    


