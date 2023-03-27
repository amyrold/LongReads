#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 27 18:23:37 2023

@author: aaronmyrold
"""
import os


# Download the .tsv into ______ folder

organism = 'Escherichia coli'
fields = '--fields accession,assminfo-sequencing-tech'
query = f"datasets summary genome taxon '{organism}' --assembly-level ‘complete’  --as-json-lines | dataformat tsv genome {fields}"
os.system()

# import the .tsv

# convert to pandas table

# drop rows that are not long read (pacbio, nanopore, ???)

# Store the accession column as a list

# iterate or query full list?
# either way, use accessions to download the appropriate genomes

# store genomes in data_raw folder