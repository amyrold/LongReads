FROM continuumio/miniconda3
RUN mkdir /root/LongReads
RUN bin/bash -c "/root/miniconda3/bin/conda install -c conda-forge biopython"
RUN bin/bash -c "/root/miniconda3/bin/conda install -c anaconda pandas"
RUN bin/bash -c "/root/miniconda3/bin/conda install -c anaconda numpy"
RUN bin/bash -c "/root/miniconda3/bin/conda install -c bioconda blast"
RUN bin/bash -c "/root/miniconda3/bin/conda install -c conda-forge ncbi-datasets-cli"
RUN bin/bash -c "/root/miniconda3/bin/conda install -c conda-forge editdistance"
RUN bin/bash -c "/root/miniconda3/bin/conda install -c conda-forge matplotlib"
RUN bin/bash -c "/root/miniconda3/bin/conda install -c bioconda entrez-direct
WORKDIR /LongReads
ADD main.py
