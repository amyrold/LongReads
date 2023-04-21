FROM tikhonovapolly/phigaro:latest
RUN bin/bash -c "/root/miniconda3/bin/conda install -c bioconda entrez-direct"
RUN bin/bash -c "/root/miniconda3/bin/conda install -c conda-forge ncbi-datasets-cli --yes"
RUN bin/bash -c "/root/miniconda3/bin/conda install -c bioconda blast --yes"
RUN bin/bash -c "/root/miniconda3/bin/conda install -c anaconda pandas"
RUN bin/bash -c "/root/miniconda3/bin/conda install -c anaconda numpy"
RUN bin/bash -c "/root/miniconda3/bin/conda install -c conda-forge biopython --yes"
RUN bin/bash -c "/root/miniconda3/bin/conda install -c conda-forge editdistance"
RUN bin/bash -c "/root/miniconda3/bin/conda install -c conda-forge matplotlib"
RUN apt-get update && apt-get install -y unzip vim
WORKDIR /LongReads
COPY main.py /LongReads/