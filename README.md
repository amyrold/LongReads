# LongReads

## Description
LongReads is a bioinformatics application that will compare the variation of the 16S region within and between strains of a specified bacterial species. Based on a user-specified genome, LongReads will download metadata and filter to genomes that were sequenced using Long Read technology (i.e. PacBio and Oxford Nanopore). It then analyzes and visualizes the intra-genomic and inter-genomic variation of the different 16S regions seen in each genome. 

# Installation Methods

There are three directories currently provided. The first, LR_base is the version you should use if you'd like to handle installation of packages on your own. There may be some troublshooting involved due to different versions of the NCBI datasets and dataformat tools. The second, LR_docker, is a self-contained image of all tools and dependencies involved in our main.py script. It currently does not support the R script and there are directions below on how to navigate that. Finally, LR_wheeler is built specifically to work on the compbio server with packages installed as they are. There is also "test" data downloaded to decrease time in the data download portion of the script. Directions will be below as well. 

## LR_wheeler process
Here, we have included the metadata and wgs.fasta files to decrease time spent downloading data. These are for 10 genomes. The rest of the script should run as intended. There are two main steps for testing this pipeline. First, run
```
python3 main.py -s 'Escherichia coli' -n 10
```
and once that finishes, run 
```
Rscript LR_stats.R
```

### Output
Once the script has finished running, four directories should be created:
* 1_raw_data
  *  metadata.tsv  - tab separated file listing the accession no and the associated sequencing technology (before filtering for long reads)
  *  longreads.tsv - tab separated file listing the accession no and the associated sequencing technology (after filtering for long reads)
  *  wgs.fasta - a multi fasta file of genomes of the specified species (after filtering for long reads)
* 2_filtered_data
  *  trim.csv - a csv of the 16S rRNA blast results, only the copy number is appended to the accession number column
  *  trim.fasta - a multi fasta file of the 16S rRNA copies
* 3_output
  *  between.csv - a csv of edit distances from the comparison of 16S copies between two genomes
  *  within.csv - a csv of edit distances from the comparison of 16S copies within one genome
  *  matrix.csv - a csv of the edit distance resultant matrix
  *  output.txt - a text file of summary statistics from the LR_stats.R script
  *  copies_per_genome.png - a plot produced from LR_stats.R showing the distribution of genomes with certain number of 16S rRNA copies
  *  inter_intra_shared_unique.png - a plot produced from LR_stats.R showing the relative abundance of shared and unique 16S copies
  *  boxplots.png - a boxplot of the distribution of edit distances based on type of variation
  
* 4_blast
  *  16S.*** files making up the 16S database for BLAST
  *  myresults.csv - a csv of the 16S rRNA blast results

## LR_docker Process
Docker is a free and open platform for developing, deploying, and running software. Docker isolates your applications from your infrastructure, allowing you to deliver software quickly. On this platform, you can manage your infrastructure the same way you manage your applications. Docker lets you package and run an application within a container, which is a loosely isolated environment.

This Docker image contains all required dependencies for this pipeline, as well as the pipeline script (main.py). 

Before building the image, you must first download docker from https://www.docker.com for your appropriate opperating system. You will also need to create and connect a docker hub account in order for the dockerfile to pull the appropriate starting image. 
Once docker is installed and you are logged in, you can continue to the next step -- cloning our repo

1. To build the docker image, first clone the repository using either
```
git clone https://github.com/amyrold/LongReads
```
or download the repo using the web/app GUI.

2. Before we can build our image, we need to pull the ubuntu image that it is based on. To do this, run:
```
sudo docker pull ubuntu
```

3. Next, build the image using:
```
sudo docker build LongReads --tag longreads:latest
```
Here, we need to make sure that "LongReads" matches the name of the cloned repo. If you cloned with a different folder name than "LongReads" simply update the docker build command to reflect that. 


4. Finally, create a Docker container with the longreads image. Here, [container name] is any user-given name.
```
sudo docker create -it --name [container name] longreads
sudo docker start -i [container name]
```
or start the interactive session automatically using
```
sudo docker run -it --name [container name] longreads
```

## LR_base Process
If docker is undesirable, we have also provided a directory with just the main script, pip requirements, and .R visualization script.
Here are the linux packages that are installed in the docker image. build-essential, ca-certs, and py-dev are all needed for biopython. They are likely installed already, but you can try installing them in case. The last two are required for the pipeline to run
```
apt-get install -y build-essential
apt-get install -y ca-certificates
apt-get install -y python3-dev
apt-get install -y ncbi-blast+
apt-get install -y ncbi-entrez-direct
```
Then, from within the LongReads/LR_base directory, run this command to download the required pip packages. 
```
pip install -r requirements.txt
```
Lastly, we need to install ncbi-datasets. For the script to run without changes, these commands must be run within the LongReads repo. 
```
curl -o datasets 'https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/linux-amd64/datasets'
curl -o dataformat 'https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/linux-amd64/dataformat'
chmod +x datasets dataformat
```


## Manual Installation
To run this program manually, you will need to download the following packages before continuing to executing the program
Here are the required dependancies to run main.py. They can all be installed via conda, pip, or apt-get. 
### Dependencies
- entrez-direct
- ncbi-datasets
- blast
- pandas
- numpy
- biopython
- editdistance

## Executing program
Once the container is up and running (or appropriate packages have been installed locally), the user can call the following command from within /LongReads/LR_base or /LongReads/LR_docker (depending on installation method) to begin the pipeline. If unspecified, the '-n' flag will download all Long Read genomes of specified bacteria. For testing the code, we reccomend using <50 genomes. 
```
python3 main.py -s [Bacterial Species of Choice] -n [Number of Genomes to download]
```
Once the main script has finished, the output files will be located in the 3_output directory in either /LongReads/LR_base or /LongReads/LR_docker (depending on installation method). You will then need to run the LR_stats.R script to create the necessary figures. Do this by running:
```
Rscript LR_stats.R
```
If you used the docker environment, we currently need to copy the 3_output directory back to your local machine before running the above command. To do this, run the following command:
```
sudo docker cp [container name]:/LongReads/3_output LR_docker
```
## Output
Once the script has finished running, four directories should be created:
* 1_raw_data
  *  metadata.tsv  - tab separated file listing the accession no and the associated sequencing technology (before filtering for long reads)
  *  longreads.tsv - tab separated file listing the accession no and the associated sequencing technology (after filtering for long reads)
  *  wgs.fasta - a multi fasta file of genomes of the specified species (after filtering for long reads)
* 2_filtered_data
  *  trim.csv - a csv of the 16S rRNA blast results, only the copy number is appended to the accession number column
  *  trim.fasta - a multi fasta file of the 16S rRNA copies
* 3_output
  *  between.csv - a csv of edit distances from the comparison of 16S copies between two genomes
  *  within.csv - a csv of edit distances from the comparison of 16S copies within one genome
  *  matrix.csv - a csv of the edit distance resultant matrix
  *  output.txt - a text file of summary statistics from the LR_stats.R script
  *  copies_per_genome.png - a plot produced from LR_stats.R showing the distribution of genomes with certain number of 16S rRNA copies
  *  inter_intra_shared_unique.png - a plot produced from LR_stats.R showing the relative abundance of shared and unique 16S copies
  *  boxplots.png - a boxplot of the distribution of edit distances based on type of variation
  
* 4_blast
  *  16S.*** files making up the 16S database for BLAST
  *  myresults.csv - a csv of the 16S rRNA blast results

## Authors
Aaron Myrold, Niru Shanbhag, Asad Shahzad, Japani Doan
