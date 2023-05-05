# LongReads

## Description
LongReads is a bioinformatics application that will compare the variation of the 16S region within and between strains of a specified bacterial species. Based on a user-specified genome, LongReads will download metadata and filter to genomes that were sequenced using Long Read technology (i.e. PacBio and Oxford Nanopore). It then analyzes and visualizes the intra-genomic and inter-genomic variation of the different 16S regions seen in each genome. 

# Installation Methods

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
Then, from within the LongReads directory, run this command to download the required pip packages. 
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
- matplotlib

## Executing program
Once the container is up and running (or appropriate packages have been installed locally), the user can call the following command from within /LongReads to begin the pipeline. If unspecified, the '-n' flag will download all Long Read genomes of specified bacteria. For testing the code, we reccomend using <50 genomes. 
```
python3 main.py -s [Bacterial Species of Choice] -n [Number of Genomes to download]
```
Once the main script has finished, the output files will be located in the 3_output directory. You will then need to run the LR_stats.R script to create the necessary figures. Do this by running:
```
```
If you used the docker environment, we currently need to copy the 3_output directory back to your local machine before running the above command. To do this, run the following command:
```
sudo docker cp [container name]:/LongReads/3_output LR_docker
```


## Authors
Aaron Myrold, Niru Shanbhag, Asad Shahzad, Japani Doan
