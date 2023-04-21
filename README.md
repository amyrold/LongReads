# LongReads

Simple overview of use/purpose.

## Description

An in-depth paragraph about your project and overview of use.

### Installation Methods

## Docker
Docker is a free and open platform for developing, deploying, and running software. Docker isolates your applications from your infrastructure, allowing you to deliver software quickly. On this platform, you can manage your infrastructure the same way you manage your applications. Docker lets you package and run an application within a container, which is a loosely isolated environment.

This Docker image contains all required dependencies for this pipeline, as well as the pipeline script (main.py). 

Before building the image, you must first download docker from https://www.docker.com for your appropriate opperating system. You will also need to create and connect a docker hub account in order for the dockerfile to pull the appropriate starting image. 
Once docker is installed and you are logged in, you can continue to the next step -- cloning our repo

To build the docker image, first clone the repository using
```
git clone https://github.com/amyrold/LongReads
```

Next, build the image using:
```
sudo docker build LongReads --tag longreads:latest
```
Here, we need to make sure that "LongReads" matches the name of the cloned repo. If you cloned with a different folder name than "LongReads" simply update the docker build command to reflect that. 


Finally, create a Docker container with the prophagedetective image using
```
sudo docker create -it --name [container name] longreads
sudo docker start -i [container name]
```
or start the interactive session automatically using
```
sudo docker run -it --name [container name] longreads
```
where 'container name' is any user-given name.

## Conda environemnt
If docker is undesirable, a conda env.yml file is provided.
First, create the conda environment
```
conda create env -f LongReads.yml
```
Then, activate the environment before continuing on to execution
```
conda activate LongReads
```

## Manual Installation
To run this program manually, you will need to download the following packages before continuing to executing the program
Here are the required dependancies to run main.py. They can all be installed via conda and some via pip.
### Dependancies
- entrez-direct
- ncbi-datasets-cli
- blast
- pandas
- numpy
- biopython
- editdistance
- matplotlib



## Executing program
Once the container is up and running, the user can call the following command from within /LongReads to begin the pipeline
```
python3 main.py
```

## Help

Any advise for common problems or issues.
```
command to run if program contains helper info
```

## Authors
Aaron Myrold
Niru Shanbhag
Asad Shahzad
Japani Doan
