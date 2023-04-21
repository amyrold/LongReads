# LongReads

Simple overview of use/purpose.

## Description

An in-depth paragraph about your project and overview of use.

## Getting Started

### Dependencies
The primary requirement for running this pipeline is docker. 

### Installing

## Docker
Docker is a free and open platform for developing, deploying, and running software. Docker isolates your applications from your infrastructure, allowing you to deliver software quickly. On this platform, you can manage your infrastructure the same way you manage your applications. Docker lets you package and run an application within a container, which is a loosely isolated environment.

This Docker image contains all required dependencies for this pipeline, as well as the pipeline script (main.py). 

To build the docker image, first clone the repository using
```
git clone https://github.com/amyrold/LongReads
```

Next, build the image using
```
sudo docker build LongReads --tag longreads:latest
```

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


### Executing program
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
