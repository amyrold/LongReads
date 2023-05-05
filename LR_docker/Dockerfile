FROM ubuntu
RUN echo 'APT::Install-Suggests "0";' >> /etc/apt/apt.conf.d/00-docker
RUN echo 'APT::Install-Recommends "0";' >> /etc/apt/apt.conf.d/00-docker
RUN DEBIAN_FRONTEND=noninteractive \
  apt update \
  && apt-get install -y build-essential \
  && apt-get install -y ca-certificates \
  && apt-get install -y curl \
  && apt-get install -y python3 \
  && apt-get install -y python3-dev \
  && apt-get install -y python3-pip \
  && apt-get install -y nano \
  && apt-get install -y ncbi-blast+ \
  && apt-get install -y ncbi-entrez-direct \
  && rm -rf /var/lib/apt/lists/*
WORKDIR /LongReads/
RUN curl -o datasets 'https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/linux-amd64/datasets'
RUN curl -o dataformat 'https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/linux-amd64/dataformat'
RUN chmod +x datasets dataformat
COPY requirements.txt /LongReads
RUN pip install -r requirements.txt
RUN rm requirements.txt
COPY main.py /LongReads/