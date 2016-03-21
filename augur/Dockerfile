# augur dockerfile

FROM ubuntu:14.04
MAINTAINER Trevor Bedford <trevor@bedford.io>
RUN apt-get -y update

# wget
RUN apt-get install -y wget

# git
RUN apt-get install -y git

# python
RUN apt-get install -y python python-dev python-pip python-virtualenv
RUN apt-get install -y python-numpy python-scipy 
RUN apt-get install -y libpng-dev libfreetype6-dev pkg-config 

# mafft
RUN apt-get install -y mafft

# fasttree
RUN apt-get install -y fasttree

# raxml
RUN apt-get install -y raxml
RUN cp /usr/bin/raxmlHPC /usr/bin/raxml

# python modules
RUN pip install biopython==1.64
RUN pip install DendroPy==3.12.0
RUN pip install seqmagick==0.5.0
RUN pip install schedule==0.3.0
RUN pip install boto==2.38.0
RUN pip install matplotlib==1.5.1
RUN pip install pandas==0.16.2
RUN pip install seaborn==0.6.0

# augur
RUN git clone https://github.com/blab/nextflu.git /nextflu
WORKDIR /nextflu/augur

# default process
CMD python -u src/Zika_process.py
