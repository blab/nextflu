FROM phusion/baseimage:0.9.13
MAINTAINER Trevor Bedford <trevor@bedford.io>
RUN apt-get -y update

# wget
RUN apt-get install -y wget

# git
RUN apt-get install -y git

# python
RUN apt-get install -y python python-dev python-pip python-virtualenv
RUN apt-get install -y python-numpy python-scipy

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

# libgit2 and python bindings
RUN apt-get install -y libgit2
RUN pip install pygit2==0.22.0

# s3cmd
RUN apt-get install -y s3cmd

# augur
RUN git clone https://github.com/blab/nextflu.git /nextflu
WORKDIR /nextflu

# default command
CMD python -u augur/src/process.py
