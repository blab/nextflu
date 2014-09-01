FROM phusion/baseimage:0.9.13
MAINTAINER Trevor Bedford <trevor@bedford.io>
RUN apt-get -y update

# wget
RUN apt-get install -y wget

# git
RUN apt-get install -y git

# headless firefox
RUN apt-get install -y firefox
RUN apt-get install -y xvfb
RUN apt-get install -y x11vnc
RUN apt-get install -y -q xfonts-100dpi xfonts-75dpi xfonts-scalable xfonts-cyrillic
ENV DISPLAY :99

# python
RUN apt-get install -y python python-dev python-pip python-virtualenv
RUN apt-get install -y python-numpy python-scipy

# mafft
RUN apt-get install -y mafft

# fasttree
RUN apt-get install -y fasttree

# raxml
RUN apt-get install -y raxml
RUN mv /usr/bin/raxmlHPC /usr/bin/raxml

# python modules
RUN pip install selenium==2.42.1
RUN pip install biopython==1.63
RUN pip install DendroPy==3.12.0
RUN pip install seqmagick==0.5.0
RUN pip install schedule==0.3.0

# s3cmd
RUN apt-get install -y s3cmd

# supervisor
RUN pip install supervisor==3.1.1

# dstat
RUN apt-get install -y dstat

# augur
RUN git clone https://github.com/blab/augur.git /augur
RUN mkdir -p /augur/log
WORKDIR /augur

# default command
CMD supervisord -c supervisord.conf

