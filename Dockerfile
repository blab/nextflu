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

# avx raxml
RUN mkdir -p /raxml
RUN curl -o /raxml/v8.1.1 https://codeload.github.com/stamatak/standard-RAxML/tar.gz/v8.1.1
RUN tar xvzf /raxml/v8.1.1 -C /raxml/
WORKDIR /raxml/standard-RAxML-8.1.1/
RUN make -f Makefile.AVX.PTHREADS.gcc
RUN mv raxmlHPC-PTHREADS-AVX /usr/bin/raxml

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
RUN git clone https://github.com/blab/augur.git /augur # fdafas
RUN mkdir -p /augur/log
WORKDIR /augur

# default command
CMD supervisord -c supervisord.conf

