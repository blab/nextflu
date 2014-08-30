FROM ubuntu:14.04
MAINTAINER Trevor Bedford <trevor@bedford.io>
RUN apt-get -y update

# wget
RUN apt-get install -y wget

# git
RUN apt-get install -y git

# headless firefox
RUN apt-get install -y firefox xvfb x11vnc
RUN apt-get install -y -q xfonts-100dpi xfonts-75dpi xfonts-scalable xfonts-cyrillic
ENV DISPLAY :99

# python
RUN apt-get install -y python python-dev python-pip python-virtualenv
RUN apt-get install -y python-numpy python-scipy

# muscle
RUN apt-get install -y muscle

# raxml
RUN apt-get install -y raxml
RUN mv /usr/bin/raxmlHPC /usr/bin/raxml

# python modules
RUN pip install selenium==2.42.1
RUN pip install biopython==1.63
RUN pip install seqmagick==0.5.0
RUN pip install ete2==2.2.1072
RUN pip install DendroPy==3.12.0
RUN pip install schedule==0.3.0

# s3 website
RUN apt-get install -y ntp
RUN apt-get install -y ruby
RUN gem install s3_website

# java (required for s3 website)
RUN apt-get update -y
RUN apt-get install -y openjdk-7-jre
RUN rm -rf /var/lib/apt/lists/*

# augur
RUN git clone https://github.com/blab/augur.git /augur
WORKDIR /augur

EXPOSE 80

# default command
CMD ["/augur/docker_run.sh"]

