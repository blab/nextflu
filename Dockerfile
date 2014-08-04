FROM ubuntu:14.04
MAINTAINER Trevor Bedford <trevor@bedford.io>
RUN apt-get -y update

# git
RUN apt-get -y install git

# headless firefox
RUN apt-get -y install firefox
RUN apt-get -y install xvfb
RUN apt-get -y install x11vnc
RUN apt-get -y install -q xfonts-100dpi xfonts-75dpi xfonts-scalable xfonts-cyrillic

# augur
RUN git clone https://github.com/blab/augur.git /augur
RUN cd /augur
WORKDIR /augur
