FROM ubuntu:14.04
MAINTAINER Trevor Bedford <trevor@bedford.io>
RUN apt-get -y update

# git
RUN apt-get -y install git

# headless firefox
RUN apt-get -y install firefox xvfb x11vnc
RUN apt-get -y install -q xfonts-100dpi xfonts-75dpi xfonts-scalable xfonts-cyrillic

# python
RUN apt-get install -y python python-dev python-pip python-virtualenv

# python modules
RUN pip install selenium==2.42.1

# augur
RUN git clone https://github.com/blab/augur.git /augur
RUN cd /augur
WORKDIR /augur

# default command
CMD ["/bin/bash"]
