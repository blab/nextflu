FROM ubuntu:14.04
MAINTAINER Trevor Bedford <trevor@bedford.io>
RUN apt-get -y update

# git
RUN apt-get -y install git

# rethinkdb
#RUN echo "deb http://download.rethinkdb.com/apt `lsb_release -cs` main" > /etc/apt/sources.list.d/rethinkdb.list
#RUN wget -O- http://download.rethinkdb.com/apt/pubkey.gpg | apt-key add -
#RUN apt-get install -y rethinkdb

# headless firefox
RUN apt-get -y install firefox xvfb x11vnc
RUN apt-get -y install -q xfonts-100dpi xfonts-75dpi xfonts-scalable xfonts-cyrillic
ENV DISPLAY :99

# python
RUN apt-get install -y python python-dev python-pip python-virtualenv

# python modules
RUN pip install selenium==2.42.1

# augur
RUN git clone https://github.com/blab/augur.git /augur
RUN cd /augur
WORKDIR /augur

# supervisor
RUN apt-get -y install supervisor
RUN mkdir -p /var/log/supervisor

# default command
CMD ["supervisord -c supervisord.conf"]

# Expose ports
#   - 8080: web UI
#   - 28015: process
#   - 29015: cluster
#EXPOSE 8080
