# auspice dockerfile

FROM ruby:2.2.4
MAINTAINER Trevor Bedford <trevor@bedford.io>
RUN apt-get -y update

# wget
RUN apt-get install -y wget

# git
RUN apt-get install -y git

# gems
RUN gem update
RUN gem install 'jekyll' -v '3.1.2'
RUN gem install 'therubyracer' -v '0.12.2'
RUN gem install 'less' -v '2.6.0'
RUN gem install 's3_website' -v '1.8.2'

# auspice
RUN git clone https://github.com/blab/nextflu.git /nextflu
WORKDIR /nextflu/auspice

# update
ADD http://www.timeapi.org/utc/now /tmp/bustcache
RUN git pull

# ports
EXPOSE 4000

# default process
CMD jekyll serve --host=0.0.0.0 --force_polling
