## Augur

Augur is Python package to forecast flu evolution.  It will

* import public sequence data
* build a phylogenetic tree from this data
* estimate clade fitnesses
* estimate clade frequencies
* project frequencies foreword

It is intended to be run in an always-on fashion, recomputing predictions daily and pushing predictions to a (static) website.

## Build

To run locally, you'll need Firefox, Python and pip.  With them installed, additional dependencies can be installed with:

	pip install -r requirements.txt
	
Alternatively, you can run across platforms using [Docker](https://www.docker.com/) and the supplied [Dockerfile](Dockerfile)

	docker pull trvrb/augur
	docker run -ti -e "GISAID_USER=user" -e "GISAID_PASS=pass" trvrb/augur /bin/bash

where `user` and `pass` are set appropriately.

Before starting Python scripts, you'll need to run:

	supervisord -c supervisord.conf

## Run

Python scripts are run in the following order

	python augur/ingest.py
	python augur/filter.py

## Process

### Ingest

Using [Selenium](https://github.com/SeleniumHQ/selenium) and Python bindings to automate web crawling. [GISAID](http://platform.gisaid.org/epi3/) requires login access.  User credentials are stored in the ENV as `GISAID_USER` and `GISAID_PASS`.

### Filter

Keeps viruses with full HA1 sequences, fully specified dates, cell passage and only one sequence per strain name.

