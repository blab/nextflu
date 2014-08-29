## Augur

Augur is Python package to forecast flu evolution.  It will

* import public sequence data
* build a phylogenetic tree from this data
* estimate clade fitnesses
* estimate clade frequencies
* project frequencies foreword

It is intended to be run in an always-on fashion, recomputing predictions daily and pushing predictions to a (static) website.

## Build

To run locally, you'll need Firefox, Python, pip, [muscle](http://www.drive5.com/muscle/), [RAxML](http://sco.h-its.org/exelixis/web/software/raxml/), ruby and [s3_website](https://github.com/laurilehmijoki/s3_website).  With them installed, additional dependencies can be installed with:

	pip install -r requirements.txt
	
Alternatively, you can run across platforms using [Docker](https://www.docker.com/) and the supplied [Dockerfile](Dockerfile)

	docker pull trvrb/augur
	docker run -ti -e "GISAID_USER=$GISAID_USER" -e "GISAID_PASS=$GISAID_PASS" -e "S3_KEY=$S3_KEY" -e "S3_SECRET=$S3_SECRET" -e "S3_BUCKET=$S3_BUCKET" trvrb/augur /bin/bash

Before starting Python scripts, you'll need to run:

	supervisord -c supervisord.conf
	
## Run

Python scripts are run in the following [`run.py`](augur/run.py).  This generates sequence and tree files, most notably `site/tree.json`.  This file is uploaded to Amazon S3 by running [`upload.py`](augur/upload.py).

## Environment

You will need a GISAID account and an Amazon S3 account.  Assumes environment variables:

* `GISAID_USER`: GISAID user name
* `GISAID_PASS`: GISAID password
* `S3_KEY`: Amazon S3 key
* `S3_SECRET`: Amazon S3 secret
* `S3_BUCKET`: Amazon S3 bucket

## Process

### Ingest

Using [Selenium](https://github.com/SeleniumHQ/selenium) and Python bindings to automate web crawling. [GISAID](http://platform.gisaid.org/epi3/) requires login access.  User credentials are stored in the ENV as `GISAID_USER` and `GISAID_PASS`.

### Filter

Keeps viruses with full HA1 sequences, fully specified dates, cell passage and only one sequence per strain name.

### Align

Align sequences with [muscle](http://www.drive5.com/muscle/) and strip to just the 987 bases of HA1.  This should take ~1.5 hours for ~15k sequences.

## Frequencies

Estimate clade frequencies using SMC particle filtering.
