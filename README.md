## Augur

Augur is Python package to track (and eventually forecast) flu evolution.  It currently

* imports public sequence data
* subsamples, cleans and aligns sequences
* builds a phylogenetic tree from this data

The program is live on Amazon EC2 with results pushed to Amazon S3.  The latest JSON-formatted flu tree is available as [`tree_streamline.json`](https://s3.amazonaws.com/augur-data/data/tree_streamline.json).  This tree is visualized at [blab.github.io/auspice/](http://blab.github.io/auspice/).

## Run

You can run across platforms using [Docker](https://www.docker.com/).  An image is up on the Docker hub repository as [trvrb/augur](https://registry.hub.docker.com/u/trvrb/augur/).  With this public image, you can immediately run augur with

	docker pull trvrb/augur
	docker run -ti -e "GISAID_USER=$GISAID_USER" -e "GISAID_PASS=$GISAID_PASS" -e "S3_KEY=$S3_KEY" -e "S3_SECRET=$S3_SECRET" -e "S3_BUCKET=$S3_BUCKET" --privileged trvrb/augur
	
This starts up [Supervisor](http://supervisord.org/) to keep augur running and other helper programs, which can be seen in the [`supervisord.conf`](supervisord.conf) control file.

To run augur, you will need a GISAID account (to pull sequences) and an Amazon S3 account (to push results).  Account information is stored as environment variables:

* `GISAID_USER`: GISAID user name
* `GISAID_PASS`: GISAID password
* `S3_KEY`: Amazon S3 key
* `S3_SECRET`: Amazon S3 secret
* `S3_BUCKET`: Amazon S3 bucket

## Develop

Full dependency information can be seen in the [`Dockerfile`](Dockerfile).  To run locally, pull the docker image with

	docker pull trvrb/augur
	
And start up a bash session with

	docker run -ti -e "GISAID_USER=$GISAID_USER" -e "GISAID_PASS=$GISAID_PASS" trvrb/augur /bin/bash
	
From here, the [build pipeline](augur/run.py) can be run with

	python augur/run.py
	
## Process notes

### Virus ingest, alignment and filtering

#### [Ingest](augur/virus_ingest.py)

Using [Selenium](https://github.com/SeleniumHQ/selenium) and Python bindings to automate downloads from [GISAID](http://platform.gisaid.org/epi3/).  GISAID requires login access.  User credentials are stored in the ENV as `GISAID_USER` and `GISAID_PASS`.

#### [Filter](augur/virus_filter.py)

Keeps viruses with full HA1 sequences, fully specified dates, cell passage and only one sequence per strain name.  Subsamples to 100 sequences per month for the last 3 before present.

#### [Align](augur/virus_align.py)

Align sequences with [mafft](http://mafft.cbrc.jp/alignment/software/).  Testing showed a much lower memory footprint than [muscle](http://www.drive5.com/muscle/).

#### [Clean](augur/virus_clean.py)

Keep only sequences that have the full 987 bases of HA1 in the alignment.

### Tree processing

#### [Infer](augur/tree_infer.py)

Using [FastTree](http://meta.microbesonline.org/fasttree/) to get a starting tree.  FastTree will build a tree for ~5000 sequences in a few minutes.  Then using [RAxML](http://sco.h-its.org/exelixis/web/software/raxml/) to refine this initial tree.  A full RAxML run on a tree with ~5000 sequences could take days or weeks, so instead RAxML is run for a fixed 1 hour and the best tree found during this search is kept.  This will always improve on FastTree.

#### [Clean](augur/tree_clean.py)

Reroot the tree based on the Beijing/32/1992 outgroup strain, collapse nodes with zero-length branches and ladderize the tree.

