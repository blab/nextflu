## Augur

Augur is the processing pipeline to track flu evolution.  It currently

* imports public sequence data
* subsamples, cleans and aligns sequences
* builds a phylogenetic tree from this data
* reports statistics about mutations and branching patterns of the tree
* infers mutation frequency trajectories through time.


## Pipeline

The entire pipeline is run with [`process.py`](src/process.py).

### Sequence download, cleaning and alignment

#### Download

Virus sequence data is manually downloaded from the [GISAID EpiFlu database](http://gisaid.org). Data from GISAID may not be disclosed outside the GISAID community. We are mindful of this and raw GISAID data has not been released publicly as part of this project. The current pipeline is designed to work specifically for HA from influenza H3N2. Save GISAID sequences as `data/gisaid_epiflu_sequence.fasta`.

#### [Filter](src/virus_filter.py)

Keeps viruses with fully specified dates, cell passage and only one sequence per strain name. Subsamples to 50 (by default) sequences per month for the last 3 (by default) years before present. Appends geographic metadata. Subsampling prefers longer sequences over shorter sequences and prefer more geographic diversity over less geographic diversity.

#### [Align](src/virus_align.py)

Aligns sequences with [mafft](http://mafft.cbrc.jp/alignment/software/).  Testing showed a much lower memory footprint than [muscle](http://www.drive5.com/muscle/).

#### [Clean](src/virus_clean.py)

Clean up alignment so that reference frame is kept intact. Remove sequences that don't conform to a rough molecular clock and remove known reassortant sequences and other outliers.

### Tree processing

#### [Infer](src/tree_infer.py)

Uses [FastTree](http://meta.microbesonline.org/fasttree/) to get a starting tree, and then refines this tree with [RAxML](http://sco.h-its.org/exelixis/web/software/raxml/). 

#### [Refine](src/tree_refine.py)

Reroot the tree based on outgroup strain, collapse nodes with zero-length branches, ladderize the tree and collect strain metadata.

### [Frequency estimation](src/bernoulli_frequency.py)

Estimate genotype and clade frequency trajectories using a Bernoulli observation model combined with a genetic drift model of process noise.

### [Streamline](src/streamline.py)

Prep and remove cruft from data files for [auspice](../auspice/) visualization.


## Develop

(This text on Docker support needs to be cleaned up).

Full dependency information can be seen in the [`Dockerfile`](Dockerfile).  To run locally, pull the docker image with

	docker pull trvrb/augur
	
And start up a bash session with

	docker run -ti -e "GISAID_USER=$GISAID_USER" -e "GISAID_PASS=$GISAID_PASS" trvrb/augur /bin/bash
	
From here, the [build pipeline](augur/run.py) can be run with

	python augur/run.py

You can run across platforms using [Docker](https://www.docker.com/).  An image is up on the Docker hub repository as [trvrb/augur](https://registry.hub.docker.com/u/trvrb/augur/).  With this public image, you can immediately run augur with

	docker pull trvrb/augur
	docker run -ti -e "GISAID_USER=$GISAID_USER" -e "GISAID_PASS=$GISAID_PASS" -e "S3_KEY=$S3_KEY" -e "S3_SECRET=$S3_SECRET" -e "S3_BUCKET=$S3_BUCKET" --privileged trvrb/augur
	
This starts up [Supervisor](http://supervisord.org/) to keep augur and helper programs running.  This uses [`supervisord.conf`](supervisord.conf) as a control file.

To run augur, you will need a GISAID account (to pull sequences) and an Amazon S3 account (to push results).  Account information is stored in environment variables:

* `GISAID_USER`: GISAID user name
* `GISAID_PASS`: GISAID password
* `S3_KEY`: Amazon S3 key
* `S3_SECRET`: Amazon S3 secret
* `S3_BUCKET`: Amazon S3 bucket
