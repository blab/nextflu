## Augur

Augur is Python package to track (and eventually forecast) flu evolution.  It currently

* imports public sequence data
* subsamples, cleans and aligns sequences
* builds a phylogenetic tree from this data
* Reports statistics about mutations and branching patterns of the tree

Results are pushed to Amazon S3.  The current JSON-formatted flu tree is available as [`tree.json`](https://s3.amazonaws.com/augur-data/auspice/tree.json).  This tree is visualized at [blab.github.io/auspice/](http://blab.github.io/auspice/).

## Run

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

## Develop

Full dependency information can be seen in the [`Dockerfile`](Dockerfile).  To run locally, pull the docker image with

	docker pull trvrb/augur
	
And start up a bash session with

	docker run -ti -e "GISAID_USER=$GISAID_USER" -e "GISAID_PASS=$GISAID_PASS" trvrb/augur /bin/bash
	
From here, the [build pipeline](augur/run.py) can be run with

	python augur/run.py
	
## Pipeline notes

### Sequence download, cleaning and alignment

#### Download

Virus sequence data is manually downloaded from the [GISAID EpiFlu database](http://gisaid.org). Data from GISAID may not be disclosed outside the GISAID community. We mindful of this and no sequence data is has been released publicly as part of this project. We are providing processed phylogenies only. Save GISAID sequences as `data/gisaid_epiflu_sequence.fasta`.

#### [Filter](augur/virus_filter.py)

Keeps viruses with full HA1 sequences, fully specified dates, cell passage and only one sequence per strain name.  Subsamples to 100 sequences per month for the last 3 years before present.

#### [Align](augur/virus_align.py)

Align sequences with [mafft](http://mafft.cbrc.jp/alignment/software/).  Testing showed a much lower memory footprint than [muscle](http://www.drive5.com/muscle/).

#### [Clean](augur/virus_clean.py)

Keep only sequences that have the full 1701 bases of HA in the alignment.

### Tree processing

#### [Infer](augur/tree_infer.py)

Using [FastTree](http://meta.microbesonline.org/fasttree/) to get a starting tree.  FastTree will build a tree for ~5000 sequences in a few minutes.  Then using [RAxML](http://sco.h-its.org/exelixis/web/software/raxml/) to refine this initial tree.  A full RAxML run on a tree with ~5000 sequences could take days or weeks, so instead RAxML is run for a fixed 1 hour and the best tree found during this search is kept.  This will always improve on FastTree. RAxML is also used to reconstruct ancestral sequences.

#### [Refine](augur/tree_refine.py)

Reroot the tree based on outgroup strain, collapse nodes with zero-length branches, ladderize the tree and calculate amino acid distances.
