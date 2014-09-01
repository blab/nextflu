## Augur

Augur is Python package to forecast flu evolution.  It will

* import public sequence data
* build a phylogenetic tree from this data
* estimate clade fitnesses
* estimate clade frequencies
* project frequencies foreword

It is intended to be run in an always-on fashion, recomputing predictions daily and pushing predictions to a (static) website.

## Build and run

You can run across platforms using [Docker](https://www.docker.com/) and the supplied [Dockerfile](Dockerfile)

	docker pull trvrb/augur
	docker run -ti -e "GISAID_USER=$GISAID_USER" -e "GISAID_PASS=$GISAID_PASS" -e "S3_KEY=$S3_KEY" -e "S3_SECRET=$S3_SECRET" -e "S3_BUCKET=$S3_BUCKET" --privileged trvrb/augur
	
This starts up [Supervisor](http://supervisord.org/) to keep augur running and other helper programs, which can be seen in the [`supervisord.conf`](supervisord.conf) control file.

To run locally, you'll need Firefox, Python, pip, [mafft](http://mafft.cbrc.jp/alignment/software/), [FastTree](http://meta.microbesonline.org/fasttree/) and some other things as well.  A complete listing can be seen in the [Dockerfile](Dockerfile).
	
The build pipeline is initiated with [`run.py`](augur/run.py), generating sequence and tree files.
	
## Environment

You will need a GISAID account and an Amazon S3 account.  Assumes environment variables:

* `GISAID_USER`: GISAID user name
* `GISAID_PASS`: GISAID password
* `S3_KEY`: Amazon S3 key
* `S3_SECRET`: Amazon S3 secret
* `S3_BUCKET`: Amazon S3 bucket

## Run on Amazon EC2

Install EC2 command line tools:

	brew install ec2-api-tools
	
Set up environment variables, see: 

	brew info ec2-api-tools

Set up key pair:

	ec2-add-keypair ec2-keypair > ~/.ec2-keypair.pem
	chmod 600 ~/.ec2-keypair.pem

Open up ports:

	ec2-authorize default -p 22
	ec2-authorize default -p 80
	
Start instance:	
	
	ec2-run-instances -k ec2-keypair -f ec2-startup.sh -t t2.micro -z us-east-1a ami-864d84ee	# ubuntu
	ec2-run-instances -k ec2-keypair -f ec2-startup.sh -t t2.micro -z us-east-1a ami-9846e2f0	# ubuntu + augur
	
SSH in:

	ssh -i ~/.ec2-keypair.pem ubuntu@ec2-xxx.amazonaws.com
	
Start augur:

	sudo docker run -d -e "GISAID_USER=$GISAID_USER" -e "GISAID_PASS=$GISAID_PASS" -e "S3_KEY=$S3_KEY" -e "S3_SECRET=$S3_SECRET" -e "S3_BUCKET=$S3_BUCKET" --privileged trvrb/augur

Terminate instance:

	ec2-terminate-instances i-xxx

## Process

### Ingest

Using [Selenium](https://github.com/SeleniumHQ/selenium) and Python bindings to automate web crawling. [GISAID](http://platform.gisaid.org/epi3/) requires login access.  User credentials are stored in the ENV as `GISAID_USER` and `GISAID_PASS`.

### Filter

Keeps viruses with full HA1 sequences, fully specified dates, cell passage and only one sequence per strain name.

### Align

Align sequences with [mafft](http://mafft.cbrc.jp/alignment/software/).  Testing showed a much lower memory footprint than [muscle](http://www.drive5.com/muscle/).

## Tree building

Using [FastTree](http://meta.microbesonline.org/fasttree/) after really attempting [RAxML](http://sco.h-its.org/exelixis/web/software/raxml/).  Time spans just don't work.  Running FastTree with double precision to distinguish single substitution branches from zero substitution branches.
