#!/usr/bin/env python
import argparse

verbose = False
lineages = ['H3N2', 'H1N1pdm', 'Vic', 'Yam']

def pull_from_dev(lineage, directory = '../auspice/data/', bucket = 'nextflu-dev'):
	"""Retrieve JSON files from S3 dev bucket"""
	"""Boto expects environmental variables AWS_ACCESS_KEY_ID and AWS_SECRET_ACCESS_KEY"""
	directory = directory.rstrip('/')+'/'

	import boto
	conn = boto.connect_s3()
	b = conn.get_bucket(bucket)
	k = boto.s3.key.Key(b)

	print "Retrieving JSONs for", lineage, "from bucket", bucket
	for postfix in ['_tree.json', '_sequences.json', '_frequencies.json', '_meta.json']:
		json = lineage + postfix
		k.key = 'data/'+json
		k.get_contents_to_filename(directory+json)
		print json,"retrieved"
	
def push_to_production(lineage, directory = '../auspice/data/', bucket = 'nextflu'):
	"""Upload JSON files to S3 production bucket"""
	"""Boto expects environmental variables AWS_ACCESS_KEY_ID and AWS_SECRET_ACCESS_KEY"""
	directory = directory.rstrip('/')+'/'

	import boto
	conn = boto.connect_s3()
	b = conn.get_bucket(bucket)
	k = boto.s3.key.Key(b)

	print "Uploading JSONs for", lineage, "to bucket", bucket
	for postfix in ['_tree.json', '_sequences.json', '_frequencies.json', '_meta.json']:
		json = lineage + postfix
		k.key = 'data/'+json
		k.set_contents_from_filename(directory+json)
		print json,"uploaded"

if __name__=="__main__":
	parser = argparse.ArgumentParser(description = "sync JSONs between local storage, dev server and production server")	
	parser.add_argument('--dev_bucket', type = str, default = "nextflu-dev", help = "dev bucket for JSON files")
	parser.add_argument('--live_bucket', type = str, default = "nextflu", help = "dev bucket for JSON files")
	parser.add_argument('--pull', action = "store_true", default = False, help = "pull local files from dev")
	parser.add_argument('--push', action = "store_true", default = False, help = "push local files to production")
	parser.add_argument('--sync', action = "store_true", default = False, help = "sync dev to production")	
	params = parser.parse_args()

	for lineage in lineages:
		print '\nLineage',lineage
		if params.pull or params.sync:
			pull_from_dev(lineage, directory = '../auspice/data/', bucket = params.dev_bucket)
		if params.push or params.sync:
			push_to_production(lineage, directory = '../auspice/data/', bucket = params.live_bucket)			
