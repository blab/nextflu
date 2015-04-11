#!/usr/bin/env python
import argparse

verbose = False
lineages = ['H3N2', 'H1N1pdm', 'Vic', 'Yam']

def pull(lineage, directory = '../auspice/data/', bucket = 'nextflu-dev'):
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
	
def push(lineage, directory = '../auspice/data/', bucket = 'nextflu-dev', cloudfront = 'E1XKGZG0ZTX4YN'):
	"""Upload JSON files to S3 production bucket"""
	"""Boto expects environmental variables AWS_ACCESS_KEY_ID and AWS_SECRET_ACCESS_KEY"""
	directory = directory.rstrip('/')+'/'

	import boto
	conn = boto.connect_s3()
	b = conn.get_bucket(bucket)
	k = boto.s3.key.Key(b)

	paths = []

	print "Uploading JSONs for", lineage, "to bucket", bucket
	for postfix in ['_tree.json', '_sequences.json', '_frequencies.json', '_meta.json']:
		json = lineage + postfix
		k.key = 'data/'+json
		k.set_contents_from_filename(directory+json)
		print json,"uploaded"
		paths.append('data/'+json)

	c = boto.connect_cloudfront()
	c.create_invalidation_request(cloudfront, paths)

if __name__=="__main__":
	parser = argparse.ArgumentParser(description = "sync JSONs between local storage, dev server and production server")	
	parser.add_argument('--dev_bucket', type = str, default = "nextflu-dev", help = "dev bucket for JSON files")
	parser.add_argument('--pro_bucket', type = str, default = "nextflu", help = "production bucket for JSON files")
	parser.add_argument('--dev_cloudfront', type = str, default = "E1XKGZG0ZTX4YN", help = "dev cloudfront distribution")
	parser.add_argument('--pro_cloudfront', type = str, default = "E1ZD3XPSGXABBI", help = "production cloudfront distribution")	
	parser.add_argument('--pull_dev', action = "store_true", default = False, help = "pull local files from dev")
	parser.add_argument('--push_dev', action = "store_true", default = False, help = "push local files to dev")
	parser.add_argument('--pull_pro', action = "store_true", default = False, help = "pull local files from pro")
	parser.add_argument('--push_pro', action = "store_true", default = False, help = "push local files to pro")	
	parser.add_argument('--sync', action = "store_true", default = False, help = "sync dev to production")	
	params = parser.parse_args()

	for lineage in lineages:
		print '\nLineage',lineage
		if params.pull_dev or params.sync:
			pull(lineage, directory = '../auspice/data/', bucket = params.dev_bucket)
		if params.push_dev:
			push(lineage, directory = '../auspice/data/', bucket = params.dev_bucket, cloudfront = params.dev_cloudfront)
		if params.pull_pro:
			pull(lineage, directory = '../auspice/data/', bucket = params.pro_bucket)
		if params.push_pro or params.sync:
			push(lineage, directory = '../auspice/data/', bucket = params.pro_bucket, cloudfront = params.pro_cloudfront)						
