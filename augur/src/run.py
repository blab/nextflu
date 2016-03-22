#!/usr/bin/env python
import os, re, time, datetime, csv, sys, argparse, subprocess
import rethinkdb as r
from vdb_download import vdb_download
from Bio import SeqIO, AlignIO

def local_count(params):
	handle = open(params.path + params.fstem + "." + params.ftype, "rU")
	records = list(SeqIO.parse(handle, "fasta"))
	return len(records)

def build(params):
	'''
	run build for single virus / resolution combination
	'''
	vdb.download_all_documents()
	vdb.output()

	process = 'src/' + params.virus + '_process.py'
	call = map(str, [params.bin, process, '--resolution', params.resolution])
	print call
	subprocess.call(call)

def tick():
	print "tick"

if __name__=="__main__":
	parser = argparse.ArgumentParser(description = "download and process")
	parser.add_argument('--bin', type = str, default = "python")	
	parser.add_argument('-v', '--virus', default='Zika', help="virus table to interact with")
	parser.add_argument('-r', '--resolution', default='', help="build params / annotation")	
	parser.add_argument('--host', default=None, help="rethink host url")
	parser.add_argument('--auth_key', default=None, help="auth_key for rethink database")
	parser.add_argument('--build', action="store_true", default=False, help ="single rebuild")
	parser.add_argument('--watch', action="store_true", default=False, help ="watch database and rebuild when updated")	
	params = parser.parse_args()

	if params.virus is None:
		params.virus = 'Zika'

	params.database = 'vdb'
	params.path = 'data/'
	params.ftype = 'fasta'	
	params.fstem = params.virus
	
	vdb = vdb_download(**params.__dict__)

	# rebuild
	if (params.build):
		build(params)

	# watch
	if (params.watch):
		while True:
			rcount = vdb.count_documents()
			lcount = local_count(params)
			print rcount, "remote documents and", lcount, "local documents"
			if (rcount > lcount):
				build(params)
			time.sleep(60)
