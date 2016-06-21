#!/usr/bin/env python
import os, re, time, datetime, csv, sys, argparse, subprocess
import rethinkdb as r
from vdb_download import vdb_download
from Bio import SeqIO, AlignIO

def local_count(params):
	count = 0
	fname = params.path + params.fstem + "." + params.ftype
	if os.path.isfile(fname):
		handle = open(fname, "rU")
		records = list(SeqIO.parse(handle, "fasta"))
		count = len(records)
	return count

def build(params):
	'''
	run build for single virus / resolution combination
	'''
	vdb.download()
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
	parser.add_argument('--zika_lineages', nargs='+', type = str,  help ="zika lineages to include")
	parser.add_argument('--zika_resolutions', nargs='+', type = str,  help ="zika resolutions to include")
	parser.add_argument('--flu_lineages', nargs='+', type = str,  help ="flu lineages to include")
	parser.add_argument('--flu_resolutions', nargs='+', type = str,  help ="flu resolutions to include")	
	parser.add_argument('--host', default=None, help="rethink host url")
	parser.add_argument('--auth_key', default=None, help="auth_key for rethink database")
	parser.add_argument('--build', action="store_true", default=False, help ="single rebuild")
	parser.add_argument('--watch', action="store_true", default=False, help ="watch database and rebuild when updated")	
	params = parser.parse_args()

	if params.zika_lineages is None:
		params.zika_lineages = ['Zika']

	if params.zika_resolutions is None:
		params.zika_resolutions = ['']	

	if params.flu_lineages is None:
		params.flu_lineages = ['H3N2', 'H1N1pdm', 'Vic', 'Yam']

	if params.flu_resolutions is None:
		params.flu_resolutions = ['3y', '6y', '12y']

	params.database = 'vdb'
	params.path = 'data/'
	params.ftype = 'fasta'	
	params.fstem = params.virus
	params.public_only = True
	params.countries = None
	
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
