#!/usr/bin/env python
import subprocess,glob, os, argparse
from Bio import SeqIO

verbose = False
patterns = {('A / H3N2', ''):'H3N2',
			('A / H1N1', 'pdm09'):'H1N1pdm',
			('B / H0N0', 'Victoria'):'Vic',
			('B / H0N0', 'Yamagata'):'Yam',
			('A / H1N1', 'seasonal'):'H1N1',
			}
lineages = ['H3N2', 'H1N1pdm', 'Vic', 'Yam']

def pull_fasta_from_s3(lineage, directory = 'data/', bucket = 'nextflu-data'):
	"""Retrieve FASTA files from S3 bucket"""
	"""Boto expects environmental variables AWS_ACCESS_KEY_ID and AWS_SECRET_ACCESS_KEY"""
	directory = directory.rstrip('/')+'/'

	import boto
	conn = boto.connect_s3()
	b = conn.get_bucket(bucket)
	k = boto.s3.key.Key(b)

	print "Retrieving FASTA for",lineage
	fasta = lineage+'_gisaid_epiflu_sequence.fasta'
	k.key = fasta
	k.get_contents_to_filename(directory+fasta)
	print fasta,"retrieved"
		
def push_fasta_to_s3(lineage, directory = 'data/', bucket = 'nextflu-data'):
	"""Upload FASTA files to S3 bucket"""
	"""Boto expects environmental variables AWS_ACCESS_KEY_ID and AWS_SECRET_ACCESS_KEY"""
	directory = directory.rstrip('/')+'/'

	import boto
	conn = boto.connect_s3()
	b = conn.get_bucket(bucket)
	k = boto.s3.key.Key(b)

	print "Uploading FASTA for",lineage
	fasta = lineage+'_gisaid_epiflu_sequence.fasta'
	k.key = fasta
	k.set_contents_from_filename(directory+fasta)
	print fasta,"uploaded"

def ammend_fasta(fname, lineage, threshold = 10, directory = 'data/'):
	directory = directory.rstrip('/')+'/'
	updated = False

	ex_fname = directory+lineage+'_gisaid_epiflu_sequence.fasta'
	existing = set()
	if not os.path.isfile(ex_fname):
		print "No existing file found for",lineage, ex_fname
	else:
		for seq in SeqIO.parse(ex_fname, 'fasta'):
			acc = int(seq.description.split('|')[-1].strip())
			existing.add(acc)
	print "Found", len(existing), 'existing for lineage', lineage 

	new_seqs = []
	seq_fname = directory + fname
	if not os.path.isfile(seq_fname):
		print "File with new sequences does not exist ", seq_fname
		return updated
	else:
		for seq in SeqIO.parse(seq_fname, 'fasta'):
			fields = map(lambda x:x.strip(), seq.description.split('|'))
			tmp_lineage = (fields[2], fields[4])
			if tmp_lineage in patterns:
				if patterns[tmp_lineage]==lineage:
					acc = int(fields[-1])
					if acc not in existing:
						new_seqs.append(seq)
			else:
				if verbose:
					print tmp_lineage,"not found"

	print "Found", len(new_seqs), 'new for lineage', lineage
	if len(new_seqs)>threshold:
		with open(ex_fname, 'a') as outfile:
			SeqIO.write(new_seqs, outfile, 'fasta')
		updated = True

	return updated


if __name__=="__main__":
	parser = argparse.ArgumentParser(description = "ammend existing files with downloaded viruses, rerun")
	parser.add_argument('--infile', type = str, default = "gisaid_epiflu_sequence.fasta")
	parser.add_argument('--bin', type = str, default = "python")
	parser.add_argument('--ATG', action = "store_true", default = False, help = "include full HA sequence starting at ATG")
	parser.add_argument('--all', action = "store_true", default = False)
	parser.add_argument('--s3', action = "store_true", default = False, help = "pull FASTA files from S3")
	parser.add_argument('--threshold', type = float, default = 10.0, help = "number of new sequences required to rerun pipeline")	
	parser.add_argument('-r', type = float, default = 1.0)
	params = parser.parse_args()

	common_args = ['--skip', 'genotype_frequencies','-r', params.r]
	if params.ATG: common_args.append('--ATG')

	for lineage in lineages:
		print '\nLineage',lineage	
		if params.s3:
			pull_fasta_from_s3(lineage, directory = 'data/', bucket = 'nextflu-data')
		if params.all:
			params.threshold = 0
		run = ammend_fasta(params.infile, lineage, threshold = params.threshold, directory = 'data/')
		if run:
			if lineage == 'H3N2':
				call = map(str, [params.bin, 'src/H3N2_process.py', '-v', 50, '-y', 3,  '--prefix', 'data/H3N2_'] + common_args)
			if lineage == 'H1N1':
				call = map(str, [params.bin, 'src/H1N1historical_process.py', '-v', 20, '--interval', 1990, 2010, '--prefix', 'data/H1N1_']+ common_args)
			if lineage == 'H1N1pdm':
				call = map(str, [params.bin, 'src/H1N1pdm_process.py', '-v', 30, '-y', 6, '--prefix', 'data/H1N1pdm_']+ common_args)
			if lineage == 'Vic':
				call = map(str, [params.bin, 'src/Vic_process.py', '-v', 30, '-y', 6, '--prefix', 'data/Vic_'] + common_args)
			if lineage == 'Yam':
				call = map(str, [params.bin, 'src/Yam_process.py', '-v', 30, '-y', 6, '--prefix', 'data/Yam_'] + common_args)
			print call
			subprocess.call(call)
			if params.s3:
				push_fasta_to_s3(lineage, directory = 'data/', bucket = 'nextflu-data')			
