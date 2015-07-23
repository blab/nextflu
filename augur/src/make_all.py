#!/usr/bin/env python
import subprocess,glob, os, argparse
from Bio import SeqIO, AlignIO
import numpy as np

verbose = False
patterns = {('A / H3N2', ''):'H3N2',
			('A / H1N1', 'pdm09'):'H1N1pdm',
			('B / H0N0', 'Victoria'):'Vic',
			('B / H0N0', 'Yamagata'):'Yam',
			('A / H1N1', 'seasonal'):'H1N1',
			('A / H7N9', ''):'H7N9',
			('A / H5N1', ''):'H5N1',
			('A / H6N1', ''):'H6N1',
			('A / H5N6', ''):'H5N6'			
			}

outgroups = {lineage:SeqIO.read('source-data/'+lineage+'_outgroup.gb', 'genbank')
			for lineage in ['H3N2', 'H1N1pdm', 'Vic', 'Yam']}

def determine_lineage(seq):
	fields = map(lambda x:x.strip(), seq.description.split('|'))
	tmp_lineage = (fields[2], fields[4])
	if tmp_lineage in patterns:
		print fields[0],"\n\tgisaid defined lineage:",tmp_lineage,'->',patterns[tmp_lineage]
		return patterns[tmp_lineage]
	else:
		scores = []
		for olineage, oseq in outgroups.iteritems():
			if (params.aligner == "seqan"):
				from seqanpy import align_overlap
				tmp_aln = align_overlap(str(oseq.seq), str(seq.seq).replace('-','').upper(),
							score_gapopen=-10, score_gapext=-1)
				tmp_aln = np.array([np.fromstring(tmp_aln[1], dtype='|S1'), np.fromstring(tmp_aln[2], dtype='|S1')])
			if (params.aligner == "mafft"):
				SeqIO.write([oseq, seq], "temp_in.fasta", "fasta")
				os.system("mafft --auto temp_in.fasta > temp_out.fasta 2>tmp")
				tmp_aln = np.array(AlignIO.read('temp_out.fasta', 'fasta'))
			scores.append((olineage, (tmp_aln[0]==tmp_aln[1]).sum()))
		scores.sort(key = lambda x:x[1], reverse=True)
		if scores[0][1]>0.85*len(seq):
			print fields[0], tmp_lineage, len(seq), "\n\t lineage based on similarity:",scores[0][0],"\n\t",scores
			return scores[0][0]
		else:
			print fields[0], tmp_lineage, len(seq), "\n\t other: best scores:",scores[0]
			return 'other'

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
	
def push_json_to_s3(lineage, resolution, directory = '../auspice/data/', bucket = 'nextflu-dev', cloudfront = 'E1XKGZG0ZTX4YN'):
	"""Upload JSON files to S3 bucket"""
	"""Boto expects environmental variables AWS_ACCESS_KEY_ID and AWS_SECRET_ACCESS_KEY"""
	directory = directory.rstrip('/')+'/'

	import boto
	conn = boto.connect_s3()
	b = conn.get_bucket(bucket)
	k = boto.s3.key.Key(b)

	paths = []	

	print "Uploading JSONs for", lineage, resolution
	for postfix in ['tree.json', 'sequences.json', 'frequencies.json', 'meta.json']:
		json = lineage + '_' + resolution + '_' + postfix
		k.key = 'data/'+json
		k.set_contents_from_filename(directory+json)
		print json,"uploaded"
		paths.append('data/'+json)

	c = boto.connect_cloudfront()
	c.create_invalidation_request(cloudfront, paths)
	
def gather_strains(lineage, directory='data/'):
	directory = directory.rstrip('/')+'/'
	fname = lineage+'_gisaid_epiflu_sequence.fasta'
	new_strains = []
	seq_fname = directory + fname
	if not os.path.isfile(seq_fname):
		print "File with new sequences does not exist ", seq_fname
	else:
		for seq in SeqIO.parse(seq_fname, 'fasta'):
			fields = map(lambda x:x.strip(), seq.description.split('|'))
			strain = fields[0]
			new_strains.append(strain)
	return new_strains

def ammend_fasta(fname, lineage, existing_strains, threshold = 10, directory = 'data/'):
	directory = directory.rstrip('/')+'/'
	updated = False
	ex_fname = directory+lineage+'_gisaid_epiflu_sequence.fasta'

	new_seqs = []
	seq_fname = directory + fname
	if not os.path.isfile(seq_fname):
		print "File with new sequences does not exist ", seq_fname
		return updated
	else:
		for seq in SeqIO.parse(seq_fname, 'fasta'):
			fields = map(lambda x:x.strip(), seq.description.split('|'))
			strain = fields[0]
			if strain not in existing_strains:
				tmp_lineage = determine_lineage(seq)
				if tmp_lineage == lineage:
					new_seqs.append(seq)

	print "Found", len(new_seqs), 'new for lineage', lineage
	if len(new_seqs)>=threshold:
		with open(ex_fname, 'a') as outfile:
			SeqIO.write(new_seqs, outfile, 'fasta')
		updated = True

	return updated

if __name__=="__main__":
	parser = argparse.ArgumentParser(description = "ammend existing files with downloaded viruses, rerun")
	parser.add_argument('--annotate', action = "store_true", default = False, help = "annotate, but don't process")	
	parser.add_argument('--infile', type = str, default = "gisaid_epiflu_sequence.fasta")
	parser.add_argument('--bin', type = str, default = "python")
	parser.add_argument('--ATG', action = "store_true", default = False, help = "include full HA sequence starting at ATG")
	parser.add_argument('--all', action = "store_true", default = False)
	parser.add_argument('--s3', action = "store_true", default = False, help = "push/pull FASTA and JSON files to/from S3")
	parser.add_argument('--fasta_bucket', type = str, default = "nextflu-data", help = "bucket for FASTA files")		
	parser.add_argument('--json_bucket', type = str, default = "nextflu-dev", help = "bucket for JSON files")	
	parser.add_argument('--threshold', type = float, default = 10.0, help = "number of new sequences required to rerun pipeline")
	parser.add_argument('--lineages', nargs='+', type = str,  help ="lineages to include")
	parser.add_argument('--resolutions', nargs='+', type = str,  help ="resolutions to include")
	parser.add_argument('--aligner', type = str, default = "mafft", help = "aligner to use, either mafft or seqan")			
	parser.add_argument('-r', type = float, default = 1.0)
	params = parser.parse_args()

	common_args = ['--skip', 'genotype_frequencies', '-r', params.r]
	if params.ATG: common_args.append('--ATG')
	
	if params.lineages is None:
		params.lineages = ['H3N2', 'H1N1pdm', 'Vic', 'Yam']
		
	if params.resolutions is None:		
#		params.resolutions = ['1y', '3y', '6y', '12y']
		params.resolutions = ['3y', '6y', '20y']

	for lineage in params.lineages:
		if params.s3:	
			pull_fasta_from_s3(lineage, directory = 'data/', bucket = params.fasta_bucket)

	existing_strains = []
	for lineage in ['H3N2', 'H1N1pdm', 'Vic', 'Yam']:
		existing_strains.extend(gather_strains(lineage, directory = 'data/'))

	for lineage in params.lineages:
		print '\n------------------------------\n'
		print 'Parsing new sequences for lineage',lineage			
		if params.all:
			params.threshold = 0
		#run = ammend_fasta(params.infile, lineage, existing_strains, threshold = params.threshold, directory = 'data/')
		if True:
			for resolution in params.resolutions:
				print '\n------------------------------\n'
				print 'Processing lineage',lineage,'with resolution',resolution		
				process = 'src/' + lineage + '_process.py'
				if resolution == '1y':
					n_viruses = 100
					n_years = 1
				if resolution == '2y':
					n_viruses = 50
					n_years = 2
				if resolution == '3y':
					n_viruses = 50
					n_years = 3
				if resolution == '5y':
					n_viruses = 20
					n_years = 5
				if resolution == '6y':
					n_viruses = 30
					n_years = 6
				if resolution == '7y':
					n_viruses = 30
					n_years = 6
				if resolution == '20y':
					n_viruses = 10
					n_years = 20
				if resolution == '10y':
					n_viruses = 10
					n_years = 10				
				if resolution == '12y':
					n_viruses = 10
					n_years = 12
				prefix = lineage + '_'
				call = map(str, [params.bin, process, '-v', n_viruses, '-y', n_years, '--start','HI',
				           		 '--prefix', prefix, '--resolution', resolution] + common_args)
				print call
				if not params.annotate:
					subprocess.call(call)
				if params.s3:
					push_fasta_to_s3(lineage, directory = 'data/', bucket = params.fasta_bucket)
					push_json_to_s3(lineage, resolution, directory = '../auspice/data/', bucket = params.json_bucket)
