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

def ammend_files(fname, lineages= ['H3N2', 'H1N1pdm', 'Vic', 'Yam'], threshold = 10, directory = 'data/'):
	directory = directory.rstrip('/')+'/'
	updated = []

	for lineage in lineages:
		print '\nLineage',lineage
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
			continue
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
			updated.append(lineage)
	return updated


if __name__=="__main__":
	parser = argparse.ArgumentParser(description = "ammend existing files with downloaded viruses, rerun")
	parser.add_argument('--infile', type = str, default = "gisaid_epiflu_sequence.fasta")
	parser.add_argument('--bin', type = str, default = "python")
	parser.add_argument('--HA1', action = "store_true", default = False)
	parser.add_argument('--all', action = "store_true", default = False)
	parser.add_argument('-r', type = float, default = 1.0)
	params = parser.parse_args()

	if params.all:	
		run_pipeline = ammend_files(params.infile, lineages = ['H3N2', 'H1N1pdm', 'Vic', 'Yam'], threshold = 1, directory = 'data/')
		run_pipeline = ['H3N2', 'H1N1pdm', 'Vic', 'Yam']
	else:
		run_pipeline = ammend_files(params.infile, lineages = ['H3N2', 'H1N1pdm', 'Vic', 'Yam'], threshold = 10, directory = 'data/')

	common_args = ['--skip', 'genotype_frequencies','-r', params.r]
	if params.HA1: common_args.append('--HA1')

	if 'H3N2' in run_pipeline:
		call = map(str, [params.bin, 'src/H3N2_process.py', '-v', 50, '-y', 3,  '--prefix', 'data/H3N2_'])
		print call
		subprocess.call(call)
	if 'H1N1' in run_pipeline:
		call = map(str, [params.bin, 'src/H1N1historical_process.py', '-v', 20, '--interval', 1990, 2010, '--prefix', 'data/H1N1_'])
		print call
		subprocess.call(call)
	if 'H1N1pdm' in run_pipeline:
		call = map(str, [params.bin, 'src/H1N1pdm_process.py', '-v', 30, '-y', 6, '--prefix', 'data/H1N1pdm_'])
		print call
		subprocess.call(call)
	if 'Vic' in run_pipeline:
		call = map(str, [params.bin, 'src/Vic_process.py', '-v', 30, '-y', 6, '--prefix', 'data/Vic_'])
		print call
		subprocess.call(call)
	if 'Yam' in run_pipeline:
		call = map(str, [params.bin, 'src/Yam_process.py', '-v', 30, '-y', 6, '--prefix', 'data/Yam_'])
		print call
		subprocess.call(call)
