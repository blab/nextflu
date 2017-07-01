#!/usr/bin/env python
import os, re, time, datetime, csv, sys, argparse, subprocess

def flu_build(lineage, resolution, passage):
	'''
	run build for single lineage / resolution / passage combination
	'''
	print '\n------------------------------\n'
	print 'Processing lineage', lineage, 'with resolution', resolution, 'and passage', passage
	if resolution == '2y':
		n_viruses = 100
		n_years = 2
	if resolution == '3y':
		n_viruses = 60
		n_years = 3
	if resolution == '6y':
		n_viruses = 24
		n_years = 6
	if resolution == '12y':
		n_viruses = 12
		n_years = 12
	process = 'src/' + lineage + '_process.py'
	HI_fname = 'data/' + lineage + '_cdc_hi_' + passage + '_titers.tsv'
	force_include = 'data/' + lineage + '_cdc_hi_' + passage + '_strains.tsv'
	call = map(str, [params.bin, process, '-v', n_viruses, '-y', n_years,
		'--resolution', resolution + '_' + passage,
		'--HI_fname', HI_fname, '--force_include', force_include,
		'--skip', 'genotype_frequencies fitness HIvalidate'])
	print call
	subprocess.call(call)

if __name__=="__main__":
	parser = argparse.ArgumentParser(description = "download and process")
	parser.add_argument('--bin', type = str, default = "python")
	parser.add_argument('--zika_lineages', nargs='+', type = str,  help ="zika lineages to include")
	parser.add_argument('--zika_resolutions', nargs='+', type = str,  help ="zika resolutions to include")
	parser.add_argument('--flu_lineages', nargs='+', type = str,  help ="flu lineages to include")
	parser.add_argument('--flu_resolutions', nargs='+', type = str,  help ="flu resolutions to include")
	parser.add_argument('--flu_titers_passages', nargs='+', type = str,  help ="flu titers passages to include")
	params = parser.parse_args()

	if params.flu_lineages is None:
		params.flu_lineages = ['H3N2', 'H1N1pdm', 'Vic', 'Yam']

	if params.flu_resolutions is None:
		params.flu_resolutions = ['2y', '3y', '6y', '12y']

	if params.flu_titers_passages is None:
		params.flu_titers_passages = ['cell', 'egg']

	for lineage in params.flu_lineages:
		for resolution in params.flu_resolutions:
			for passage in params.flu_titers_passages:
				flu_build(lineage, resolution, passage)
