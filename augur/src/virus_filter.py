# filter viruses after ingest, criteria based on metadata
#  - viruses equal to or longer than 1701 bases
#  - viruses with exact dates
#  - viruses that are not egg-passaged
#  - a single sequence per virus strain, taken as first sequence in list
# outputs to virus_filter.json

import os, re, time, datetime, csv, sys
from io_util import *
sys.path.append('../source-data')
from H3N2_outgroup_and_vaccine import outgroup, vaccine_strains

def parse_gisaid(fasta):
	"""Parse FASTA file from GISAID with default header formating"""
	viruses = []
	try:
		handle = open(fasta, 'r')
	except IOError:
		print fasta + " not found"
	else:
		for record in SeqIO.parse(handle, "fasta"):
			words = record.description.replace(">","").replace(" ","").split('|')
			strain = words[0]
			accession = words[1]
			passage = words[3]
			date = words[5]
			seq = str(record.seq).upper()
			v = {
				"strain": strain,
				"date": date,
				"accession": accession,
				"db": "GISAID",
				"seq": seq
			}
			if passage != "":
				v['passage'] = passage
			viruses.append(v)
		handle.close()

	return viruses

def sort_length(viruses):
	"""Sort by length, but randomize viruses of a given length"""
	from random import shuffle
	shuffle(viruses)
	return sorted(viruses, key = lambda v: len(v['seq']), reverse = True)

def fix_strain_names(viruses):
	for v in viruses:
		v['strain'] = v['strain'].replace('\'','').replace('(','').replace(')','').replace('H3N2','').replace('Human','').replace('human','').replace('//','/')

def filter_strain_names(viruses):
	filtered_viruses = filter(lambda v: re.match(r'^A/', v['strain']) != None, viruses)
	return filtered_viruses

def filter_length(viruses):
	return filter(lambda v: len(v['seq']) >= 987, viruses)

def filter_date(viruses):
	return filter(lambda v: re.match(r'\d\d\d\d-\d\d-\d\d', v['date']) != None, viruses)

def filter_passage(viruses):
	round_one = filter(lambda v: re.match(r'^E\d+', v.get('passage',''), re.I) == None, viruses)
	return filter(lambda v: re.match(r'^Egg', v.get('passage',''), re.I) == None, round_one)

def filter_unique(viruses):
	"""Keep only the first isolate of a strain"""
	filtered_viruses = {}
	for v in viruses:
		label = v['strain'].lower() 
		if not label in filtered_viruses:
			filtered_viruses[label] = v
	return filtered_viruses.values()

def append_country_and_region(viruses):
	"""Label viruses with geographic location based on strain name"""
	"""Location is to the level of country of administrative division when available"""
	reader = csv.DictReader(open("source-data/geo_synonyms.tsv"), delimiter='\t')		# list of dicts
	label_to_country = {}
	for line in reader:
		label_to_country[line['label'].lower()] = line['country']
	for v in viruses:
		label = re.match(r'^A/([^/]+)/', v['strain']).group(1).lower()					# check first for whole geo match
		v['country'] = 'Unknown'
		if label in label_to_country:
			v['country'] = label_to_country[label]
		else:
			label = re.match(r'^A/([^\-^\/]+)[\-\/]', v['strain']).group(1).lower()		# check for partial geo match
			if label in label_to_country:
				v['country'] = label_to_country[label]

	reader = csv.DictReader(open("source-data/geo_regions.tsv"), delimiter='\t')		# list of dicts
	country_to_region = {}
	for line in reader:
		country_to_region[line['country']] = line['region']
	for v in viruses:
		v['region'] = 'Unknown'
		if v['country'] in country_to_region:
			v['region'] = country_to_region[v['country']]
	
	return filter(lambda v: v['region'] != 'Unknown', viruses)

def get_virus_tuples(viruses):
	'''
	make dictionary of lists of viruses belonging to a certain date and region
	'''
	from collections import defaultdict
	virus_tuples = defaultdict(list)
	for v in viruses:
		vdate = datetime.datetime.strptime(v['date'], '%Y-%m-%d').date()
		virus_tuples[(vdate.year, vdate.month, v['region'])].append(v)

	return virus_tuples

def streamline(viruses, years_back, viruses_per_month):
	"""Subsample x viruses per month"""
	"""Take from beginning of list - this will prefer longer sequences"""
	"""Take viruses 1 per region in a cycle to get geographic diversity"""
	"""But pad with additional viruses from populous regions if necessary"""

	virus_tuples = get_virus_tuples(viruses)

	filtered_viruses = []
	first_year = datetime.datetime.today().year - years_back
	first_month = datetime.datetime.today().month
	regions = [v['region'] for v in viruses]
	regions = list(set(regions))

	print "Filtering between " + str(first_month) + "/" + str(first_year) + " and today"
	print "Selecting " + str(viruses_per_month) + " viruses per month"
	y = first_year
	for m in range(first_month,13):
		filtered_viruses.extend(select_viruses(virus_tuples, y, m, viruses_per_month, regions))
	for y in range(first_year+1,datetime.datetime.today().year+1):
		for m in range(1,13):
			filtered_viruses.extend(select_viruses(virus_tuples, y, m, viruses_per_month, regions))
	return filtered_viruses

def select_viruses(virus_tuples, y, m, viruses_per_month, regions):
	'''
	select viruses_per_month strains as evenly as possible from all regions
	'''
	from itertools import izip_longest
	select_set = []
	for representative in izip_longest(*[virus_tuples[(y,m,r)] for r in regions], fillvalue = None):
		select_set.extend([v in representative if v is not None])
	return select_set[:viruses_per_month]


def main(in_fname=None, years_back=3, viruses_per_month=50):

	print "--- Filter at " + time.strftime("%H:%M:%S") + " ---"

	if in_fname is None: in_fname = 'data/gisaid_epiflu_sequence.fasta'
	viruses = parse_gisaid(in_fname)
	print str(len(viruses)) + " initial viruses"

	# sort by sequence length
	viruses = sort_length(viruses)

	# fix strain names
	fix_strain_names(viruses)

	# add vaccine strains
	viruses = vaccine_strains + viruses
	print str(len(viruses)) + " with vaccine strains"

	# filter strain names
	viruses = filter_strain_names(viruses)
	print str(len(viruses)) + " with proper strain names"

	# filter short sequences
	viruses = filter_length(viruses)
	print str(len(viruses)) + " without truncation"

	# filter imprecise dates
	viruses = filter_date(viruses)
	print str(len(viruses)) + " with precise dates"

	# filter passage history
	viruses = filter_passage(viruses)
	print str(len(viruses)) + " without egg passage"

	# filter to unique strains
	viruses = filter_unique(viruses)
	print str(len(viruses)) + " with unique strain names"
	
	# append geo information
	viruses = append_country_and_region(viruses)
	print str(len(viruses)) + " with geographic information"	

	# reduce to manageable volume
	viruses = streamline(viruses, years_back, viruses_per_month)
	print str(len(viruses)) + " after streamlining"

	# add outgroup
	add_outgroup(viruses)
	print str(len(viruses)) + " with outgroup"
	out_fname = 'data/virus_filter.json'
	write_json(viruses, out_fname)
	return out_fname
	
if __name__ == "__main__":
	main()