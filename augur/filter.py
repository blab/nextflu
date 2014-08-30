# filter viruses after ingest, criteria based on metadata
#  - viruses longer than 987 bases
#  - viruses with exact dates
#  - viruses that are not egg-passaged
#  - a single sequence per virus strain, taken as first sequence in list
# outputs to virus_filter.json

import os, re, time, datetime
from io_util import *

YEARS_BACK = 3
VIRUSES_PER_MONTH = 100

def fix_strain_names(viruses):
	for v in viruses:
		v['strain'] = v['strain'].replace('\'','')
		v['strain'] = v['strain'].replace('(H3N2)','')
		
def filter_length(viruses):
	return filter(lambda v: len(v['seq']) >= 987, viruses)

def filter_date(viruses):
	return filter(lambda v: re.match(r'\d\d\d\d-\d\d-\d\d', v['date']) != None, viruses)
	
def filter_passage(viruses):
	round_one = filter(lambda v: re.match(r'^E\d+', v.get('passage',''), re.I) == None, viruses)
	return filter(lambda v: re.match(r'^Egg', v.get('passage',''), re.I) == None, round_one)
	
def add_outgroup(viruses):
	viruses.append({
		'strain': 'A/Beijing/32/1992',
		'db': 'IRD',
		'accession': 'U26830',
		'date': '1992-01-01',
		'seq': 'ATGAAGACTATCATTGCTTTGAGCTACATTTTATGTCTGGTTTTCGCTCAAAAACTTCCCGGAAATGACAACAGCACAGCAACGCTGTGCCTGGGACATCATGCAGTGCCAAACGGAACGCTAGTGAAAACAATCACGAATGATCAAATTGAAGTGACTAATGCTACTGAGCTGGTTCAGAGTTCCTCAACAGGTAGAATATGCGACAGTCCTCACCGAATCCTTGATGGAAAAAACTGCACACTGATAGATGCTCTATTGGGAGACCCTCATTGTGATGGCTTCCAAAATAAGGAATGGGACCTTTTTGTTGAACGCAGCAAAGCTTACAGCAACTGTTACCCTTATGATGTACCGGATTATGCCTCCCTTAGGTCACTAGTTGCCTCATCAGGCACCCTGGAGTTTATCAATGAAGACTTCAATTGGACTGGAGTCGCTCAGGATGGGGGAAGCTATGCTTGCAAAAGGGGATCTGTTAACAGTTTCTTTAGTAGATTGAATTGGTTGCACAAATCAGAATACAAATATCCAGCGCTGAACGTGACTATGCCAAACAATGGCAAATTTGACAAATTGTACATTTGGGGGGTTCACCACCCGAGCACGGACAGAGACCAAACCAGCCTATATGTTCGAGCATCAGGGAGAGTCACAGTCTCTACCAAAAGAAGCCAACAAACTGTAACCCCGAATATCGGGTCTAGACCCTGGGTAAGGGGTCAGTCCAGTAGAATAAGCATCTATTGGACAATAGTAAAACCGGGAGACATACTTTTGATTAATAGCACAGGGAATCTAATTGCTCCTCGGGGTTACTTCAAAATACGAAATGGGAAAAGCTCAATAATGAGGTCAGATGCACCCATTGGCACCTGCAGTTCTGAATGCATCACTCCAAATGGAAGCATTCCCAATGACAAACCTTTTCAAAATGTAAACAGGATCACATATGGGGCCTGCCCCAGATATGTTAAGCAAAACACT'
	})

def filter_unique(viruses):
	filtered_viruses = []
	strains = set()	
	for v in viruses:
		if not v['strain'].lower() in strains:
			strains.add(v['strain'].lower())
			filtered_viruses.append(v)
	return filtered_viruses
	
def streamline(viruses):
	filtered_viruses = []
	first_year = datetime.datetime.today().year - YEARS_BACK
	first_month = datetime.datetime.today().month
	print "Filtering between " + str(first_month) + "/" + str(first_year) + " and today"
	print "Selecting " + str(VIRUSES_PER_MONTH) + " viruses per month"
	y = first_year
	for m in range(first_month,13):
		filtered_viruses.extend(select_viruses(viruses, y, m))	
	for y in range(first_year+1,datetime.datetime.today().year+1):
		for m in range(1,13):
			filtered_viruses.extend(select_viruses(viruses, y, m))
	return filtered_viruses
	
def select_viruses(viruses, y, m):
	count = 0
	filtered = []
	for v in viruses:
		date = datetime.datetime.strptime(v['date'], '%Y-%m-%d').date()
		if y == date.year and m == date.month:
			filtered.append(v)
			count += 1
			if count == VIRUSES_PER_MONTH:
				break
	return filtered

def main():

	print "--- Filter at " + time.strftime("%H:%M:%S") + " ---"

	viruses = read_json('virus_ingest.json')
	print str(len(viruses)) + " initial viruses"
	
	# fix strain names
	fix_strain_names(viruses)

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
	
	# reduce to manageable volume
	viruses = streamline(viruses)
	print str(len(viruses)) + " after streamlining"	
	
	# add outgroup
	add_outgroup(viruses)
	print str(len(viruses)) + " with outgroup"
	
	write_json(viruses, 'virus_filter.json')

if __name__ == "__main__":
    main()