# filter viruses after ingest, criteria based on metadata
#  - viruses longer than 987 bases
#  - viruses with exact dates
#  - viruses that are not egg-passaged
#  - a single sequence per virus strain, taken as first sequence in list
# outputs to virus_filter.json

import os, re, time, datetime
from io_util import *

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
	for y in range(2011,2015):
		for m in range(1,13):
			count = 0
			for v in viruses:
				date = datetime.datetime.strptime(v['date'], '%Y-%m-%d').date()
				if y == date.year and m == date.month:
					filtered_viruses.append(v)
					count += 1
					if count == 1:
						break
	return filtered_viruses
		
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