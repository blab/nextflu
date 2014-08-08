# filter viruses after ingest, criteria based on metadata
#  - viruses longer than 987 bases
#  - viruses with exact dates
#  - viruses that are not egg-passaged
#  - a single sequence per virus strain, taken as first sequence in list
# outputs to v_filter.json

import os, re, json

def read_viruses():
	try:
		handle = open('v_ingest.json', 'r')  
	except IOError:
		pass
	else:	
  		viruses = json.load(handle)
  		handle.close()
	return viruses

def filter_length(viruses):
	return filter(lambda v: len(v['nt']) >= 987, viruses)

def filter_date(viruses):
	return filter(lambda v: re.match(r'\d\d\d\d-\d\d-\d\d', v['date']) != None, viruses)
	
def filter_passage(viruses):
	round_one = filter(lambda v: re.match(r'E\d+', v.get('passage',''), re.I) == None, viruses)
	return filter(lambda v: re.match(r'Egg', v.get('passage',''), re.I) == None, round_one)
	
def add_outgroup(viruses):
	viruses.append({
		'strain': 'A/Beijing/32/1992',
		'db': 'ird',
		'accession': 'U26830',
		'date': '1990-01-01',
		'nt': 'ATGAAGACTATCATTGCTTTGAGCTACATTTTATGTCTGGTTTTCGCTCAAAAACTTCCCGGAAATGACAACAGCACAGCAACGCTGTGCCTGGGACATCATGCAGTGCCAAACGGAACGCTAGTGAAAACAATCACGAATGATCAAATTGAAGTGACTAATGCTACTGAGCTGGTTCAGAGTTCCTCAACAGGTAGAATATGCGACAGTCCTCACCGAATCCTTGATGGAAAAAACTGCACACTGATAGATGCTCTATTGGGAGACCCTCATTGTGATGGCTTCCAAAATAAGGAATGGGACCTTTTTGTTGAACGCAGCAAAGCTTACAGCAACTGTTACCCTTATGATGTACCGGATTATGCCTCCCTTAGGTCACTAGTTGCCTCATCAGGCACCCTGGAGTTTATCAATGAAGACTTCAATTGGACTGGAGTCGCTCAGGATGGGGGAAGCTATGCTTGCAAAAGGGGATCTGTTAACAGTTTCTTTAGTAGATTGAATTGGTTGCACAAATCAGAATACAAATATCCAGCGCTGAACGTGACTATGCCAAACAATGGCAAATTTGACAAATTGTACATTTGGGGGGTTCACCACCCGAGCACGGACAGAGACCAAACCAGCCTATATGTTCGAGCATCAGGGAGAGTCACAGTCTCTACCAAAAGAAGCCAACAAACTGTAACCCCGAATATCGGGTCTAGACCCTGGGTAAGGGGTCAGTCCAGTAGAATAAGCATCTATTGGACAATAGTAAAACCGGGAGACATACTTTTGATTAATAGCACAGGGAATCTAATTGCTCCTCGGGGTTACTTCAAAATACGAAATGGGAAAAGCTCAATAATGAGGTCAGATGCACCCATTGGCACCTGCAGTTCTGAATGCATCACTCCAAATGGAAGCATTCCCAATGACAAACCTTTTCAAAATGTAAACAGGATCACATATGGGGCCTGCCCCAGATATGTTAAGCAAAACACT'
	})

def filter_unique(viruses):
	filtered_viruses = []
	strains = set()	
	for v in viruses:
		if not v['strain'] in strains:
			strains.add(v['strain'])
			filtered_viruses.append(v)
	return filtered_viruses
	
def write_viruses(viruses):
	try:
		handle = open('v_filter.json', 'w') 
	except IOError:
		pass
	else:				
		json.dump(viruses, handle, indent=2)
		handle.close()		

def main():

	print "--- Virus filter ---"

	viruses = read_viruses()
	print str(len(viruses)) + " initial viruses"

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
	
	# add outgroup
	add_outgroup(viruses)
	print str(len(viruses)) + " with outgroup"
	
	write_viruses(viruses)

if __name__ == "__main__":
    main()