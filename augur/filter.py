# after ingesting sequence data, run this to remove extra sequences, leaving only
#  - viruses longer than 900 bases
#  - viruses with exact dates
#  - viruses that are not egg-passaged
#  - a single sequence per virus strain, taken as first sequence in list

import os, re, json

def filter_length(viruses):
	return filter(lambda v: len(v['nt']) > 900, viruses)

def filter_date(viruses):
	return filter(lambda v: re.match(r'\d\d\d\d-\d\d-\d\d', v['date']) != None, viruses)
	
def filter_passage(viruses):
	round_one = filter(lambda v: re.match(r'E\d+', v.get('passage',''), re.I) == None, viruses)
	return filter(lambda v: re.match(r'Egg', v.get('passage',''), re.I) == None, round_one)

def filter_unique(viruses):
	filtered_viruses = []
	strains = set()	
	for v in viruses:
		if not v['strain'] in strains:
			strains.add(v['strain'])
			filtered_viruses.append(v)
	return filtered_viruses

def main():

	print "--- Virus filter ---"

	try:
		handle = open('v_ingest.json', 'r')  
	except IOError:
		pass
	else:	
  		viruses = json.load(handle)
  		handle.close()
	print str(len(viruses)) + " initial viruses"

	# filter short sequences
	viruses = filter_length(viruses)
	print str(len(viruses)) + " with full HA1"
	
	# filter imprecise dates
	viruses = filter_date(viruses)
	print str(len(viruses)) + " with precise dates"
	
	# filter passage history
	viruses = filter_passage(viruses)
	print str(len(viruses)) + " without egg passage"
	
	# filter to unique strains
	viruses = filter_unique(viruses)
	print str(len(viruses)) + " with unique strain names"
	
	try:
		handle = open('v_filter.json', 'w') 
	except IOError:
		pass
	else:				
		json.dump(viruses, handle, indent=2)
		handle.close()	

if __name__ == "__main__":
    main()