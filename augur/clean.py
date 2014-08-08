# clean sequences after alignment, criteria based on sequences

import os, re, json, datetime
from scipy import stats
import numpy

OUTGROUP = 'A/Beijing/32/1992'

def read_viruses():
	try:
		handle = open('v_align.json', 'r')  
	except IOError:
		pass
	else:	
  		viruses = json.load(handle)
  		handle.close()
	return viruses
	
def clean_gaps(viruses):
	return [v for v in viruses if not '-' in v['nt']]

def date_from_virus(virus):
	return datetime.datetime.strptime(virus['date'], '%Y-%m-%d').date()
	
def times_from_outgroup(viruses):
	outgroup_date = date_from_virus(filter(lambda v: v['strain'] == OUTGROUP, viruses)[0])
	return map(lambda v: (date_from_virus(v) - outgroup_date).days, viruses)
	
def distance(seq_a, seq_b):
	dist = 0
	for (a, b) in zip(list(seq_a), list(seq_b)):
		if a != b:
			dist += 1
	return dist

def distances_from_outgroup(viruses):
	outgroup_seq = filter(lambda v: v['strain'] == OUTGROUP, viruses)[0]['nt']
	return map(lambda v: distance(outgroup_seq, v['nt']), viruses)

def clean_distances(viruses):
	times = times_from_outgroup(viruses)
	distances = distances_from_outgroup(viruses)
	slope, intercept, r_value, p_value, std_err = stats.linregress(times, distances)
	print "  slope: " + str(slope*365)
	print "  r: " + str(r_value)
	residuals = map(lambda t, d: (intercept + slope*t) - d, times, distances)
	r_sd = numpy.std(residuals)	
	print "  residuals sd: " + str(r_sd)	
	new_viruses = []
	for (v,r) in zip(viruses,residuals):		# filter viruses more than 5 sds up or down
		if r > -5 * r_sd and r < 5 * r_sd:
			new_viruses.append(v)
	return new_viruses			
		
	
def write_viruses(viruses):
	try:
		handle = open('v_clean.json', 'w') 
	except IOError:
		pass
	else:				
		json.dump(viruses, handle, indent=2)
		handle.close()			

def main():

	print "--- Virus clean ---"
	
	viruses = read_viruses()
	print str(len(viruses)) + " initial viruses"
	
	# clean gapped sequences
	viruses = clean_gaps(viruses)
	print str(len(viruses)) + " with complete HA1"
	
	# clean sequences by distance	
	viruses = clean_distances(viruses)
	print str(len(viruses)) + " with clock"	
	
	write_viruses(viruses)	

if __name__ == "__main__":
    main()