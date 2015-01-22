# assign epitope fitness to each node in the phylogeny

import time
from io_util import *
from tree_util import *
from date_util import *
from seq_util import *
import numpy as np
from itertools import izip
from collections import defaultdict

def append_nonepitope_sites(viruses):
	for virus in viruses:
		sites_ne = nonepitope_sites(translate(virus['seq']))
		virus['sites_ne'] = sites_ne
		
def remove_nonepitope_sites(viruses):
	for virus in viruses:
		virus.pop("sites_ne", None)
		
def remove_nonepitope_distances(viruses):
	for virus in viruses:
		virus.pop("distance_ne", None)			
		
def most_frequent(char_list):
	d = defaultdict(int)
	for i in char_list:
		d[i] += 1
	return sorted(d.iteritems(), key=lambda x: x[1], reverse=True)[0][0]
		
def consensus_nonepitope(viruses):
	"""Return consensus non-epitope sequence"""
	consensus = ""
	length = len(viruses[0]['sites_ne'])
	for i in range(0, length):
		column = [v['sites_ne'][i] for v in viruses]
		consensus += most_frequent(column)
	return consensus
	
def distance_to_consensus(virus, consensus_ne):
	"""Return distance of virusA to virusB by comparing non-epitope sites"""
	virus_ne = virus['sites_ne']
	ne_distance = sum(a != b for a, b in izip(virus_ne, consensus_ne))
	return ne_distance
	
def compute(viruses):
	"""Append non-epitope distances to each virus"""
	print "Computing epitope distances"
	consensus = consensus_nonepitope(viruses)
	for virus in viruses:
		distance = distance_to_consensus(virus, consensus)
		virus['distance_ne'] = distance
		print virus['strain'] + ": " + str(virus['distance_ne'])
	
def normalize(viruses):
	"""Normalizing non-epitope distances to give non-epitope fitness"""
	print "Normalizing non-epitope distances"
	distances = [v['distance_ne'] for v in viruses]
	mean = np.mean(distances)
	sd = np.std(distances)
	for virus in viruses:
		virus['fitness_ne'] = -1 * ( ( virus['distance_ne'] - mean) / sd )
		print virus['strain'] + ": " + str(virus['fitness_ne'])

def main():
	print "--- Non-epitope fitness at " + time.strftime("%H:%M:%S") + " ---"

	viruses = read_json('data/virus_epitope.json')
				
	append_nonepitope_sites(viruses)
	compute(viruses)
#	normalize(viruses)		
	remove_nonepitope_sites(viruses)
#	remove_nonepitope_distances(viruses)
	
	write_json(viruses, "data/virus_nonepitope.json")	
	
if __name__ == "__main__":
    main()
