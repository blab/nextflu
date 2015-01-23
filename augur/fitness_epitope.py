# assign epitope fitness to each node in the phylogeny

import time
from io_util import *
from tree_util import *
from date_util import *
from seq_util import *
import numpy as np
from itertools import izip

def append_epitope_sites(viruses):
	for virus in viruses:
		sites_ep = epitope_sites(virus['seq'])
		virus['sites_ep'] = sites_ep
		
def remove_epitope_sites(viruses):
	for virus in viruses:
		virus.pop("sites_ep", None)
		
def remove_epitope_distances(viruses):
	for virus in viruses:
		virus.pop("distance_ep", None)			
		
def pairwise_distance(virusA, virusB):
	"""Return distance of virusA to virusB by comparing epitope sites"""
	epA = virusA['sites_ep']
	epB = virusB['sites_ep']
	ep_distance = sum(a != b for a, b in izip(epA, epB))
	return ep_distance

def population_distance(virus, population):
	distances = [pairwise_distance(virus, v) for v in population]
	return np.mean(distances)
	
def compute(viruses):
	"""Append epitope distances to each virus"""
	print "Computing epitope distances"
	for virus in viruses:
		distance = population_distance(virus, viruses)
		virus['distance_ep'] = distance
		print virus['strain'] + ": " + str(virus['distance_ep'])
	
def normalize(viruses):
	"""Normalizing epitope distances to give epitope fitness"""
	print "Normalizing epitope distances"
	distances = [v['distance_ep'] for v in viruses]
	mean = np.mean(distances)
	sd = np.std(distances)
	for virus in viruses:
		virus['fitness_ep'] = (virus['distance_ep'] - mean) / sd
		print virus['strain'] + ": " + str(virus['fitness_ep'])

def main():
	print "--- Epitope fitness at " + time.strftime("%H:%M:%S") + " ---"

	viruses = read_json('data/virus_clean.json')
					
	append_epitope_sites(viruses)
	compute(viruses)
#	normalize(viruses)		
	remove_epitope_sites(viruses)	
#	remove_epitope_distances(viruses)
	
	write_json(viruses, "data/virus_epitope.json")	
	
if __name__ == "__main__":
    main()
