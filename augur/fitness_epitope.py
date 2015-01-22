# assign epitope fitness to each node in the phylogeny

import time
from io_util import *
from tree_util import *
from date_util import *
from seq_util import *
import numpy
from itertools import izip

IMMUNITY_SCALING = 0.023

def append_epitope_sites(viruses):
	for virus in viruses:
		epi = epitope_sites(translate(virus['seq']))
		virus['epi'] = epi
		
def remove_epitope_sites(viruses):
	for virus in viruses:
		virus.pop("epi", None)	
		
def pairwise_distance(virusA, virusB, scaling):
	"""Return distance of virusA to virusB by comparing epitope sites"""
	"""Scaling is decrease in cross-immunity per epitope amino acid substitution"""
	"""Default value for scaling is 0.023"""
	epiA = virusA['epi']
	epiB = virusB['epi']
	epi_distance = sum(a != b for a, b in izip(epiA, epiB))
	return epi_distance * scaling

def population_distance(virus, population, scaling):
	distances = [pairwise_distance(virus, v, scaling) for v in population]
	return numpy.mean(distances)

def main():
	print "--- Epitope fitness at " + time.strftime("%H:%M:%S") + " ---"

	viruses = read_json('data/virus_clean.json')
				
	append_epitope_sites(viruses)

	for virus in viruses:
		distance = population_distance(virus, viruses, IMMUNITY_SCALING)
		virus['fitness_ep'] = distance
		
	remove_epitope_sites(viruses)	
	write_json(viruses, "data/virus_epitope.json")	
	
if __name__ == "__main__":
    main()
