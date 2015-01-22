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

	tree = read_json('data/tree_clean.json')
	nodes = [n for n in all_descendants(tree)]
	
	tips = []
	for node in nodes:
		if 'strain' in node: 
			tips.append(node)		
			
	append_epitope_sites(tips)

	for tip in tips:
		distance = population_distance(tip, tips, IMMUNITY_SCALING)
		tip['fitness_ep'] = round(distance, 5)
		
	remove_epitope_sites(tips)	
	write_json(tree, "data/tree_epitope.json")	
	
if __name__ == "__main__":
    main()
