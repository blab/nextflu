import numpy as np
from Bio import Seq, AlignIO
import time
from io_util import read_json
from io_util import write_json
from tree_util import json_to_dendropy, dendropy_to_json
from itertools import izip

import dendropy

def load_mutational_tolerance():
	fname = 'data/MD_H1_combined_prefs_noStart_noStop_rescaled.txt'
	with open(fname) as f:
		aa = map(lambda x:x.split('_')[1], f.readline().strip().split()[3:])
	sites = np.loadtxt(fname, usecols=[0], dtype=int)
	wt_aa = np.loadtxt(fname, usecols=[1], dtype='S1')
	aa_prob = np.loadtxt(fname, usecols=range(3,23), dtype=float)
	return aa, sites, wt_aa, aa_prob

def calc_fitness_tolerance(aa_seq, aa_prob, aa, indices, beta = 1.0):
	'''
	determine the indices of aligned amino acids and sum the logged probabilities
	'''
	H3_aa_indices = np.array([aa.index(aa_seq[p]) if aa_seq[p] in aa else -1 for p in indices])
	return np.sum((H3_aa_indices!=-1)*np.log(aa_prob[(np.arange(len(indices)), H3_aa_indices)])*beta)

def assign_fitness_tolerance(nodes, epitope_mask=None, attr='tol'):
	'''
	loops over all viruses, translates their sequences and calculates the virus fitness
	takes a dendropy tree or a list of nodes
	'''
	aa, sites, wt_aa, aa_prob = load_mutational_tolerance()
	aln = AlignIO.read('source-data/H1_H3_v2.fasta', 'fasta')
	# returns true whenever neither of the sequences have a gap
	aligned = (np.array(aln)!='-').min(axis=0)
	# map alignment positions to sequence positions, subset to aligned amino acids
	indices = {}
	for seq in aln:
		indices[seq.name] = (np.cumsum(np.fromstring(str(seq.seq), dtype='S1')!='-')-1)[aligned]

	# make a reduced set of amino-acid probabilities that only contains aligned positions
	# column ii in aa_prob[ii,:] now corresponds to position indices['H3'][ii] in the H3 HA
	aa_prob=aa_prob[indices['H1'],:]
	# attach another column for non-canonical amino acids
	aa_prob = np.hstack((aa_prob, 1e-5*np.ones((aa_prob.shape[0],1))))
	print("AA preference matrix:",aa_prob.shape)
	if epitope_mask is not None: # strip epitope positions
	#	nonepi_pos = np.where(epitope_mask=='0')[0]
	#	nonepi_indices = np.in1d(indices['H3'], nonepi_pos)
		# determine those position in the H3 index set that correspond to nonepitope positions
		# nonepi_indices[ii]==True if indices['H3'][ii] is a non epitope position		
		nonepi_indices = []
		for pos in indices['H3']:
			if epitope_mask[pos][0] == '0':
				nonepi_indices.append(True)
			if epitope_mask[pos][0] == '1':
				nonepi_indices.append(False)								
		# remove epitope positions from the aa matrix and the indices				
		aa_prob=aa_prob[nonepi_indices,:]
	#	indices['H3'] = indices['H3'][nonepi_indices]
		indices['H3'] = [x for x, y in izip(indices['H3'], nonepi_indices) if y == True]
		# column ii in aa_prob[ii,:] again corresponds to position indices['H3'][ii] in the H3 HA
		# but only nonepitope indices are in indices['H3']
		#assert all(epitope_mask[indices['H3']]=='0')
	if isinstance(nodes, list):
		for node in nodes:
			tolerance = calc_fitness_tolerance(Seq.translate(node['seq']), aa_prob, aa, indices['H3'])
			node.__setattr__(attr, tolerance)
	elif isinstance(nodes, dendropy.Tree):
		for node in nodes:
			tolerance = calc_fitness_tolerance(node.aa_seq['SigPep']+node.aa_seq['HA1']+node.aa_seq['HA2'], aa_prob, aa, indices['H3'])
			node.__setattr__(attr, tolerance)

def main(in_fname='data/tree_refine.json', tree=True):

	print "--- Mutational tolerance at " + time.strftime("%H:%M:%S") + " ---"
	viruses = read_json(in_fname)
	if tree:
		viruses = json_to_dendropy(viruses)

	assign_fitness_tolerance(viruses)

	if tree:
		out_fname = "data/tree_tolerance.json"
		write_json(dendropy_to_json(viruses.seed_node), out_fname)
	else:
		out_fname = "data/virus_tolerance.json"
		write_json(viruses, out_fname)
	return out_fname, viruses

if __name__=='__main__':
	out_fname, viruses = main()
