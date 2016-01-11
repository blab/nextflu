import dendropy, time
import numpy as np
from itertools import izip
from scipy.stats import linregress
from seq_util import translate
from io_util import read_json
from io_util import write_json
from tree_util import json_to_dendropy
from tree_util import dendropy_to_json
from fitness_tolerance import load_mutational_tolerance, calc_fitness_tolerance

class fitness_predictors(object):

	def __init__(self, predictor_names = ['ep', 'lb', 'dfreq'], **kwargs):
		if "epitope_masks_fname" in kwargs and "epitope_mask_version" in kwargs:
			self.setup_epitope_mask(epitope_masks_fname = kwargs["epitope_masks_fname"], epitope_mask_version = kwargs["epitope_mask_version"])
		else:
			self.setup_epitope_mask()
		self.predictor_names = predictor_names
	
	def setup_predictor(self, tree, pred, timepoint):
		if pred == 'lb':
			self.calc_LBI(tree, tau = 0.0005, transform = lambda x:x)
		if pred == 'ep':
			self.calc_epitope_distance(tree)
		if pred == 'epH':
			self.calc_epitope_history(tree, timepoint)			
		if pred == 'ne':
			self.calc_nonepitope_distance(tree)
		if pred == 'ne_star':
			self.calc_nonepitope_star_distance(tree)
		if pred == 'tol':
			self.calc_tolerance(tree)
		#if pred == 'dfreq':
			# do nothing
		#if pred == 'cHI':
			# do nothing

	def setup_epitope_mask(self, epitope_masks_fname = 'source-data/H3N2_epitope_masks.tsv', epitope_mask_version = 'wolf'):
		self.epitope_mask = ""
		epitope_map = {}
		with open(epitope_masks_fname) as f:
			for line in f:
				(key, value) = line.split()
				epitope_map[key] = value
		if epitope_mask_version in epitope_map:
			self.epitope_mask = epitope_map[epitope_mask_version]

	def epitope_sites(self, aa):
		sites = []
		for a, m in izip(aa, self.epitope_mask):
			if m == '1':
				sites.append(a)
		return ''.join(sites)

	def nonepitope_sites(self, aa):
		sites = []
		for a, m in izip(aa, self.epitope_mask):
			if m == '0':
				sites.append(a)
		return ''.join(sites)

	def receptor_binding_sites(self, aa):
		sp = 16
		aaa = np.fromstring(aa, 'S1')
		receptor_binding_list = map(lambda x:x+sp-1, [145, 155, 156, 158, 159, 189, 193])
		return ''.join(aaa[receptor_binding_list])	

	def epitope_distance(self, aaA, aaB):
		"""Return distance of sequences aaA and aaB by comparing epitope sites"""
		epA = self.epitope_sites(aaA)
		epB = self.epitope_sites(aaB)
		distance = sum(a != b for a, b in izip(epA, epB))
		return distance

	def nonepitope_distance(self, aaA, aaB):
		"""Return distance of sequences aaA and aaB by comparing non-epitope sites"""
		neA = self.nonepitope_sites(aaA)
		neB = self.nonepitope_sites(aaB)
		distance = sum(a != b for a, b in izip(neA, neB))
		return distance

	def rbs_distance(self, aaA, aaB):
		"""Return distance of sequences aaA and aaB by comparing receptor binding sites (Koel sites)"""
		rbsA = self.receptor_binding_sites(aaA)
		rbsB = self.receptor_binding_sites(aaB)
		distance = sum(a != b for a, b in izip(rbsA, rbsB))
		return distance

	def calc_epitope_distance(self, tree, attr='ep', ref = None):
		'''
		calculates the distance at epitope sites of any tree node  to ref
		tree   --   dendropy tree
		attr   --   the attribute name used to save the result
		'''
		if ref == None:
			ref = translate(tree.seed_node.seq)
		for node in tree.postorder_node_iter():
			if not hasattr(node, 'aa'):
				node.aa = translate(node.seq)
			node.__setattr__(attr, self.epitope_distance(node.aa, ref))
			
	def calc_epitope_history(self, tree, timepoint, attr='epH'):
		'''
		calculates the distance at epitope sites of any tree node to all nodes before timepoint
		these comparison nodes are a proxy for the host immune landscape at that time
		tree   --   dendropy tree
		attr   --   the attribute name used to save the result
		'''
		comparison_nodes = [] 
		for node in tree.postorder_node_iter():
			if not hasattr(node, 'aa'):
				node.aa = translate(node.seq)
			node.__setattr__(attr, 0)
			if node.is_leaf():
				if node.num_date <= timepoint:
					comparison_nodes.append(node)
		for node in tree.postorder_node_iter():
			if node.is_leaf():
				mean_distance = 0
				count = 0
				for comp_node in comparison_nodes:
					mean_distance += self.epitope_distance(node.aa, comp_node.aa)
					count += 1
				epitope_history = mean_distance / float(count)
				node.__setattr__(attr, epitope_history)			

	def calc_rbs_distance(self, tree, attr='rb', ref = None):
		'''
		calculates the distance at receptor binding sites of any tree node to ref
		tree   --   dendropy tree
		attr   --   the attribute name used to save the result
		'''
		if ref == None:
			ref = translate(tree.seed_node.seq)
		for node in tree.postorder_node_iter():
			if not hasattr(node, 'aa'):
				node.aa = translate(node.seq)
			node.__setattr__(attr, self.rbs_distance(node.aa, ref))

	def calc_tolerance(self, tree, attr='tol'):
		'''
		calculates the distance at epitope sites of any tree node  to ref
		tree   --   dendropy tree
		attr   --   the attribute name used to save the result
		'''
		from Bio import AlignIO
		aa, sites, wt_aa, aa_prob = load_mutational_tolerance()
		aln = AlignIO.read('source-data/H1_H3.fasta', 'fasta')
		# returns true whenever either of the sequences have a gap
		aligned = (np.array(aln)!='-').min(axis=0)
		# map alignment positions to sequence positions, subset to aligned amino acids
		indices = {}
		for seq in aln:
			indices[seq.name] = (np.cumsum(np.fromstring(str(seq.seq), dtype='S1')!='-')-1)[aligned]

		# make a reduced set of amino-acid probabilities that only contains aligned positions
		aa_prob=aa_prob[indices['H1'],:]
		# attach another column for non-canonical amino acids

		aa_prob = np.hstack((aa_prob, 1e-5*np.ones((aa_prob.shape[0],1))))	
		
		for node in tree.postorder_node_iter():
			if not hasattr(node, 'aa'):
				node.aa = translate(node.seq)
			node.__setattr__(attr, calc_fitness_tolerance(node.aa, aa_prob, aa, indices['H3']))


	def calc_nonepitope_distance(self, tree, attr='ne', ref = None):
		'''
		calculates the distance at nonepitope sites of any tree node to ref
		tree   --   dendropy tree
		attr   --   the attribute name used to save the result
		'''
		if ref == None:
			ref = translate(tree.seed_node.seq)
		for node in tree.postorder_node_iter():
			if not hasattr(node, 'aa'):
				node.aa = translate(node.seq)
			node.__setattr__(attr, self.nonepitope_distance(node.aa, ref))

	def calc_nonepitope_star_distance(self, tree, attr='ne_star', seasons = []):
		'''
		calculates the distance at nonepitope sites of any tree node to ref
		tree   --   dendropy tree
		attr   --   the attribute name used to save the result
		'''
		for node in tree.postorder_node_iter():
			if len(node.season_tips) and node!=tree.seed_node:
				if not hasattr(node, 'aa'):
					node.aa = translate(node.seq)
				tmp_node = node.parent_node
				cur_season = min(node.season_tips.keys())
				prev_season = seasons[max(0,seasons.index(cur_season)-1)]
				while True:
					if tmp_node!=tree.seed_node:
						if prev_season in tmp_node.season_tips and len(tmp_node.season_tips[prev_season])>0:
							break
						else:
							tmp_node=tmp_node.parent_node
					else:
						break
				if not hasattr(tmp_node, 'aa'):
					tmp_node.aa = translate(tmp_node.seq)
				node.__setattr__(attr, self.nonepitope_distance(node.aa, tmp_node.aa))
			else:
				node.__setattr__(attr, np.nan)				


	def calc_LBI(self, tree, attr = 'lb', tau=0.0005, transform = lambda x:x):
		'''
		traverses the tree in postorder and preorder to calculate the
		up and downstream tree length exponentially weighted by distance.
		then adds them as LBI
		tree -- dendropy tree for whose node the LBI is being computed
		attr	 -- the attribute name used to store the result
		'''
		# traverse the tree in postorder (children first) to calculate msg to parents
		for node in tree.postorder_node_iter():
			node.down_polarizer = 0
			node.up_polarizer = 0
			for child in node.child_nodes():
				node.up_polarizer += child.up_polarizer
			bl =  node.edge_length/tau
			node.up_polarizer *= np.exp(-bl)
			if node.alive: node.up_polarizer += tau*(1-np.exp(-bl))

		# traverse the tree in preorder (parents first) to calculate msg to children
		for node in tree.preorder_internal_node_iter():
			for child1 in node.child_nodes():
				child1.down_polarizer = node.down_polarizer
				for child2 in node.child_nodes():
					if child1!=child2:
						child1.down_polarizer += child2.up_polarizer

				bl =  child1.edge_length/tau
				child1.down_polarizer *= np.exp(-bl)
				if child1.alive: child1.down_polarizer += tau*(1-np.exp(-bl))

		# go over all nodes and calculate the LBI (can be done in any order)
		for node in tree.postorder_node_iter():
			tmp_LBI = node.down_polarizer
			for child in node.child_nodes():
				tmp_LBI += child.up_polarizer
			node.__setattr__(attr, transform(tmp_LBI))

def main(tree_fname = 'data/tree_refine.json'):

	print "--- Testing predictor evaluations ---"
	tree =  json_to_dendropy(read_json(tree_fname))

	print "Calculating epitope distances"
	calc_epitope_distance(tree)

	print "Calculating nonepitope distances"
	calc_nonepitope_distance(tree)

	print "Calculating LBI"
#	calc_LBI(tree)

	print "Writing decorated tree"
	out_fname = "data/tree_predictors.json"
	write_json(dendropy_to_json(tree.seed_node), out_fname)
	return out_fname

if __name__=='__main__':
	main()
