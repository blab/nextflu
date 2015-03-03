import argparse
import numpy as np
import dendropy
from collections import defaultdict
from date_util import numerical_date
from datetime import date
from itertools import izip
from tree_titer import *
from fitness_predictors import *

min_freq = 0.1
max_freq = 0.9
min_tips = 100
pc=1e-2
regularization = 1e-5

class fitness_model(object):

	def __init__(self, predictors = ['lb', 'ep', 'ne_star'], verbose=0,**kwargs):
		'''
		parameters:
		tree -- tree of sequences for which a fitness model is to be determined
		'''
		self.verbose=verbose
		self.seasons = [ (numerical_date(date(year=y, month = 10, day = 1)), 
							numerical_date(date(year = y+1, month = 4, day=1)))
						for y in xrange(int(self.time_interval[0]), int(self.time_interval[1]))]

		self.predictors = predictors
		self.predictor_functions = {'lb':self.calc_LBI, 
									'ne_star':self.calc_nonepitope_star_distance,
									'HI':self.calc_HI}
		self.LBI_tau = 0.005
		self.LBI_trans = lambda x:x**1.0
		self.current_season = self.seasons[-1]

	def calc_tip_counts(self):
		'''
		goes over all nodes and assigns tips to seasons (or any interval, really)
		then counts tips in every season for all internal nodes and calculates the
		fraction of viruses the clade defined by the internal node accounts for in 	
		a given season, i.e. the frequency 
		'''
		# count tips
		leaf_count = 0
		self.tip_aln = []
		for node in self.tree.postorder_node_iter():
			tmp_list = defaultdict(list)
			for child in node.child_nodes():
				for season, count in child.tips_by_season.iteritems():
					tmp_list[season].extend(count)

			node.tips_by_season = {}
			for season, strain_list in tmp_list.iteritems():
				node.tips_by_season[season] = np.array(strain_list, dtype=int)

			if node.is_leaf():
				tmp_season_list = [s for s in self.seasons if node.num_date>=s[0] and node.num_date<s[1]]
				if len(tmp_season_list)==1:
					node.tips_by_season[tmp_season_list[0]]=[leaf_count]
				leaf_count+=1
				self.tip_aln.append(np.fromstring(node.seq, 'S1'))
				node.tip_index = leaf_count

		self.tip_aln = np.array(self.tip_aln)
		# short cut to total number of tips per seaons
		total_counts = {s:len(strain_list) for s, strain_list in self.tree.seed_node.tips_by_season.iteritems()} 
		if self.verbose>1:
			for d,c in sorted(total_counts.items()): 
				print "number of tips in", d[0], '--', d[1],':',c

		# calculate frequencies
		for node in self.tree.postorder_node_iter():
			node.frequencies = defaultdict(float)
			for season, strain_list in node.tips_by_season.iteritems():
				node.frequencies[season]=float(len(strain_list))/(total_counts[season]+1e-10)

	def calc_predictors(self):
		print "calculating predictors"
		for pred in self.predictors:
			if pred in self.predictor_functions:
				self.predictor_functions[pred]()

	def select_nodes_in_season(self, season):
		self.current_season = season
		for node in self.tree.postorder_node_iter():
			if season in node.tips_by_season and len(node.tips_by_season[season])>0:
				node.alive=True
			else:
				node.alive=False

	def calc_all_predictors(self):
		self.predictor_arrays={}
		self.season_af = {}
		for node in self.tree.postorder_node_iter():
			node.predictors = {}
		for s in self.seasons:
			self.season_af[s] = np.zeros((len(self.nuc_alphabet), self.tip_aln.shape[1]))
			for ni, nuc in enumerate(self.nuc_alphabet):
				self.season_af[s][ni,:] = (self.tip_aln[self.tree.seed_node.tips_by_season[s],:]==nuc).mean(axis=0)

			tmp_preds = []
			t0=time.time()
			if self.verbose: print "calculating predictors for season ", s[0], '--', s[1],
			self.select_nodes_in_season(s)
			self.calc_predictors()
			if self.verbose: print np.round(time.time()-t0,2), 'seconds'
			for node in self.tree.postorder_node_iter():
				if node.alive:
					node.predictors[s] = np.array([node.__getattribute__(pred) 
				                              for pred in self.predictors])
					if node.is_leaf():
						tmp_preds.append(node.predictors[s])
				else:
					node.predictors[s]=None
					if node.is_leaf():
						tmp_preds.append(np.zeros(len(self.predictors)))
			self.predictor_arrays[s]=np.array(tmp_preds)

	def standardize_predictors(self):
		self.season_means = []
		self.season_std = []
		if self.verbose: print "standardize predictors for season"
		for s in self.seasons:
			self.season_means.append(self.predictor_arrays[s][self.tree.seed_node.tips_by_season[s],:].mean(axis=0))
			self.season_std.append(self.predictor_arrays[s][self.tree.seed_node.tips_by_season[s],:].std(axis=0))

		self.global_std = np.mean(self.season_std, axis=0)

		for s, m, stddev in izip(self.seasons, self.season_means, self.season_std):
			# keep this for internal nodes
			for node in self.tree.postorder_node_iter():
				if node.predictors[s] is not None:
					#node.predictors[s] = (node.predictors[s]-m)/self.global_std
					node.predictors[s] = (node.predictors[s]-m)/stddev

			self.predictor_arrays[s]-=m
			#self.predictor_arrays[s]/=self.global_std
			self.predictor_arrays[s]/=stddev

	def model_fit(self, params):
		# walk through season pairs
		seasonal_errors = []
		#import pdb; pdb.set_trace()
		for s,t in self.fit_test_season_pairs:		
			# normalize strain frequencies
			pred_af = self.fitness_biased_af(params, s)
			seasonal_errors.append(self.allele_frequency_distance(pred_af, self.season_af[t]))
		mean_error = np.mean(seasonal_errors)
		if any(np.isnan(seasonal_errors)+np.isinf(seasonal_errors)):
			mean_error = 1e10
		self.last_fit = mean_error
		if self.verbose>2: print params, self.last_fit
		return mean_error + regularization*np.sum(params**2)
		
	def fitness(self, params, pred):
		return np.sum(params*pred, axis=-1)

	def fitness_biased_af(self, params, season):
		pred_af = np.zeros_like(self.season_af[season])
		ind = self.tree.seed_node.tips_by_season[season]
		pred_freq = np.exp(self.fitness(params, self.predictor_arrays[season][ind,:]))
		for ni, nuc in enumerate(self.nuc_alphabet):
			pred_af[ni,:] = ((self.tip_aln[ind,:]==nuc).T*pred_freq).sum(axis=1)
		pred_af/=pred_af.sum(axis=0)
		return pred_af

	def allele_frequency_distance(self, af1, af2):
		return 1.0 - (af1*af2).sum(axis=0).mean(axis=-1)

	def minimize_error(self):
		from scipy.optimize import fmin_powell as minimizer
		if self.verbose:		
			print "initial function value:", self.model_fit(self.params)
			print "initial parameters:", self.params
		self.params = minimizer(self.model_fit, self.params, disp = self.verbose>1)
		if self.verbose:
			print "final function value:", self.model_fit(self.params)		
			print "final parameters:", self.params, '\n'		


	def learn_parameters(self, niter = 10):

		print "fitting parameters of the fitness model\n"

		params_stack = []

		if self.verbose:
			print "null parameters"
		self.params = 0*np.ones(len(self.predictors))  # initial values
		self.minimize_error()
		params_stack.append((self.last_fit, self.params))
		
		for ii in xrange(niter):
			if self.verbose:
				print "iteration:", ii+1
			self.params = 0.5*np.random.randn(len(self.predictors)) #0*np.ones(len(self.predictors))  # initial values
			self.minimize_error()
			params_stack.append((self.last_fit, self.params))

		self.params = params_stack[np.argmin([x[0] for x in params_stack])][1]
		self.model_fit(self.params)
		if self.verbose:
			print "best after",niter,"iterations\nfunction value:", self.last_fit
			print "fit parameters:"
			for pred, val in izip(self.predictors, self.params):
				print pred,':', val


	def assign_fitness(self, season):
		if self.verbose: print "calculating predictors for the last season"

		#FIXME: standardize predictors
		for node in self.tree.postorder_node_iter():
			if node.predictors[season] is not None:
				node.fitness = self.fitness(self.params, node.predictors[season])
			else:
				node.fitness = 0.0

	def predict(self, niter = 10):
		self.calc_tip_counts()
		self.calc_all_predictors()
		self.standardize_predictors()
		self.learn_parameters(niter = niter)
		self.assign_fitness(self.seasons[-1])

	######################################################################
	### fitness predictors
	######################################################################

	def calc_nonepitope_star_distance(self, attr='ne_star'):
		'''
		calculates the distance at nonepitope sites of any tree node to ref
		tree   --   dendropy tree
		attr   --   the attribute name used to save the result
		'''
		if not hasattr(self.tree, "nonepitope_star_distance_assigned") \
						or self.tree.nonepitope_star_distance_assigned==False:
			for node in self.tree.postorder_node_iter():
				if len(node.tips_by_season) and node!=self.tree.seed_node:
					tmp_node = node.parent_node
					cur_season = min(node.tips_by_season.keys())
					prev_season = self.seasons[max(0,self.seasons.index(cur_season)-1)]
					while True:
						if tmp_node!=self.tree.seed_node:
							if prev_season in tmp_node.tips_by_season and len(tmp_node.tips_by_season[prev_season])>0:
								break
							else:
								tmp_node=tmp_node.parent_node
						else:
							break
					node.__setattr__(attr, self.nonepitope_distance(node.aa_seq, tmp_node.aa_seq))
				else:
					node.__setattr__(attr, np.nan)				
			self.tree.nonepitope_star_distance_assigned=True



	def calc_LBI(self, attr = 'lb'):
		'''
		traverses the tree in postorder and preorder to calculate the
		up and downstream tree length exponentially weighted by distance.
		then adds them as LBI
		tree -- dendropy tree for whose node the LBI is being computed
		attr	 -- the attribute name used to store the result
		'''
		# traverse the tree in postorder (children first) to calculate msg to parents
		for node in self.tree.postorder_node_iter():
			node.down_polarizer = 0
			node.up_polarizer = 0
			for child in node.child_nodes():
				node.up_polarizer += child.up_polarizer
			bl =  node.edge_length/self.LBI_tau
			node.up_polarizer *= np.exp(-bl)
			if node.alive: node.up_polarizer += self.LBI_tau*(1-np.exp(-bl))

		# traverse the tree in preorder (parents first) to calculate msg to children
		for node in self.tree.preorder_internal_node_iter():
			for child1 in node.child_nodes():
				child1.down_polarizer = node.down_polarizer
				for child2 in node.child_nodes():
					if child1!=child2:
						child1.down_polarizer += child2.up_polarizer

				bl =  child1.edge_length/self.LBI_tau
				child1.down_polarizer *= np.exp(-bl)
				if child1.alive: child1.down_polarizer += self.LBI_tau*(1-np.exp(-bl))

		# go over all nodes and calculate the LBI (can be done in any order)
		for node in self.tree.postorder_node_iter():
			tmp_LBI = node.down_polarizer
			for child in node.child_nodes():
				tmp_LBI += child.up_polarizer
			node.__setattr__(attr, self.LBI_trans(tmp_LBI))

	def calc_HI(self, attr = 'HI'):
		print "estimating HI for season",self.current_season
		self.map_HI_to_tree(method='nnl1reg', lam=5, cutoff_date = self.current_season[1])
		for node in self.tree.postorder_node_iter():
			node.__setattr__(attr, node.cHI)

def test(params):
	from io_util import read_json
	from tree_util import json_to_dendropy, to_Biopython, color_BioTree_by_attribute
	from Bio import Phylo
	tree_fname='data/tree_refine_10y_50v.json'
	tree_fname='data/tree_refine.json'
	tree =  json_to_dendropy(read_json(tree_fname))
	fm = fitness_model(tree, predictors = params['predictors'], verbose=2)
	fm.predict(niter = params['niter'])
	#btree = to_Biopython(tree)
	#color_BioTree_by_attribute(btree, 'fitness')
	#Phylo.draw(btree, label_func=lambda x:'')
	return fm

def main(params):
	import time
	from io_util import read_json
	from io_util import write_json	
	from tree_util import json_to_dendropy, dendropy_to_json
	
	print "--- Start fitness model optimization at " + time.strftime("%H:%M:%S") + " ---"

	tree_fname='data/tree_refine.json'
	tree =  json_to_dendropy(read_json(tree_fname))
	fm = fitness_model(tree, predictors = params['predictors'], verbose=1)
	fm.predict(niter = params['niter'])
	out_fname = "data/tree_fitness.json"
	write_json(dendropy_to_json(tree.seed_node), out_fname)
	return out_fname

if __name__=="__main__":
	parser = argparse.ArgumentParser(description='Optimize predictor coefficients')
	parser.add_argument('-n', '--niter', type = int, default=10, help='number of replicate optimizations')
	parser.add_argument("-t", "--test", help="run test", action="store_true")
	parser.add_argument('-p', '--predictors', default=default_predictors, help='predictors to optimize', nargs='+')
	params = parser.parse_args().__dict__
	if params['test']:
		fm = test(params)
	else:
		main(params)

