import argparse
import numpy as np
import dendropy
from collections import defaultdict
from datetime import date
from itertools import izip
from fitness_predictors import *

ymin = 2005
ymax = 2015
min_freq = 0.1
max_freq = 0.9
min_tips = 10
pc=1e-2
default_predictors = ['lb', 'ep', 'ne_start']

class fitness_model(object):

	def __init__(self, tree, predictors = None, verbose=0):
		'''
		parameters:
		tree -- tree of sequences for which a fitness model is to be determined
		'''
		self.tree = tree
		self.verbose=verbose
		self.seasons = [ (date(year=y, month = 10, day = 1), date(year = y+1, month = 4, day=1)) 
						for y in xrange(ymin, ymax)]

		self.predictors = []
		for p in predictors:
			if p == 'lb':
				self.predictors.append(('lb',calc_LBI, {'tau':0.0005, 'transform':lambda x:x}))
			if p == 'ep':
				self.predictors.append(('ep',calc_epitope_distance,{}))
			if p == 'ne':
				self.predictors.append(('ne',calc_nonepitope_distance,{}))
			if p == 'ne_star':
				self.predictors.append(('ne_star',calc_nonepitope_star_distance,{"seasons":self.seasons}))				
			if p == 'tol':
				self.predictors.append(('tol',calc_tolerance,{}))				

	def calc_tip_counts(self):
		'''
		goes over all nodes and assigns tips to seasons (or any interval, really)
		then counts tips in every season for all internal nodes and calculates the
		fraction of viruses the clade defined by the internal node accounts for in 	
		a given season, i.e. the frequency 
		'''
		# count tips
		leaf_count = 0
		for node in self.tree.postorder_node_iter():
			tmp_list = defaultdict(list)
			for child in node.child_nodes():
				for season, count in child.tips.iteritems():
					tmp_list[season].extend(count)

			node.tips = {}
			for season, strain_list in tmp_list.iteritems():
				node.tips[season] = np.array(strain_list, dtype=int)

			if node.is_leaf():
				node_date = date(*map(int, node.date.split('-')))
				tmp_season_list = [s for s in self.seasons if node_date>=s[0] and node_date<s[1]]
				if len(tmp_season_list)==1:
					node.tips[tmp_season_list[0]]=[leaf_count]
				leaf_count+=1
				node.tip_index = leaf_count

		# short cut to total number of tips per seaons
		total_counts = {s:len(strain_list) for s, strain_list in self.tree.seed_node.tips.iteritems()} 
		if self.verbose>1:
			for d,c in sorted(total_counts.items()): 
				print "number of tips in", d[0].strftime('%Y-%m-%d'), '--', \
											d[1].strftime('%Y-%m-%d'),':',c

		# calculate frequencies
		for node in self.tree.postorder_node_iter():
			node.frequencies = defaultdict(float)
			for season, strain_list in node.tips.iteritems():
				if season in total_counts:
					node.frequencies[season]=float(len(strain_list))/(total_counts[season]+1e-10)
				else:
					node.frequencies[season]=0


	def calc_predictors(self):
		for pred, func, kwargs in self.predictors:
			# calculate the predictors for all nodes of the tree and save as node.attr
			func(self.tree, attr = pred, **kwargs)

	def select_nodes_in_season(self, season):
		for node in self.tree.postorder_node_iter():
			if season in node.tips and len(node.tips[season])>0:
				node.alive=True
			else:
				node.alive=False

	def calc_all_predictors(self):
		self.predictor_arrays={}
		for node in self.tree.postorder_node_iter():
			node.predictors = {}
		for s in self.seasons:
			tmp_preds = []
			t0=time.time()
			if self.verbose: print "calculating predictors for season ", s[0].strftime("%Y-%m-%d"), '--', s[1].strftime("%Y-%m-%d"),
			self.select_nodes_in_season(s)
			self.calc_predictors()
			if self.verbose: print np.round(time.time()-t0,2), 'seconds'
			for node in self.tree.postorder_node_iter():
				if node.alive:
					node.predictors[s] = np.array([node.__getattribute__(pred[0]) 
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
			self.season_means.append(self.predictor_arrays[s][self.tree.seed_node.tips[s],:].mean(axis=0))
			self.season_std.append(self.predictor_arrays[s][self.tree.seed_node.tips[s],:].std(axis=0))

		self.global_std = np.mean(self.season_std, axis=0)

		for s, m in izip(self.seasons, self.season_means):
			# keep this for internal nodes
			for node in self.tree.postorder_node_iter():
				if node.predictors[s] is not None:
					node.predictors[s] = (node.predictors[s]-m)/self.global_std

			self.predictor_arrays[s]-=m
			self.predictor_arrays[s]/=self.global_std


	def select_clades_for_fitting(self):
		# short cut to total number of tips per seaons
		total_counts = {s:len(strain_list) for s, strain_list in self.tree.seed_node.tips.iteritems()}
		# prune seasons where few observations were made, only consecutive pairs
		# with sufficient tip count are retained
		self.fit_test_season_pairs = [(s,t) for s,t in izip(self.seasons[:-1], self.seasons[1:]) 
							if total_counts[s]>min_tips and total_counts[t]>min_tips]
		
		# for each pair of seasons with sufficient tip counts, 
		# select clades in a certain frequency window
		# keep track of these season / clade combinations in the dict clades_for_season
		self.clades_for_season = {}
		for s,t in self.fit_test_season_pairs:
			self.clades_for_season[(s,t)] = []
			for node in self.tree.postorder_node_iter():
				if node.frequencies[s]>=min_freq and node.frequencies[s]<max_freq:
					if max([c.frequencies[s] for c in node.child_nodes()])<node.frequencies[s]:
						self.clades_for_season[(s,t)].append(node)
					else:
						if self.verbose: print "clade fully contained in daughter clade"


	def model_fit(self, params):
		# walk through season pairs
		seasonal_errors = []
		for s,t in self.fit_test_season_pairs:		
			# normalize strain frequencies
			total_strain_freq = np.exp(self.fitness(params, self.predictor_arrays[s][self.tree.seed_node.tips[s],:])).sum()
		
			# project clades forward according to strain makeup
			clade_errors = []
			test_clades = self.clades_for_season[(s,t)]
			for clade in test_clades:
				initial_freq = clade.frequencies[s]
				obs_freq = clade.frequencies[t]
				pred_freq = np.sum(np.exp(self.fitness(params, self.predictor_arrays[s][clade.tips[s],:])))/total_strain_freq
				clade_errors.append(np.absolute(pred_freq - obs_freq))
			seasonal_errors.append(np.mean(clade_errors))
		mean_error = np.mean(seasonal_errors)
		if any(np.isnan(seasonal_errors)+np.isinf(seasonal_errors)):
			mean_error = 1e10
		self.last_fit = mean_error
		if self.verbose>2: print params, self.last_fit
		return mean_error
		
	def fitness(self, params, pred):
		return np.sum(params*pred, axis=-1)


	def minimize_error(self):
		from scipy.optimize import fmin as minimizer
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
			self.params = 0.1+0.5*np.random.randn(len(self.predictors)) #0*np.ones(len(self.predictors))  # initial values
			self.minimize_error()
			params_stack.append((self.last_fit, self.params))

		self.params = params_stack[np.argmin([x[0] for x in params_stack])][1]
		self.model_fit(self.params)
		if self.verbose:
			print "best after",niter,"iterations\nfunction value:", self.last_fit
			print "fit parameters:"
			for pred, val in izip(self.predictors, self.params):
				print pred[0],':', val


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
		self.select_clades_for_fitting()
		self.learn_parameters(niter = niter)
		self.assign_fitness(self.seasons[-1])

def test(params):
	from io_util import read_json
	from tree_util import json_to_dendropy, to_Biopython, color_BioTree_by_attribute
	from Bio import Phylo
	tree_fname='data/tree_refine_10y_50v.json'
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

