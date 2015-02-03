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

class fitness_model(object):

	def __init__(self, tree, predictors = None, verbose=0):
		'''
		parameters:
		tree -- tree of sequences for which a fitness model is to be determined
		'''
		self.tree = tree
		self.verbose=verbose
		if predictors is None:
			self.predictors = [
								('lb',calc_LBI, {'tau':0.0005, 'transform':lambda x:x}), 
								('ep',calc_epitope_distance,{}), 
								('ne',calc_nonepitope_distance,{})
								]
		else:
			self.predictors = predictors

		self.seasons = [ (date(year=y, month = 10, day = 1), date(year = y+1, month = 4, day=1)) 
						for y in xrange(ymin, ymax)]

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
		if self.verbose: print "selecting predictors for fitting"
		# short cut to total number of tips per seaons
		total_counts = {s:len(strain_list) for s, strain_list in self.tree.seed_node.seed_node.tips.iteritems()}
		# prune seasons where few observations were made, only consecutive pairs
		# with sufficient tip count are retained
		self.fit_test_season_pairs = [(s,t) for s,t in izip(self.seasons[:-1], self.seasons[1:]) 
							if total_counts[s]>min_tips and total_counts[t]>min_tips]

		# for each pair of seasons with sufficient tip counts, 
		# select clades in a certain frequency window. 
		self.freq_and_predictors = []
		for s,t in self.fit_test_season_pairs:
			tmp_freq = []
			tmp_pred = []
			for node in self.tree.postorder_node_iter():
				if node.frequencies[s]>=min_freq and node.frequencies[s]<max_freq:
					if node.frequencies[s]>1.01*np.max([c.frequencies[s] for c in node.child_nodes()]):
						tmp_freq.append([node.frequencies[s], node.frequencies[t]])
						tmp_pred.append(node.predictors[s])
			self.freq_and_predictors.append((np.array(tmp_freq), np.array(tmp_pred)))


	def model_fit_by_season(self, params):
		'''
		this function should work for a params = [1, ..., k] and a pred = n x k matrix 
		'''
		return [np.mean( np.absolute(freq[:,1] - freq[:,0]*np.exp(self.fitness(params, pred))) ) 		# mean absolute error of clade frequencies
				for freq, pred in self.freq_and_predictors]
#		return [np.mean((freq[:,1] - freq[:,0]*np.exp(self.fitness(params, pred)))**2) 					# mean squared errors of clade frequencies
#				for freq, pred in self.freq_and_predictors]
#		return [np.sum((np.log((freq[:,1]+pc)/(freq[:,0]+pc))-self.fitness(params, pred))**2)			# sum of squared errors of log fold changes
#				for freq, pred in self.freq_and_predictors]

	def model_fit(self, params):
		self.last_fit = np.mean(self.model_fit_by_season(params))
		if self.verbose>1: print params, self.last_fit
		return np.sum(self.last_fit)

	def prep_clades_for_fitting_tree(self):
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
					self.clades_for_season[(s,t)].append(node)
					
	def model_fit_tree(self, params):
		# walk through season pairs
		seasonal_errors = []
		for s,t in self.fit_test_season_pairs:		
			# normalize strain frequencies
			total_strain_freq = np.exp(self.fitness_tree(params, self.predictor_arrays[s][self.tree.seed_node.tips[s],:])).sum()
		
			# project clades forward according to strain makeup
			clade_errors = []
			test_clades = self.clades_for_season[(s,t)]
			for clade in test_clades:
				initial_freq = clade.frequencies[s]
				obs_freq = clade.frequencies[t]
				pred_freq = np.sum(np.exp(self.fitness_tree(params, self.predictor_arrays[s][clade.tips[s],:])))/total_strain_freq
				clade_errors.append(np.absolute(pred_freq - obs_freq))				
			seasonal_errors.append(np.mean(clade_errors))
		mean_error = np.mean(seasonal_errors)
		self.last_fit = mean_error
		if self.verbose>2: print params, self.last_fit
		return mean_error
		
	def fitness_tree(self, params, pred):
		return np.sum(params*pred, axis=-1)														
			
	def learn_parameters_tree(self):
		from scipy.optimize import fmin as minimizer
		self.params = 0*np.ones(len(self.predictors))  # initial values
		if self.verbose: 
			print "fitting parameters of the fitness model\ninitial function value:", self.model_fit_tree(self.params)
			
		self.params = minimizer(self.model_fit_tree, self.params, disp = self.verbose>1) # minimization
		if self.verbose:
			print "final function value:", self.last_fit
			print "fit parameters:"
			for pred, val in izip(self.predictors, self.params):
				print pred[0],':', val				

	def assign_fitness_tree(self, season):
		if self.verbose: print "calculating predictors for the last season"

		#FIXME: standardize predictors
		for node in self.tree.postorder_node_iter():
			if node.predictors[season] is not None:
				node.fitness = self.fitness_tree(self.params, node.predictors[season])
			else:
				node.fitness = 0.0

	def fitness(self, params, pred):
		return params[0]+np.sum(params[1:]*pred, axis=-1)

	def learn_parameters(self):
		from scipy.optimize import fmin
		# we seem to need an affine contribution since the normalization depends a lot on
		# whether only leafs or all internal clades are included in the mean. 
		self.params = 0*np.ones(1+len(self.predictors))  # initial values
		if self.verbose: 
			print "fitting parameters of the fitness model\ninitial function value:", self.model_fit(self.params)

		self.params = fmin(self.model_fit, self.params, disp = self.verbose>1) # minimization
		if self.verbose: 
			print "final function value:", self.last_fit
			print "fit parameters:"
			print "offset:", self.params[0]
			for pred, val in izip(self.predictors, self.params[1:]):
				print pred[0],':', val

	def assign_fitness(self, season):
		if self.verbose: print "calculating predictors for the last season"

		#FIXME: standardize predictors
		for node in self.tree.postorder_node_iter():
			if node.predictors[season] is not None:
				node.fitness = self.fitness(self.params, node.predictors[season])
			else:
				node.fitness = 0.0

	def predict(self):
		self.calc_tip_counts()
		self.calc_all_predictors()
		self.standardize_predictors()
		self.select_clades_for_fitting()
		self.learn_parameters()
		self.assign_fitness(self.seasons[-1])

	def predict_tree(self):
		self.calc_tip_counts()
		self.calc_all_predictors()
		self.standardize_predictors()
		self.prep_clades_for_fitting_tree()
		self.learn_parameters_tree()
		self.assign_fitness_tree(self.seasons[-1])

def test():
	from io_util import read_json
	from tree_util import json_to_dendropy, to_Biopython, color_BioTree_by_attribute
	from Bio import Phylo
	tree_fname='data/tree_refine_10y_50v.json'
	tree =  json_to_dendropy(read_json(tree_fname))
	fm = fitness_model(tree, verbose=3)
	fm.predict_tree()
	#btree = to_Biopython(tree)
	#color_BioTree_by_attribute(btree, 'fitness')
	#Phylo.draw(btree, label_func=lambda x:'')
	return fm

def main():
	from io_util import read_json
	from io_util import write_json	
	from tree_util import json_to_dendropy, dendropy_to_json

	tree_fname='data/tree_refine.json'
	tree =  json_to_dendropy(read_json(tree_fname))
	fm = fitness_model(tree, verbose=1)
	fm.predict_tree()
	out_fname = "data/tree_fitness.json"
	write_json(dendropy_to_json(tree.seed_node), out_fname)
	return out_fname

if __name__=="__main__":
	fm = test()
	#main()
