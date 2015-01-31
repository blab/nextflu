import numpy as np
import dendropy
from collections import defaultdict
from datetime import date
from itertools import izip
from fitness_predictors import *

ymin = 2005
ymax = 2015
min_freq = 0.1
max_freq = 0.5
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
			self.predictors = [('lb',calc_LBI, {'tau':0.0005, 'transform':lambda x:x**1.0}), 
								('ep',calc_epitope_distance,{}), 
								('ne',calc_nonepitope_distance,{})]
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
		for node in self.tree.postorder_node_iter():
			node.tip_counts = defaultdict(int)
			for child in node.child_nodes():
				for season, count in child.tip_counts.iteritems():
					node.tip_counts[season]+=count
			if node.is_leaf():
				node_date = date(*map(int, node.date.split('-')))
				tmp_season_list = [s for s in self.seasons if node_date>=s[0] and node_date<s[1]]
				if len(tmp_season_list)==1:
					node.tip_counts[tmp_season_list[0]]=1

		# short cut to total number of tips per seaons
		total_counts = self.tree.seed_node.tip_counts
		if self.verbose>1:
			for d,c in sorted(total_counts.items()): 
				print "number of tips in", d[0].strftime('%Y-%m-%d'), '--', \
											d[1].strftime('%Y-%m-%d'),':',c
		# prune seasons where few observations were made, only consecutive pairs
		# with sufficient tip count are retained
		self.fit_test_season_pairs = [(s,t) for s,t in izip(self.seasons[:-1], self.seasons[1:]) 
							if total_counts[s]>min_tips and total_counts[t]>min_tips]

		# calculate frequencies
		for node in self.tree.postorder_node_iter():
			node.frequencies = defaultdict(float)
			for season, count in node.tip_counts.iteritems():
				node.frequencies[season]=float(count)/(total_counts[season]+1e-10)


	def calc_predictors(self):
		for pred, func, kwargs in self.predictors:
			# calculate the predictors for all nodes of the tree and save as node.attr
			func(self.tree, attr = pred, **kwargs)

	def make_flat_lists_per_season(self):
		# for each pair of seaons with sufficient tip counts, 
		# make a list of clades that are to be used in model fitting
		self.clades_in_season = []  # list of clades for each season...
		for s,t in self.fit_test_season_pairs:
			tmp_clades = []
			for node in self.tree.postorder_node_iter():
				if node.frequencies[s]>=min_freq and node.frequencies[s]<max_freq:
					if node.frequencies[s]>1.01*np.max([c.frequencies[s] for c in node.child_nodes()]):
						tmp_clades.append(node)
			self.clades_in_season.append(tmp_clades)

		# calculate the predictors for the clades used for model fitting
		self.freq_and_predictors = []
		for (s,t), clades in izip(self.fit_test_season_pairs, self.clades_in_season):
			if self.verbose>1:
				print "calculating predictors for ", s[0].strftime('%Y-%m-%d'),'--',s[1].strftime('%Y-%m-%d')
			for node in self.tree.postorder_node_iter():
				if node.tip_counts[s]>0:
					node.alive=True
				else:
					node.alive=False
			self.calc_predictors()
			tmp_freq = [[x.frequencies[s], x.frequencies[t]] for x in clades]
			#tmp_pred = [[0 for pred in self.predictors] for x in clades]
			tmp_pred = [[x.__getattribute__(pred[0]) for pred in self.predictors] for x in clades]
			self.freq_and_predictors.append([np.array(tmp_freq), np.array(tmp_pred)])


	def standardize_predictors(self):
		self.season_means = []
		self.season_std = []
		for freq, pred in self.freq_and_predictors:
			self.season_means.append(pred.mean(axis=0))
			self.season_std.append(pred.std(axis=0))
			pred-=self.season_means[-1]
		self.time_averaged_std = np.array(self.season_std).mean(axis=0)
		for freq_pred in self.freq_and_predictors:
			freq_pred[1]/=self.time_averaged_std



	def model_fit_by_season(self, params):
		'''
		this function should work for a params = [1, ..., k] and a pred = n x k matrix 
		'''
#		return [np.sum((freq[:,1] - freq[:,0]*np.exp(self.fitness(params, pred)))**2) 
#				for freq, pred in self.freq_and_predictors]
		return [np.sum((np.log((freq[:,1]+pc)/(freq[:,0]+pc))-self.fitness(params, pred))**2) 
				for freq, pred in self.freq_and_predictors]

	def model_fit(self, params):
		self.last_fit = np.sum(self.model_fit_by_season(params))
		return np.sum(self.last_fit)

	def fitness(self, params, pred):
		return np.sum(params*pred, axis=-1)

	def learn_parameters(self):
		from scipy.optimize import fmin
		if self.verbose: print "fitting parameters of the fitness model"
		self.params = 0*np.ones(len(self.predictors))  # initial values
		self.params = fmin(self.model_fit, self.params, disp = self.verbose>0) # minimzation. need to see what kind of minimization is useful 
		if self.verbose: 
			print "fit parameters:"
			for pred, val in izip(self.predictors, self.params):
				print pred[0],':', val

	def assign_fitness(self):
		s =  self.seasons[-1]
		if self.verbose: print "calculating predictors for the last season"
		for node in self.tree.postorder_node_iter():
			if node.tip_counts[s]>0:
				node.alive=True
			else:
				node.alive=False
		self.calc_predictors()
		#FIXME: standardize predictors
		for node in self.tree.postorder_node_iter():
			node.fitness = self.fitness(self.params, np.array([node.__getattribute__(pred[0]) 
																for pred in self.predictors]))


	def predict(self):
		self.calc_tip_counts()
		self.make_flat_lists_per_season()
		self.standardize_predictors()
		self.learn_parameters()
		self.assign_fitness()



if __name__=="__main__":
	from io_util import *
	from tree_util import *
	from Bio import Phylo
	tree_fname='data/tree_refine.json'
	tree =  json_to_dendropy(read_json(tree_fname))
	fm = fitness_model(tree, verbose=2)
	fm.predict()


	btree = to_Biopython(tree)
	color_BioTree_by_attribute(btree, 'fitness')
	Phylo.draw(btree, label_func=lambda x:'')

