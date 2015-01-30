import numpy as np
import dendropy
from collections import defaultdict

ymin = 2000
ymax = 2015

class fitness_model(object):

	def __init__(self, tree, predictors = None):
		'''
		parameters:
		tree -- tree of sequences for which a fitness model is to be determined
		'''
		self.tree = tree
		if predictors is None:
			self.predictors = [('lb',calc_LBI), ('ep',calc_epitope_distance), ('ne',calc_nonepitope_distance)]
		else:
			self.predictors = predictors

		from datetime import date
		self.seasons = [ (date(year=y, month = 10, day = 1), date(year = y+1, month = 4, day=1)) 
						for y in xrange(ymin, ymax)]

	def calc_tip_counts(self):
		for node in self.T.postorder_node_iter():
			node.tip_counts = defaultdict(int)
			for child in node.child_nodes():
				for season, count in child.tip_counts.iteritems():
					node.tip_counts[season]+=count
			if node.is_leaf():
				node_date = date(*map(int, node.date.split('-')))
				tmp_season_list = [s for s in self.seasons if node_date>=s[0] and node_date<s[1]]
				if len(tmp_season_list)==1:
					node.tip_counts[tmp_season_list[0]]=1

		total_counts = self.T.seed_node.tip_counts
		self.test_fit_season_pairs = [(s,n) for s,n in izip(self.seasons[:-1], self.seasons[1:] 
							if total_counts[s]>min_tips and total_counts[n]>min_tips]

		for node in self.T.postorder_node_iter():
			node.frequencies = defaultdict(float)
			for season, count in child.tip_counts.iteritems():
				node.frequencies[season]=float(count)/total_counts[season]


	def calc_predictors(self):
		for pred, func in self.predictors:
			# calculate the predictors for all nodes of the tree and save as node.attr
			func(self.tree.seed_node, attr = pred)

	def make_flat_lists_per_season(self):
		self.clades_in_season = []  # list of clades for each season...
		for s,t in self.test_fit_season_pairs:
			tmp_clades = []
			for node in self.T.postorder_node_iter():
				if node.frequencies[s]>=cut_off_frequency:
					tmp_clades.append(node)

		self.freq_and_predictors = []
		for (s,t), clades in izip(self.test_fit_season_pairs, self.clades_in_season):
			tmp_freq = [[x.frequencies[s], x.frequencies[t]] for x in clades]
			tmp_pred = [[x.__getattr__[pred[0]] for pred in self.predictors] for x in clades]
			self.freq_and_predictors.append((np.array(tmp_freq), np.array(tmp_pred))


	def standardize_predictors(self):
		self.season_means = []
		self.season_variances = []
		for freq, pred in self.freq_and_predictors:
			self.season_means.append(pred(axis=0))
			self.season_variances.append(pred(axis=0))
			pred-=self.season_means[-1]
			## TODO VARIANCE STANDARDIZATION

	def model_fit_by_season(self, params):
		'''
		this function should work for a params = [1, ..., k] and a pred = n x k matrix 
		'''
		return [ np.sum((np.log(freq[1]/freq[0]) - self.fitness(params, pred))**2) for freq, pred in self.seasons]

	def model_fit(self, params):
		self.last_fit = np.sum(self.model_fit_by_season(params))
		return np.sum(self.last_fit)

	def fitness(self, params, pred):
		return np.sum(params*pred, axis=1)

	def learn_parameters(self):
		from scipy.optimize import fmin
		self.params = np.array([1,1,1])  # initial values
		self.params = fmin(self.model_fit, self.params) # minimzation. need to see what kind of minimization is useful 

