import numpy as np
import dendropy


class fitness_model(object):

	def __init__(self, tree, predictors = None):
		'''
		parameters:
		tree -- tree of sequences for which a fitness model is to be determined
		'''
		self.tree = tree
		if predictors is None:
			self.predictors = {'lb':calc_LBI, 'ep':calc_epitope_distance, 'ne':calc_nonepitope_distance}
		else:
			self.predictors = predictors

	def calc_predictors(self):
		for pred, func in self.predictors.iteritems():
			# calculate the predictors for all nodes of the tree and save as node.attr
			func(selt.tree.seed_node, attr = pred)

	def make_flat_lists_per_season(self):
		self.clades_in_season = []  # list of clades for each season...
		
		self.freq_and_predictors = []
		for clades in self.seasons:
			tmp_freq = [[x.freq, x.next_freq] for x in clades]
			tmp_pred = [[x.__getattr__[pred] for pred in self.predictors] for x in clades]
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

