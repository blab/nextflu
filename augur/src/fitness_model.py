import argparse
import numpy as np
from scipy.interpolate import interp1d
import dendropy
from collections import defaultdict
from datetime import date
from date_util import calendar_date
from itertools import izip
from fitness_predictors import *

min_freq = 0.1
max_freq = 0.99
min_tips = 10
pc=1e-2
regularization = 1e-3
default_predictors = ['lb', 'ep', 'ne_star']

def dummy(tree, attr='dummy'):
	return None

class fitness_model(object):

	def __init__(self, predictor_input = ['ep', 'lb', 'dfreq'], verbose=0, **kwargs):
		'''
		parameters:
		tree -- tree of sequences for which a fitness model is to be determined
		predictor_input -- list of predictors to fit or dict of predictors to coefficients / std deviations
		'''
		self.verbose = verbose
		self.estimate_coefficients = True

		if isinstance(predictor_input, dict):
			predictors = predictor_input.keys()
			self.estimate_coefficients = False
		else:
			predictors = predictor_input
		if "estimate_fitness_model" in self.kwargs:
			if self.kwargs["estimate_fitness_model"]:
				self.estimate_coefficients = True

		self.seasons = [ (date(year=y, month = 10, day = 1), date(year = y+1, month = 4, day=1)) 
						for y in xrange(int(self.time_interval[0])+1, int(self.time_interval[1]))]
						
		self.initial_times = [2013.0, 2014.0]
		self.final_times = [2014.0, 2015.0]
		
		final_date = calendar_date(self.time_interval[1])
		final_interval = (date.fromordinal(final_date.toordinal()-180), final_date)
		
		if self.estimate_coefficients:
			self.seasons.append(final_interval)
		else:
			self.seasons = [final_interval]

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
			if p == 'dfreq':
				self.predictors.append(('dfreq',dummy,{}))
			if p == 'cHI':
				self.predictors.append(('cHI',dummy,{}))
		
		self.model_params = 0*np.ones(len(self.predictors))
		if isinstance(predictor_input, dict):
			self.model_params = np.array([predictor_input[k][0] for k in predictors])
			
		self.global_sds = 0*np.ones(len(self.predictors))
		if isinstance(predictor_input, dict):
			self.global_sds = np.array([predictor_input[k][1] for k in predictors])

	def prep_nodes(self):
		self.nodes = [node for node in self.tree.postorder_node_iter()]
		self.tips = [node for node in self.nodes if node.is_leaf()]
		self.rootnode = self.tree.seed_node

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
				for season, count in child.season_tips.iteritems():
					tmp_list[season].extend(count)

			node.season_tips = {}
			for season, strain_list in tmp_list.iteritems():
				node.season_tips[season] = np.array(strain_list, dtype=int)

			if node.is_leaf():
				node_date = date(*map(int, node.date.split('-')))
				tmp_season_list = [s for s in self.seasons if node_date>=s[0] and node_date<s[1]]
				if len(tmp_season_list)==1:
					node.season_tips[tmp_season_list[0]]=[leaf_count]
				leaf_count+=1
				node.tip_index = leaf_count

		# short cut to total number of tips per seaons
		total_counts = {s:len(strain_list) for s, strain_list in self.tree.seed_node.season_tips.iteritems()} 
		if self.verbose>1:
			for d,c in sorted(total_counts.items()): 
				print "number of tips in", d[0].strftime('%Y-%m-%d'), '--', \
											d[1].strftime('%Y-%m-%d'),':',c

		# calculate frequencies
		for node in self.tree.postorder_node_iter():
			node.season_frequencies = defaultdict(float)
			for season, strain_list in node.season_tips.iteritems():
				if season in total_counts:
					node.season_frequencies[season]=float(len(strain_list))/(total_counts[season]+1e-10)
				else:
					node.season_frequencies[season]=0

	def calc_node_frequencies(self):
		'''
		goes over all nodes and calculates frequencies at timepoints based on previously calculated frequency trajectories
		'''
		for node in self.nodes:		
			interpolation = interp1d(self.rootnode.pivots, node.freq['global'], kind='linear', bounds_error=False)
			node.initial_freqs = defaultdict(float)
			for time in self.initial_times:
				node.initial_freqs[time] = interpolation(time)
			node.final_freqs = defaultdict(float)				
			for time in self.final_times:
				node.final_freqs[time] = interpolation(time)
		# freq_arrays list *all* tips for each initial timepoint
		self.initial_freq_arrays={}
		for time in self.initial_times:		
			tmp_freqs = []			                              
			for tip in self.tips:
				tmp_freqs.append(tip.initial_freqs[time])
			self.initial_freq_arrays[time] = np.array(tmp_freqs)		
				

	def calc_predictors(self):
		for pred, func, kwargs in self.predictors:
			# calculate the predictors for all nodes of the tree and save as node.attr
			if pred!='dfreq':
				func(self.tree, attr = pred, **kwargs)

	def select_nodes_in_season(self, timepoint):
		# used by fitness_predictors:calc_LBI
		# TODO: fix me for continous time model
		cutoff = 0.01
		for node in self.nodes:
			#if season in node.season_tips and len(node.season_tips[season])>0:		
			if node.initial_freqs[timepoint] > cutoff:
				node.alive=True
			else:
				node.alive=False

	def calc_time_censored_tree_frequencies(self):
		print("fitting clade frequencies for timepoints")
		region = "global_fit"
		freq_cutoff = 25.0
		total_pivots = 12
		pivots_fit = 2
		freq_window = 1.0
		from date_util import numerical_date
		for node in self.nodes:
			node.fit_frequencies = {}
			node.freq_slope = {}
		for time in self.initial_times:
			time_interval = [time - freq_window, time]
			pivots = np.linspace(time_interval[0], time_interval[1], total_pivots)
			self.estimate_tree_frequencies(pivots=pivots, threshold = 40, regions=None,
								region_name = region, time_interval=time_interval)
			for node in self.nodes:
				if node.logit_freq[region] is not None:
					node.fit_frequencies[time] = np.minimum(freq_cutoff, np.maximum(-freq_cutoff,node.logit_freq[region]))
				else:
					node.fit_frequencies[time] = node.parent_node.fit_frequencies[time]
				try:
					slope, intercept, rval, pval, stderr = linregress(pivots[pivots_fit:], node.fit_frequencies[time][pivots_fit:])
					node.freq_slope[time] = slope
				except:
					import ipdb; ipdb.set_trace()
		# reset pivots in tree to global pivots
		self.rootnode.pivots = self.pivots



	def calc_all_predictors(self, estimate_frequencies = True):
		if estimate_frequencies and 'dfreq' in [x[0] for x in self.predictors]:
			self.calc_time_censored_tree_frequencies()
		# predictor_arrays list *all* tips for each timepoint
		self.predictor_arrays={}
		for node in self.nodes:
			node.predictors = {}
		for time in self.initial_times:
			if self.verbose: print "calculating predictors for time", time
			self.select_nodes_in_season(time)
			self.calc_predictors()
			for node in self.nodes:
				if 'dfreq' in [x[0] for x in self.predictors]: node.dfreq = node.freq_slope[time]
				node.predictors[time] = np.array([node.__getattribute__(pred[0]) 
			                              for pred in self.predictors])
			tmp_preds = []			                              
			for tip in self.tips:
				tmp_preds.append(tip.predictors[time])
			self.predictor_arrays[time]=np.array(tmp_preds)

	def standardize_predictors(self):
		self.predictor_means = {}
		self.predictor_sds = {}
		if self.verbose: print "standardizing predictors"
		for time in self.initial_times:
			values = self.predictor_arrays[time]
			weights = self.initial_freq_arrays[time]
			means = np.average(values, weights=weights, axis=0)
			variances = np.average((values-means)**2, weights=weights, axis=0)
			sds = np.sqrt(variances)
			self.predictor_means[time] = means
			self.predictor_sds[time] = sds

		if self.estimate_coefficients:
			self.global_sds = np.mean(self.predictor_sds.values(), axis=0)

		for time in self.initial_times:
			for node in self.nodes:
				if node.predictors[time] is not None:
					node.predictors[time] = (node.predictors[time]-self.predictor_means[time]) / self.global_sds
			self.predictor_arrays[time] -= self.predictor_means[time]
			self.predictor_arrays[time] /= self.global_sds		


	def select_clades_for_fitting(self):
		# for each initial time point, select clades that are within the specified frequency window
		# keep track in the dict initial_clades that maps timepoint to clade list
		self.initial_clades = {}
		for time in self.initial_times:
			self.initial_clades[time] = []
			for node in self.nodes:
				if  node.initial_freqs[time] >= min_freq and \
					node.initial_freqs[time] <= max_freq and \
					node.initial_freqs[time] < node.parent_node.initial_freqs[time]:
					self.initial_clades[time].append(node)


	def clade_fit(self, params):
		# walk through season pairs
		seasonal_errors = []
		self.pred_vs_true = []
		for s,t in self.fit_test_season_pairs:		
			# normalize strain frequencies
			total_strain_freq = np.exp(self.fitness(params, self.predictor_arrays[s][self.tree.seed_node.season_tips[s],:])).sum()
		
			# project clades forward according to strain makeup
			clade_errors = []
			tmp_pred_vs_true = []
			test_clades = self.clades_for_season[(s,t)]
			for clade in test_clades:
				initial_freq = clade.season_frequencies[s]
				obs_freq = clade.season_frequencies[t]
				pred_freq = np.sum(np.exp(self.fitness(params, self.predictor_arrays[s][clade.season_tips[s],:])))/total_strain_freq
				tmp_pred_vs_true.append((initial_freq, obs_freq, pred_freq))
				clade_errors.append(np.absolute(pred_freq - obs_freq))
			seasonal_errors.append(np.mean(clade_errors))
			self.pred_vs_true.append(np.array(tmp_pred_vs_true))
		mean_error = np.mean(seasonal_errors)
		if any(np.isnan(seasonal_errors)+np.isinf(seasonal_errors)):
			mean_error = 1e10
		self.last_fit = mean_error
		if self.verbose>2: print params, self.last_fit
		return mean_error + regularization*np.sum(params**2)

	def weighted_af(self, seqs, weights):
		af = np.zeros((4, seqs.shape[1]))
		for ni, nuc in enumerate('ACGT'):
			af[ni] += (weights*(seqs==nuc).T).sum(axis=1)/weights.sum()
		return af

	def af_fit(self, params):
		# walk through season pairs
		seasonal_errors = []
		self.pred_vs_true = []
		for s,t in self.fit_test_season_pairs:		
			weights = np.exp(self.fitness(params, self.predictor_arrays[s][self.tree.seed_node.season_tips[s],:]))
			pred_af = self.weighted_af(self.seqs[s],weights)
			#seasonal_errors.append(np.mean(np.sum((pred_af-self.af[t])**2, axis=0), axis=0))
			future_diameter = 0.5*np.sum(np.sum(self.af[t]*(1-self.af[t]), axis=0), axis=0)
			seasonal_errors.append(np.sum(np.sum(pred_af*(1-self.af[t]), axis=0), axis=0)-future_diameter)
			good_ind = self.af[s]*(1-self.af[s])>0.05
			self.pred_vs_true.append(np.array(zip(self.af[s][good_ind], self.af[t][good_ind], pred_af[good_ind])))

		mean_error = np.mean(seasonal_errors)
		if any(np.isnan(seasonal_errors)+np.isinf(seasonal_errors)):
			mean_error = 1e10
		self.last_fit = mean_error
		if self.verbose>2: print params, self.last_fit
		return mean_error + regularization*np.sum(params**2)
		
	def fitness(self, params, pred):
		return np.sum(params*pred, axis=-1)

	def minimize_clade_error(self):
		from scipy.optimize import fmin as minimizer
		if self.verbose:		
			print "initial function value:", self.clade_fit(self.model_params)
			print "initial parameters:", self.model_params
		self.model_params = minimizer(self.clade_fit, self.model_params, disp = self.verbose>1)
		if self.verbose:
			print "final function value:", self.clade_fit(self.model_params)		
			print "final parameters:", self.model_params, '\n'		

	def prep_af(self):
		# TODO: fix me for continous time model
		if not hasattr(self,'variable_nuc'):
			self.determine_variable_positions()
		self.seqs = {}
		self.af = {}
		fit_aln = np.zeros((len(self.viruses), len(self.variable_nuc)), dtype='S1')
		for leaf in self.tree.leaf_iter():
			fit_aln[leaf.tip_index] = np.fromstring(leaf.seq, 'S1')[self.variable_nuc]
		for s in self.seasons:
			self.seqs[s] = fit_aln[self.tree.seed_node.season_tips[s]]
			self.af[s] = self.weighted_af(self.seqs[s], np.ones(len(self.seqs[s])))

	def minimize_af_error(self):
		from scipy.optimize import fmin as minimizer
		if self.verbose:		
			print "initial function value:", self.af_fit(self.model_params)
			print "initial parameters:", self.model_params
		self.model_params = minimizer(self.af_fit, self.model_params, disp = self.verbose>1)
		if self.verbose:
			print "final function value:", self.af_fit(self.model_params)		
			print "final parameters:", self.model_params, '\n'		


	def learn_parameters(self, niter = 10, fit_func = "clade"):
		if fit_func=='clade':
			minimize_error=self.minimize_clade_error
			fit_func=self.clade_fit
		elif fit_func=="af":
			minimize_error=self.minimize_af_error
			fit_func=self.af_fit
		else:
			print("fit function", fit_func,"does not exist")
			raise NotImplementedError

		print "fitting parameters of the fitness model\n"

		params_stack = []

		if self.verbose:
			print "null parameters"
		self.model_params = 0*np.ones(len(self.predictors))  # initial values
		minimize_error()
		params_stack.append((self.last_fit, self.model_params))
		
		for ii in xrange(niter):
			if self.verbose:
				print "iteration:", ii+1
			self.model_params = 0.5*np.random.randn(len(self.predictors)) #0*np.ones(len(self.predictors))  # initial values
			minimize_error()
			params_stack.append((self.last_fit, self.model_params))

		self.model_params = params_stack[np.argmin([x[0] for x in params_stack])][1]
		fit_func(self.model_params)
		if self.verbose:
			print "best after",niter,"iterations\nfunction value:", self.last_fit
			print "fit parameters:"
			for pred, val in izip(self.predictors, self.model_params):
				print pred[0],':', val


	def assign_fitness(self, season):
		if self.verbose: print "calculating predictors for the last season"

		#FIXME: standardize predictors
		for node in self.tree.postorder_node_iter():
			if node.predictors[season] is not None:
				node.fitness = self.fitness(self.model_params, node.predictors[season])
			else:
				node.fitness = 0.0
				
		weights = np.exp(self.fitness(self.model_params, self.predictor_arrays[season][self.tree.seed_node.season_tips[season],:]))
		pred_af = self.weighted_af(self.seqs[season], weights)

		for node in self.tree.postorder_node_iter():
			if node.predictors[season] is not None:
				seq = np.fromstring(node.seq, 'S1')[self.variable_nuc]
				seq_indicators = self.weighted_af(np.array([seq]), np.ones(1))		
				node.pred_distance = np.sum(np.sum(pred_af*(1-seq_indicators), axis=0), axis=0)
			else:
				node.pred_distance = 0.0

	def predict(self, niter = 10, estimate_frequencies = True):
		self.prep_nodes()
		self.calc_tip_counts()
		self.calc_node_frequencies()
		self.calc_all_predictors(estimate_frequencies = estimate_frequencies)
		self.standardize_predictors()
		self.select_clades_for_fitting()
#		self.prep_af()
		raise ValueError('End of redo')		
		if self.estimate_coefficients:
			self.learn_parameters(niter = niter, fit_func = "clade")
		self.assign_fitness(self.seasons[-1])

	def validate_prediction(self):
		import matplotlib.pyplot as plt
		from scipy.stats import spearmanr 
		fig, axs = plt.subplots(1,4, figsize=(10,5))
		for season, pred_vs_true in izip(self.fit_test_season_pairs, self.pred_vs_true):
			# 0: initial, 1: observed, 2: predicted
			axs[0].scatter(pred_vs_true[:,1], pred_vs_true[:,2])
			axs[1].scatter(pred_vs_true[:,1]/pred_vs_true[:,0], 
						   pred_vs_true[:,2]/pred_vs_true[:,0], c=pred_vs_true[0])
			for s, o, p  in pred_vs_true:
				axs[2].arrow(s,s, o-s,p-s)
			axs[3].scatter(pred_vs_true[:,0], 
						   (pred_vs_true[:,2]+0.01)/(pred_vs_true[:,1]+0.01))				

		# pred_vs_true is initial, observed, predicted
		tmp = np.vstack(self.pred_vs_true)
		print("Spearman's rho, null",spearmanr(tmp[:,0], tmp[:,1]))
		print("Spearman's rho, raw",spearmanr(tmp[:,1], tmp[:,2]))
		print("Spearman's rho, rel",spearmanr(tmp[:,1]/tmp[:,0], 
											  tmp[:,2]/tmp[:,0]))
	
		growth_list = [pred > initial for (initial, obs, pred) in tmp if obs > initial]
		print ("Correct at predicting growth", growth_list.count(True) / float(len(growth_list)))

		decline_list = [pred < initial for (initial, obs, pred) in tmp if obs < initial]
		print ("Correct at predicting decline", decline_list.count(True) / float(len(decline_list)))

		axs[0].set_ylabel('predicted')
		axs[0].set_xlabel('observed')
		axs[1].set_ylabel('predicted/initial')
		axs[1].set_xlabel('observed/initial')
		axs[1].set_yscale('linear')
		axs[1].set_xscale('linear')
		axs[2].set_ylabel('predicted')
		axs[2].set_xlabel('observed')
		axs[2].set_ylim(-0.1, 1.1)
		axs[2].set_xlim(-0.1, 1.1)
		axs[3].set_ylabel('predicted / observed')
		axs[3].set_xlabel('initial')
		axs[3].set_yscale('log')		
		plt.tight_layout()

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

