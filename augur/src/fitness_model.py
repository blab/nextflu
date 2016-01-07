import argparse
import numpy as np
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
		
		final_date = calendar_date(self.time_interval[1])
		final_interval = (date.fromordinal(final_date.toordinal()-180), final_date)
		
		if self.estimate_coefficients:
			self.seasons.append(final_interval)
		else:
			self.seasons = [final_interval]

		self.predictors = []
		for p in predictors:
			self.setup_epitope_mask()
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
			
		self.global_std = 0*np.ones(len(self.predictors))
		if isinstance(predictor_input, dict):
			self.global_std = np.array([predictor_input[k][1] for k in predictors])

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


	def calc_predictors(self):
		for pred, func, kwargs in self.predictors:
			# calculate the predictors for all nodes of the tree and save as node.attr
			if pred!='dfreq':
				func(self.tree, attr = pred, **kwargs)

	def select_nodes_in_season(self, season):
		for node in self.tree.postorder_node_iter():
			if season in node.season_tips and len(node.season_tips[season])>0:
				node.alive=True
			else:
				node.alive=False

	def calc_time_censcored_tree_frequencies(self):
		print("fitting clade frequencies for seasons")
		region = "global_fit"
		freq_cutoff = 25.0
		total_pivots = 12
		pivots_fit = 2
		freq_window = 0.0
		from date_util import numerical_date
		for n in self.tree.preorder_node_iter():
			n.fit_frequencies = {}
			n.freq_slope = {}
		for s in self.seasons:
			time_interval = [numerical_date(s[0]) - freq_window, numerical_date(s[1])]
			pivots = np.linspace(time_interval[0], time_interval[1], total_pivots)
			n_nodes = len(self.tree.seed_node.season_tips[s])
			self.estimate_tree_frequencies(pivots=pivots, threshold = 40, regions=None,
								region_name = region, time_interval=time_interval)
			for n in self.tree.preorder_node_iter():
				if n.logit_freq[region] is not None:
					n.fit_frequencies[s] = np.minimum(freq_cutoff, np.maximum(-freq_cutoff,n.logit_freq[region]))
				else:
					n.fit_frequencies[s] = n.parent_node.fit_frequencies[s]
				try:
					slope, intercept, rval, pval, stderr = linregress(pivots[pivots_fit:], n.fit_frequencies[s][pivots_fit:])
					n.freq_slope[s] = slope
				except:
					import ipdb; ipdb.set_trace()
		# reset pivots in tree to global pivots
		self.tree.seed_node.pivots = self.pivots



	def calc_all_predictors(self, estimate_frequencies = True):
		if estimate_frequencies and 'dfreq' in [x[0] for x in self.predictors]:
			self.calc_time_censcored_tree_frequencies()
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
				if 'dfreq' in [x[0] for x in self.predictors]: node.dfreq = node.freq_slope[s]
				node.predictors[s] = np.array([node.__getattribute__(pred[0]) 
			                              for pred in self.predictors])
				if node.is_leaf():
					tmp_preds.append(node.predictors[s])
			self.predictor_arrays[s]=np.array(tmp_preds)

	def standardize_predictors(self):
		self.season_means = []
		self.season_std = []
		if self.verbose: print "standardize predictors for season"
		for s in self.seasons:
			self.season_means.append(self.predictor_arrays[s][self.tree.seed_node.season_tips[s],:].mean(axis=0))
			self.season_std.append(self.predictor_arrays[s][self.tree.seed_node.season_tips[s],:].std(axis=0))

		if self.estimate_coefficients:
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
		total_counts = {s:len(strain_list) for s, strain_list in self.tree.seed_node.season_tips.iteritems()}
		# prune seasons where few observations were made, only consecutive pairs
		# with sufficient tip count are retained
		self.fit_test_season_pairs = [(s,t) for s,t in izip(self.seasons[:-2], self.seasons[1:-1]) 
							if total_counts[s]>min_tips and total_counts[t]>min_tips]
		
		# for each pair of seasons with sufficient tip counts, 
		# select clades in a certain frequency window
		# keep track of these season / clade combinations in the dict clades_for_season
		self.clades_for_season = {}
		for s,t in self.fit_test_season_pairs:
			self.clades_for_season[(s,t)] = []
			for node in self.tree.postorder_node_iter():
				if  node.season_frequencies[s]>=min_freq and \
					node.season_frequencies[s]<max_freq and\
					node.season_frequencies[s]<node.parent_node.season_frequencies[s]:
					#node.aa_muts['HA1']!='':
					self.clades_for_season[(s,t)].append(node)


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


	def learn_parameters(self, niter = 10, fit_func = "af"):
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
		self.calc_tip_counts()
		self.calc_all_predictors(estimate_frequencies = estimate_frequencies)
		self.standardize_predictors()
		self.select_clades_for_fitting()
		self.prep_af()
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

