# estimates clade frequencies using SMC

from scipy.interpolate import interp1d
import time
from io_util import *
from tree_util import *
from seq_util import *
from date_util import *
import numpy as np

pc=1e-4
dfreq_pc = 1e-2

clade_designations = { "3C3.a":[(128,'A'),(142,'G'), (159,'S')],
					   "3C3":[(128,'A'),(142,'G'), (159,'F')],
					   "3C2.a":[(144,'S'), (159,'Y'), (225,'D'), (311,'H'),(489,'N')],
					   "3C2":[(144,'S'), (159,'F'), (225,'D'), (311,'H'),(489,'N')],
						}

def running_average(obs, ws):
	'''
	calculates a running average
	obs 	--	observations
	ws 		--	winodw size (number of points to average)
	'''
	tmp_vals = np.convolve(np.ones(ws, dtype=float)/ws, obs, mode='same')
	# fix the edges. using mode='same' assumes zeros outside the range
	tmp_vals[:ws//2]*=float(ws)/np.arange(ws//2,ws)
	tmp_vals[-ws//2+1:]*=float(ws)/np.arange(ws-1,ws//2,-1.0)
	return tmp_vals

def fix_freq(freq, pc):
	'''
	restricts frequencies to the interval [pc, 1-pc]
	'''
	return np.minimum(1-pc, np.maximum(pc,freq))

def get_pivots(start, stop):
	return np.arange(2012, 2015.3, 1.0/12)

def logit_transform(freq):
	return np.log(freq/(1-freq))

def logit_inv(logit_freq):
	tmp_freq = np.exp(logit_freq)
	return tmp_freq/(1.0+tmp_freq)

def pq(p):
	return p*(1-p)

class frequency_estimator(object):
	'''
	estimates a smooth frequency trajectory given a series of time stamped
	0/1 observations. The most likely set of frequencies at specified pivot values
	is deterimned by numerical minimization. Likelihood consist of a bernoulli sampling
	term as well as a term penalizing rapid frequency shifts. this term is motivated by 
	genetic drift, i.e., sampling variation.
	'''

	def __init__(self, observations, pivots = None, stiffness = 20.0, logit=False, verbose = 0):
		self.tps = np.array([x[0] for x in observations])
		self.obs = np.array([x[1]>0 for x in observations])
		self.stiffness = stiffness
		self.interolation_type = 'cubic'
		self.logit = logit
		self.verbose=verbose
		# make sure they are searchsorted
		tmp = np.argsort(self.tps)
		self.tps = self.tps[tmp]
		self.obs = self.obs[tmp]

		if pivots is None:
			self.pivots = get_pivots(self.tps[0], self.tps[1])
		elif np.isscalar(pivots):
			self.pivot_tps = np.linspace(self.tps[0], self.tps[-1], pivots)
		else:
			self.pivots = pivots

		# generate a useful initital case from a running average of the counts
		ws=100
		tmp_vals = running_average(self.obs, ws)
		tmp_interpolator = interp1d(self.tps, tmp_vals)
		if self.logit:
			freq = tmp_interpolator(self.pivot_tps)
			self.pivot_freq = logit_transform(fix_freq(freq, pc))
		else:
			self.pivot_freq = tmp_interpolator(self.pivot_tps)
		if self.verbose:
			print "Initial pivots:", self.pivot_freq


	def stiffLH(self, pivots):
		if self.logit: # if logit, convert to frequencies
			freq = logit_inv(pivots)
		else:
			freq = pivots
		dfreq = np.diff(freq)
		# return wright fisher diffusion likelihood for frequency change. 
		return -0.25*self.stiffness*np.sum(dfreq**2/np.diff(self.pivot_tps)/pq(fix_freq(freq[:-1],dfreq_pc)))


	def logLH(self, pivots):
		freq = interp1d(self.pivot_tps, pivots, kind=self.interolation_type)
		if self.logit: # if logit, convert to frequencies
			estfreq = fix_freq(logit_inv(freq(self.tps)), pc)
		else:
			estfreq = fix_freq(freq(self.tps), pc)
		stiffness_LH = self.stiffLH(pivots)
		bernoulli_LH = np.sum(np.log(estfreq[self.obs])) + np.sum(np.log((1-estfreq[~self.obs])))
		LH = stiffness_LH + bernoulli_LH 
		if self.verbose>2: print "LH:",bernoulli_LH,stiffness_LH
		if self.logit:
			return -LH/len(self.obs)+0.0001*np.mean(pivots**2) # penalize too large or two small pivots
		else:
			return -LH/len(self.obs)+100000*(np.sum((pivots<0)*np.abs(pivots))+np.sum((pivots>1)*np.abs(pivots-1)))


	def learn(self):
		from scipy.optimize import fmin_powell as minimizer
		self.pivot_freq = minimizer(self.logLH, self.pivot_freq)
		self.frequency_estimate = interp1d(self.pivot_tps, self.pivot_freq, kind=self.interolation_type, bounds_error=False)


def estimate_clade_frequency(node, all_dates, tip_to_date_index):
	dt = 0.5 # time interval to include before the first after the last sample from a clade
	first_date = all_dates[tip_to_date_index[node.tips[0]]]
	last_date = all_dates[tip_to_date_index[node.tips[-1]]]

	# add dt on both ends
	start_index = max(0,np.searchsorted(all_dates, first_date-dt))
	stop_index = min(np.searchsorted(all_dates, last_date+dt), all_dates.shape[0]-1)

	# extract time points and the subset of observations that fall in the clade.
	tps = all_dates[start_index:stop_index]
	obs = np.in1d(tps, all_dates[tip_to_date_index[node.tips]])
	last_negative = np.argmin(obs - np.linspace(0,0.5, obs.shape[0]))
	stop_date = tps[last_negative]+dt
	new_stop_index = min(np.searchsorted(all_dates, stop_date), all_dates.shape[0]-1)
	if stop_index>new_stop_index:
		tps = tps[:new_stop_index-stop_index] 
		obs = obs[:new_stop_index-stop_index] 

	# make n pivots a year
	pivots = get_pivots(tps[0], tps[1])
	fe = frequency_estimator(zip(tps, obs), pivots=pivots, stiffness=2.0, logit=True)
	fe.learn()
	# return the final interpolation object
	return fe.frequency_estimate

def estimate_tree_frequencies(tree):
	'''
	loop over nodes of the tree and estimate frequencies of all clade above a certain size
	'''
	all_dates = []
	# loop over all nodes, make time ordered lists of tips
	for node in tree.postorder_node_iter():
		tmp_tips = []
		if node.is_leaf():
			all_dates.append(node.num_date)
			tmp_tips.append((node.tip_index, node.num_date))			
		for child in node.child_nodes():
			tmp_tips.extend(child.tips)
		node.tips = np.array([x for x in sorted(tmp_tips, key = lambda x:x[1] )])
	# erase the dates from the tip lists and cast to int such that they can be used for indexing
	for node in tree.postorder_node_iter():
		node.tips = np.array(node.tips[:,0], dtype=int)

	# sort the dates and provide a reverse ordering as a mapping of tip indices to dates
	all_dates = np.array(all_dates)
	leaf_order = np.argsort(all_dates)
	reverse_order = np.argsort(leaf_order)
	all_dates = all_dates[leaf_order]

	for node in tree.postorder_node_iter():
		if node.tips.shape[0]>50: # only do this for reasonably large clades
			print "# leafs:", node.tips.shape
			print "date:", node.date
			node.freq_est = estimate_clade_frequency(node, all_dates, reverse_order)
		else:
			node.freq_est=None


def estimate_genotype_frequency(tree, gt, time_interval=None):
	'''
	estimate the frequency of a particular genotype specified 
	gt   --		[(position, amino acid), ....]
	'''
	all_dates = []
	observations = []
	for node in tree.leaf_iter():
		is_gt = all([node.aa_seq[pos]==aa for pos, aa in gt])
		if time_interval is not None:
			good_time = (node.num_date>= time_interval[0]) and (node.num_date<time_interval[1])
		else:
			good_time=True
		if good_time:
			all_dates.append(node.num_date)
			observations.append(is_gt)
		leaf_count+=1

	all_dates = np.array(all_dates)
	leaf_order = np.argsort(all_dates)
	tps = all_dates[leaf_order]
	obs = np.array(observations)[leaf_order]

	# define pivots and estimate
	pivots = get_pivots(tps[0], tps[1])
	print "number of pivots:", len(pivots)
	fe = frequency_estimator(zip(tps, obs), npivots=npivots, stiffness=5.0, logit=True, verbose = 0)
	fe.learn()
	return fe.frequency_estimate, (tps,obs)

def determine_major_genotypes(tree, HA1=True, positions=None, time_interval=None):
	from collections import defaultdict
	gt_counts = defaultdict(int)
	leaf_count = 0
	for node in tree.leaf_iter():
		if time_interval is not None:
			good_time = (node.num_date>= time_interval[0]) and (node.num_date<time_interval[1])
		else:
			good_time=True
		if good_time:
			if positions is not None:
					gt_counts["".join([node.aa_seq[pos] for pos in positions])]+=1
			else:
				if HA1:
					gt_counts[get_HA1(node.aa_seq)]+=1
				else:
					gt_counts[node.aa_seq]+=1
		leaf_count+=1

	return gt_counts

def genotype_to_mutations(ref_aa_seq, gt):
	from itertools import izip
	return ",".join([a+str(pos+1)+b for pos, (a,b) in enumerate(izip(ref_aa_seq, gt)) if a!=b])

def mutation_counts(gt_counts):
	from collections import defaultdict
	mut_counts = defaultdict(int)
	for gt, val in gt_counts.iteritems():
		for mut in gt.split(','):
			mut_counts[mut] += val
	return mut_counts



def test():
	import matplotlib.pyplot as plt
	tps = np.sort(100 * np.random.uniform(size=100))
	freq = [0.1]
	logit = True
	stiffness=100
	s=-0.02
	for dt in np.diff(tps):
		freq.append(freq[-1]*np.exp(-s*dt)+np.sqrt(2*np.max(0,freq[-1]*(1-freq[-1]))*dt/stiffness)*np.random.normal())
	obs = np.random.uniform(size=tps.shape)<freq
	fe = frequency_estimator(zip(tps, obs), pivots=10, stiffness=stiffness, logit=logit)
	fe.learn()
	plt.figure()
	plt.plot(tps, freq, 'o', label = 'actual frequency')
	freq = fe.frequency_estimate(fe.tps)
	if logit: freq = logit_inv(freq)
	plt.plot(fe.tps, freq, '-', label='interpolation')
	plt.plot(tps, (2*obs-1)*0.05, 'o')
	plt.plot(tps[obs], 0.05*np.ones(np.sum(obs)), 'o', c='r', label = 'observations')
	plt.plot(tps[~obs], -0.05*np.ones(np.sum(1-obs)), 'o', c='g')
	plt.plot(tps, np.zeros_like(tps), 'k')
	ws=20
	r_avg = running_average(obs, ws)
	plt.plot(fe.tps[ws/2:-ws/2+1], np.convolve(np.ones(ws, dtype=float)/ws, obs, mode='valid'), 'r', label = 'running avg')
	plt.plot(fe.tps, r_avg, 'k', label = 'running avg')
	plt.legend(loc=2)


def main():
	# load tree
	import matplotlib.pyplot as plt
	from io_util import read_json
	from tree_util import json_to_dendropy, to_Biopython, color_BioTree_by_attribute
	from Bio import Phylo
	time_interval=(2012, 2015.5)

	tree_fname='data/tree_refine_10y_50v.json'
	koel_sites =  [160, 170, 171, 173, 174, 204, 208]
	tree =  json_to_dendropy(read_json(tree_fname))
	gt_counts = determine_major_genotypes(tree, time_interval=time_interval, positions =koel_sites, HA1=False)
	tree.seed_node.aa_seq = translate(tree.seed_node.seq)

	for clade_name, clade_gt in clade_designations.iteritems():
		print clade_name, clade_gt
		if True:
			freq, (tps, obs) = estimate_genotype_frequency(tree, [(pos+15, aa) for pos, aa in clade_gt], time_interval)
			#
			if np.mean(obs)<0.99:
				grid_tps = np.linspace(time_interval[0], time_interval[1], 100)
				freq_est = np.exp(freq(grid_tps))
				freq_est = freq_est/(1+freq_est)
				plt.plot(grid_tps, freq_est, label=clade_name, lw=2)
				#r_avg = running_average(obs, 100)
				#plt.plot(tps, r_avg)
	plt.legend()
	ticloc = np.arange(time_interval[0], int(time_interval[1])+1,1)
	plt.xticks(ticloc, map(str, ticloc))
	plt.xlim([time_interval[0], time_interval[1]+1])
#
#	gt_as_muts_counts = {genotype_to_mutations(get_HA1(tree.seed_node.aa_seq), gt):val for gt,val in gt_counts.iteritems()}
#	mut_counts = mutation_counts(gt_as_muts_counts)
#	for mut, val in mut_counts.iteritems():
#		if val>50:
#			freq, (tps, obs) = estimate_genotype_frequency(tree, [(int(mut[1:-1])+15,mut[-1])], time_interval)
#			#
#			if np.mean(obs)<0.99:
#				grid_tps = np.linspace(time_interval[0], time_interval[1], 100)
#				freq_est = np.exp(freq(grid_tps))
#				freq_est = freq_est/(1+freq_est)
#				plt.plot(grid_tps, freq_est, label=mut)
#				r_avg = running_average(obs, 100)
#				plt.plot(tps, r_avg)
#	# return the final interpolation object

#	estimate_tree_frequencies(tree)
#	for node in tree.postorder_node_iter():
#		if node.freq_est is not None:
#			tps = np.linspace(node.freq_est.x[0], node.freq_est.x[-1],100)
#			logit_freq = np.exp(node.freq_est(tps))
#			plt.plot(tps, logit_freq/(1+logit_freq))
#	return tree, gt_counts

if __name__=="__main__":
	test()
	#tree,gt_counts = main()


