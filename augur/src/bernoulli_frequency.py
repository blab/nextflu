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
time_interval = (2012.0, 2015.0)
flu_stiffness = 10.0
pivots_per_year = 6.0

clade_designations = { "3C3.a":[(128,'A'),(142,'G'), (159,'S')],
					   "3C3":[(128,'A'),(142,'G'), (159,'F')],
					   "3C2.a":[(144,'S'), (159,'Y'), (225,'D'), (311,'H'),(489,'N')],
					   "3C2":[(144,'S'), (159,'F'), (225,'D'), (311,'H'),(489,'N')],
						}

cols  = np.array([(166,206,227),(31,120,180),(178,223,138),(51,160,44),(251,154,153),(227,26,28),(253,191,111),(255,127,0),(202,178,214),(106,61,154)], dtype=float)/255
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
	freq[np.isnan(freq)]=pc
	return np.minimum(1-pc, np.maximum(pc,freq))

def get_pivots(start=None, stop=None):
	return np.arange(np.floor(time_interval[0]*pivots_per_year), np.ceil(time_interval[1]*pivots_per_year)+0.5, 1.0)/pivots_per_year

def get_extrapolation_pivots(start=None, dt=0.5):
	return np.arange(np.floor(time_interval[1]*pivots_per_year), np.ceil((dt+time_interval[1])*pivots_per_year)+0.5, 1.0)/pivots_per_year


def logit_transform(freq):
	return np.log(freq/(1-freq))

def logit_inv(logit_freq):
	tmp_freq = np.exp(logit_freq)
	return tmp_freq/(1.0+tmp_freq)

def pq(p):
	return p*(1-p)

def extrapolation(freq_interp,x):
	def ep(freq_interp, x):
		if x>freq_interp.x[-1]:
			return freq_interp.y[-1] + (freq_interp.y[-1] - freq_interp.y[-2])/(freq_interp.x[-1] - freq_interp.x[-2]) * (x-freq_interp.x[-1])
		elif x<freq_interp.x[0]:
			return freq_interp.y[0] + (freq_interp.y[1] - freq_interp.y[0])/(freq_interp.x[1] - freq_interp.x[0]) * (x-freq_interp.x[0])
		else:
			return float(freq_interp(x))

	if np.isscalar(x): 
		return ep(freq_interp,x)
	else:
		return [ep(freq_interp,tmp_x) for tmp_x in x]	

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
		self.interolation_type = 'linear'
		self.logit = logit
		self.verbose=verbose
		# make sure they are searchsorted
		tmp = np.argsort(self.tps)
		self.tps = self.tps[tmp]
		self.obs = self.obs[tmp]

		if pivots is None:
			self.pivot_tps = get_pivots(self.tps[0], self.tps[1])
		elif np.isscalar(pivots):
			self.pivot_tps = np.linspace(self.tps[0], self.tps[-1], pivots)
		else:
			self.pivot_tps = pivots

		# generate a useful initital case from a running average of the counts
		ws=40
		tmp_vals = running_average(self.obs, ws)
		tmp_interpolator = interp1d(self.tps, tmp_vals, bounds_error=False)
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


def estimate_clade_frequencies(node, all_dates, tip_to_date_index):
	# extract time points and the subset of observations that fall in the clade.
	tps = all_dates[tip_to_date_index[node.parent_node.tips]]
	start_index = max(0,np.searchsorted(tps, time_interval[0]))
	stop_index = min(np.searchsorted(tps, time_interval[1]), all_dates.shape[0]-1)
	tps = tps[start_index:stop_index]
	obs = np.in1d(tps, all_dates[tip_to_date_index[node.tips]])
	if obs.sum()>50:
		print node.num_date, len(node.tips)

		# make n pivots a year
		pivots = get_pivots(tps[0], tps[1])
		fe = frequency_estimator(zip(tps, obs), pivots=pivots, stiffness=2.0, logit=True)
		fe.learn()

		node.freq = node.parent_node.freq * logit_inv(fe.pivot_freq)
		node.logit_freq = logit_transform(node.freq)
		for child in node.child_nodes():
			estimate_clade_frequencies(child, all_dates, tip_to_date_index)
	else:
		node.freq=None
		node.logit_freq=None

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

	tree.seed_node.pivots = get_pivots(time_interval[0], time_interval[1])
	tree.seed_node.freq = np.ones_like(tree.seed_node.pivots)
	for child in tree.seed_node.child_nodes():
		estimate_clade_frequencies(child, all_dates, reverse_order)


def estimate_genotype_frequency(tree, gt, time_interval=None, region = None):
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
		if region is not None:
			good_region = node.region==region
		else:
			good_region = True

		if good_time and good_region:
			all_dates.append(node.num_date)
			observations.append(is_gt)

	all_dates = np.array(all_dates)
	leaf_order = np.argsort(all_dates)
	tps = all_dates[leaf_order]
	obs = np.array(observations)[leaf_order]
	# define pivots and estimate
	pivots = get_pivots(tps[0], tps[1])
	print "number of pivots:", len(pivots)
	fe = frequency_estimator(zip(tps, obs), pivots=pivots, stiffness=flu_stiffness, logit=True, verbose = 0)
	fe.learn()
	return fe.frequency_estimate, (tps,obs)


def determine_clade_frequencies(tree):
	'''
	loop over different clades and determine their frequencies
	returns a dictionary with clades:frequencies
	'''
	import matplotlib.pyplot as plt
	xpol_pivots = get_extrapolation_pivots(time_interval[1], dt=0.5)
	clade_frequencies = {"pivots":list(get_pivots(time_interval[0], time_interval[1])),
						 "xpol_pivots":list(xpol_pivots)}

	for ci, (clade_name, clade_gt) in enumerate(clade_designations.iteritems()):
		print clade_name, clade_gt
		freq, (tps, obs) = estimate_genotype_frequency(tree, [(pos+15, aa) for pos, aa in clade_gt], time_interval)
		xpol_freq = extrapolation(freq, xpol_pivots)
		clade_frequencies[clade_name] = [list(freq.y), list(logit_inv(freq.y))]
		clade_frequencies['xpol_'+clade_name] = [list(xpol_freq), list(logit_inv(xpol_freq))]

		grid_tps = np.linspace(time_interval[0], time_interval[1], 100)
		plt.plot(grid_tps, logit_inv(freq(grid_tps)), label=clade_name, lw=2, c=cols[ci%len(cols)])
		plt.plot(xpol_pivots, logit_inv(xpol_freq),lw=2, ls='--', c=cols[ci%len(cols)])
		r_avg = running_average(obs, 100)
		plt.plot(tps, r_avg, c=cols[ci%len(cols)])
	plt.legend()
	ticloc = np.arange(time_interval[0], int(time_interval[1])+1,1)
	plt.xticks(ticloc, map(str, ticloc))
	plt.xlim([time_interval[0], time_interval[1]+1])
	return clade_frequencies


def determine_mutation_frequencies(tree):
	import matplotlib.pyplot as plt
	from collections import defaultdict
	from itertools import izip

	mut_counts = defaultdict(int)
	ref_seq = tree.seed_node.aa_seq
	total_leaf_count = 0
	for node in tree.leaf_iter():
		if (node.num_date>= time_interval[0]) and (node.num_date<time_interval[1]):
			total_leaf_count+=1
			for pos, (a,b) in enumerate(izip(ref_seq, node.aa_seq)):
				if a!=b: mut_counts[(pos, b)]+=1

	plt.figure()
	xpol_pivots = get_extrapolation_pivots(time_interval[1], dt=0.5)
	mutation_frequencies = {"pivots":list(get_pivots(time_interval[0], time_interval[1])),
						 "xpol_pivots":list(xpol_pivots)}
	for mi, (mut, count) in enumerate(mut_counts.iteritems()):
		if count>50 and count<total_leaf_count-50:
			print mut, count
			freq, (tps, obs) = estimate_genotype_frequency(tree, [mut], time_interval)
			mutation_frequencies[str(mut[0]-15)+mut[1]] = [list(freq.y), list(logit_inv(freq.y))]
			xpol_freq = extrapolation(freq, xpol_pivots)
			mutation_frequencies['xpol_'+str(mut[0]-15)+mut[1]] = [list(xpol_freq), list(logit_inv(xpol_freq))]

			grid_tps = np.linspace(time_interval[0], time_interval[1], 100)
			plt.plot(grid_tps, logit_inv(freq(grid_tps)), label=str(mut[0]-15)+mut[1], lw=2, c=cols[mi%len(cols)])
			plt.plot(xpol_pivots, logit_inv(xpol_freq), lw=2, ls='--', c=cols[mi%len(cols)])
			r_avg = running_average(obs, 100)
			plt.plot(tps, r_avg, c=cols[mi%len(cols)])
	plt.legend()
	ticloc = np.arange(time_interval[0], int(time_interval[1])+1,1)
	plt.xticks(ticloc, map(str, ticloc))
	plt.xlim([time_interval[0], time_interval[1]+1])
	return mutation_frequencies




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

	tree_fname='data/tree_refine.json'
	tree =  json_to_dendropy(read_json(tree_fname))

	print "--- "+"determining clade frequencies "  + time.strftime("%H:%M:%S") + " ---"
	clade_frequencies = determine_clade_frequencies(tree)
	out_fname = 'data/clade_frequencies.json'
	write_json(clade_frequencies, out_fname)
	plt.savefig('data/clade_frequencies.pdf')

	print "--- "+"determining mutation frequencies "  + time.strftime("%H:%M:%S") + " ---"
	mutation_frequencies = determine_mutation_frequencies(tree)
	out_fname = 'data/mutation_frequencies.json'
	write_json(mutation_frequencies, out_fname)
	plt.savefig('data/mutation_frequencies.pdf')

	print "--- "+"adding frequencies to tree "  + time.strftime("%H:%M:%S") + " ---"
	out_fname = 'data/tree_frequencies.json'
	estimate_tree_frequencies(tree)
	write_json(dendropy_to_json(tree.seed_node), out_fname)

if __name__=="__main__":
	#test()
	main()


