# estimates clade frequencies using SMC

from scipy.interpolate import interp1d
import time
from io_util import *
from tree_util import *
from date_util import *
import numpy as np

pc=1e-2

class frequency_estimator(object):

	def __init__(self, observations, npivots = 10, stiffness = 20.0, logit=False):
		self.tps = np.array([x[0] for x in observations])
		self.obs = np.array([x[1]>0 for x in observations])
		self.stiffness = stiffness
		self.interolation_type = 'cubic'
		self.logit = logit
		# make sure they are searchsorted
		tmp = np.argsort(self.tps)
		self.tps = self.tps[tmp]
		self.obs = self.obs[tmp]

		# generate a useful initital case from a running average of the counts
		ws=10
		self.pivot_tps = np.linspace(self.tps[0], self.tps[-1], npivots)
		tmp_vals = np.convolve(np.ones(ws, dtype=float)/ws, self.obs, mode='same')
		# fix the edges. using mode='same' assumes zeros outside the range
		tmp_vals[:ws-1]*=ws/np.arange(1, ws)
		tmp_vals[-ws+1:]*=ws/np.arange(ws-1,0,-1)
		# calculate interpolated frequences as initial pivots
		tmp_interpolator = interp1d(self.tps, tmp_vals)
		if self.logit:
			freq = tmp_interpolator(self.pivot_tps)
			self.pivot_freq = np.log(np.maximum(pc,freq)/np.maximum(pc, 1-freq))
		else:
			self.pivot_freq = tmp_interpolator(self.pivot_tps)


	def stiffLH(self, pivots):
		if self.logit: # if logit, convert to frequencies
			logit_freq = np.exp(pivots)
			freq = logit_freq/(1+logit_freq)
			dfreq = np.diff(pivots)
		else:
			dfreq = np.diff(pivots)
			freq = pivots
		# return wright fisher diffusion likelihood for frequency change. 
		return -0.25*self.stiffness*np.sum(dfreq**2/np.diff(self.pivot_tps)/
											(freq[:-1]+pc)*(1-freq[:-1]+pc))


	def logLH(self, pivots):
		freq = interp1d(self.pivot_tps, pivots, kind=self.interolation_type)
		if self.logit: # if logit, convert to frequencies
			logit_freq = np.exp(freq(self.tps))
			estfreq = logit_freq/(1+logit_freq)
		else:
			estfreq = freq(self.tps)
		LH = self.stiffLH(pivots) + np.sum(np.log(np.maximum(estfreq[self.obs],pc))) + np.sum(np.log(np.maximum(1-estfreq[~self.obs], pc)))
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

	# make six pivots a year
	npivots = int(6*(tps[-1]-tps[0]))
	fe = frequency_estimator(zip(tps, obs), npivots=npivots, stiffness=2.0, logit=True)
	fe.learn()
	# return the final interpolation object
	return fe.frequency_estimate

def estimate_tree_frequencies(tree):
	all_dates = []
	leaf_count = 0
	# loop over all nodes, make time ordered lists of tips
	for node in tree.postorder_node_iter():
		if not node.is_leaf():
			node.date = min([c.date for c in node.child_nodes()])
		node.num_date = numerical_date(string_to_date(node.date))+1e-7*leaf_count
		tmp_tips = []
		if node.is_leaf():
			all_dates.append(node.num_date)
			node.tip_index = leaf_count
			leaf_count+=1
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



def test():
	import matplotlib.pyplot as plt
	tps = np.sort(100 * np.random.uniform(size=100))
	freq = [0.1]
	stiffness=1000
	s=-0.02
	for dt in np.diff(tps):
		freq.append(freq[-1]*np.exp(-s*dt)+np.sqrt(2*np.max(0,freq[-1]*(1-freq[-1]))*dt/stiffness)*np.random.normal())
	obs = np.random.uniform(size=tps.shape)<freq
	fe = frequency_estimator(zip(tps, obs), npivots=10, stiffness=stiffness)
	fe.learn()
	plt.figure()
	plt.plot(tps, freq, 'o', label = 'actual frequency')
	freq = fe.frequency_estimate(fe.tps)
	plt.plot(fe.tps, freq, '-', label='interpolation')
	plt.plot(tps, (2*obs-1)*0.05, 'o')
	plt.plot(tps[obs], 0.05*np.ones(np.sum(obs)), 'o', c='r', label = 'observations')
	plt.plot(tps[~obs], -0.05*np.ones(np.sum(1-obs)), 'o', c='g')
	plt.plot(tps, np.zeros_like(tps), 'k')
	ws=20
	plt.plot(fe.tps[ws/2:-ws/2+1], np.convolve(np.ones(ws, dtype=float)/ws, obs, mode='valid'), 'r', label = 'running avg')
	plt.legend(loc=2)

def main():
	# load tree
	import matplotlib.pyplot as plt
	from io_util import read_json
	from tree_util import json_to_dendropy, to_Biopython, color_BioTree_by_attribute
	from Bio import Phylo
	tree_fname='data/tree_refine_10y_50v.json'
	tree =  json_to_dendropy(read_json(tree_fname))
	estimate_tree_frequencies(tree)
	for node in tree.postorder_node_iter():
		if node.freq_est is not None:
			tps = np.linspace(node.freq_est.x[0], node.freq_est.x[-1],100)
			logit_freq = np.exp(node.freq_est(tps))
			plt.plot(tps, logit_freq/(1+logit_freq))
	return tree

if __name__=="__main__":
	#test()
	main()