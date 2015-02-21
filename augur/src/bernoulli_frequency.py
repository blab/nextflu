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
time_interval = (2012.0, 2015.1)
flu_stiffness = 10.0
pivots_per_year = 12.0
inertia = 0.7    # fraction of previous frequency changes that is carried over
window_size = 20 # smooting window
tol = 1e-4
reg = 1e-6
debug = False

clade_designations = { "3c3.a":[(128,'A'), (142,'G'), (159,'S')],
					   "3c3":  [(128,'A'), (142,'G'), (159,'F')],
					   "3c2.a":[(144,'S'), (159,'Y'), (225,'D'), (311,'H'),(489,'N')],
					   "3c2":  [(144,'N'), (159,'F'),(225,'N'), (489,'N')],
						}

region_names = ['Europe', 'India', 'NorthAmerica', 'SouthAmerica', 'Africa', 
			'JapanKorea', 'Oceania', 'China', 'WestAsia', 'SoutheastAsia']

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
	logit_freq[logit_freq<-20]=-20
	logit_freq[logit_freq>20]=20
	tmp_freq = np.exp(logit_freq)
	return tmp_freq/(1.0+tmp_freq)

def pq(p):
	return p*(1-p)

def logit_regularizer(logit_freqs):
	return reg*np.mean(np.abs(8-np.abs(logit_freqs))) # penalize too large or too small pivots

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
		return np.array([ep(freq_interp,tmp_x) for tmp_x in x])

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
			self.final_pivot_tps = get_pivots(self.tps[0], self.tps[1])
		elif np.isscalar(pivots):
			self.final_pivot_tps = np.linspace(self.tps[0], self.tps[-1], pivots)
		else:
			self.final_pivot_tps = pivots

	def initial_guess(self, pivots, ws=50):
		# generate a useful initital case from a running average of the counts
		tmp_vals = running_average(self.obs, ws)
		tmp_interpolator = interp1d(self.tps, tmp_vals, bounds_error=False, fill_value = -1)
		pivot_freq = tmp_interpolator(pivots)
		pivot_freq[pivots<=tmp_interpolator.x[0]] = tmp_vals[0]
		pivot_freq[pivots>=tmp_interpolator.x[-1]] = tmp_vals[-1]
		pivot_freq = fix_freq(pivot_freq, pc)
		if self.logit:
			self.pivot_freq = logit_transform(pivot_freq)

		return pivot_freq


	def stiffLH(self, pivots):
		if self.logit: # if logit, convert to frequencies
			freq = logit_inv(pivots)
		else:
			freq = pivots
		dfreq = np.diff(freq)
		dt = np.diff(self.pivot_tps)
		tmp_freq = fix_freq(freq,dfreq_pc)
		# return wright fisher diffusion likelihood for frequency change. 
		# return -0.25*self.stiffness*np.sum(dfreq**2/np.diff(self.pivot_tps)/pq(fix_freq(freq[:-1],dfreq_pc)))
		return -0.25*self.stiffness*(np.sum((dfreq[1:] - inertia*dfreq[:-1])**2/(dt[1:]*pq(tmp_freq[1:-1])))
									+dfreq[0]**2/(dt[0]*pq(tmp_freq[0])))


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
			return -LH/len(self.obs) + logit_regularizer(pivots)
		else:
			return -LH/len(self.obs)+100000*(np.sum((pivots<0)*np.abs(pivots))+np.sum((pivots>1)*np.abs(pivots-1)))


	def learn(self):
		from scipy.optimize import fmin_powell as minimizer
		self.final_pivot_freq = self.initial_guess(self.final_pivot_tps, ws=2*(min(50,len(self.obs))//2))
		if self.verbose:
			print "Initial pivots:", self.final_pivot_freq
		steps= [4,2,1]
		for step in steps:
			# subset the pivots, if the last point is not included, attach it
			self.pivot_tps = self.final_pivot_tps[::step]
			if self.pivot_tps[-1]!=self.final_pivot_tps[-1]:
				self.pivot_tps = np.concatenate((self.pivot_tps, self.final_pivot_tps[-1:]))
				attached_endpoint=True
			else:
				attached_endpoint=False				

			if step==steps[0]: # in first iteration, take intial frequency from the running average calc above
				self.pivot_freq = self.final_pivot_freq[::step]
				if attached_endpoint:
					self.pivot_freq = np.concatenate((self.pivot_freq, self.final_pivot_freq[-1:]))
			else: #otherwise interpolate from previous iterations
				self.pivot_freq = self.frequency_estimate(self.pivot_tps)

			# determine the optimal pivot freqquencies
			self.pivot_freq = minimizer(self.logLH, self.pivot_freq, ftol = tol, xtol = tol, disp = self.verbose>0)
			# instantiate an interpolation object based on the optimal frequency pivots
			self.frequency_estimate = interp1d(self.pivot_tps, self.pivot_freq, kind=self.interolation_type, bounds_error=False)
			if self.verbose: print "neg logLH using",len(self.pivot_tps),"pivots:", self.logLH(self.pivot_freq)
		self.final_pivot_freq=self.pivot_freq			


def estimate_sub_frequencies(node, all_dates, tip_to_date_index, threshold=50, region_name="global"):
	# extract time points and the subset of observations that fall in the clade.
	tps = all_dates[tip_to_date_index[node.tips]]
	start_index = max(0,np.searchsorted(tps, time_interval[0]))
	stop_index = min(np.searchsorted(tps, time_interval[1]), all_dates.shape[0]-1)
	tps = tps[start_index:stop_index]
	# we estimate frequencies of subclades, they will be multiplied by the 
	# frequency of the parent node and corrected for the frequency of sister clades 
	# already fit
	frequency_left = np.array(node.freq[region_name])
	if frequency_left is None:
		import pdb; pdb.set_trace()
	ci=0
	# need to resort, since the clade size order might differs after subsetting to regions
	children_by_size = sorted(node.child_nodes(), key = lambda x:len(x.tips), reverse=True)
	for child in children_by_size[:-1]: # clades are ordered by decreasing size
		if len(child.tips)<threshold: # skip tiny clades
			break
		else:
			if debug: print child.num_date, len(child.tips)

			obs = np.in1d(tps, all_dates[tip_to_date_index[child.tips]])

			# make n pivots a year, interpolate frequencies
			# FIXME: adjust stiffness to total number of observations in a more robust manner
			pivots = get_pivots(tps[0], tps[1])
			fe = frequency_estimator(zip(tps, obs), pivots=pivots, stiffness=flu_stiffness*len(all_dates)/2000.0, logit=True)
			fe.learn()

			try:
				# assign the frequency vector to the node
				child.freq[region_name] = frequency_left * logit_inv(fe.pivot_freq)
				child.logit_freq[region_name] = logit_transform(child.freq[region_name])
			except:
				import pdb; pdb.set_trace()

			# update the frequency remaining to be explained and prune explained observations
			frequency_left *= (1.0-logit_inv(fe.pivot_freq))
			tps_left = np.ones_like(tps,dtype=bool)
			tps_left[obs]=False # prune observations from clade
			tps = tps[tps_left]
		ci+=1

	# if the above loop finished assign the frequency of the remaining clade to the frequency_left
	if ci>0 and ci==len(node.child_nodes())-1:
		last_child = node.child_nodes()[-1]
		last_child.freq[region_name] = frequency_left
		last_child.logit_freq[region_name] = logit_transform(last_child.freq[region_name])
	else:  # broke out of loop because clades too small. 
		for child in children_by_size[ci:]: # assign freqs of all remaining clades to None.
			child.freq[region_name] = None
			child.logit_freq[region_name] = None

	# recursively repeat for subclades
	for child in node.child_nodes():
		estimate_sub_frequencies(child, all_dates, tip_to_date_index, threshold, region_name)

def estimate_tree_frequencies(tree, threshold = 20, regions=None, region_name = None):
	'''
	loop over nodes of the tree and estimate frequencies of all clade above a certain size
	'''
	all_dates = []
	# loop over all nodes, make time ordered lists of tips, restrict to the specified regions
	tip_index_region_specific = 0
	for node in tree.postorder_node_iter():
		tmp_tips = []
		if node.is_leaf():
			if regions is None or node.region in regions:
				all_dates.append(node.num_date)
				tmp_tips.append((tip_index_region_specific, node.num_date))
				tip_index_region_specific +=1
		for child in node.child_nodes():
			tmp_tips.extend(child.tips)
		node.tips = np.array([x for x in sorted(tmp_tips, key = lambda x:x[1] )])
		if not hasattr(node, "freq"): node.freq = {}
		if not hasattr(node, "logit_freq"): node.logit_freq = {}
		if not hasattr(node, "virus_count"): node.virus_count = {}

	# erase the dates from the tip lists and cast to int such that they can be used for indexing
	for node in tree.postorder_node_iter():
		if len(node.tips.shape)==2:
			node.tips = np.array(node.tips[:,0], dtype=int)
		else:
			node.tips = np.array([], dtype=int)

	# sort the dates and provide a reverse ordering as a mapping of tip indices to dates
	all_dates = np.array(all_dates)
	leaf_order = np.argsort(all_dates)
	reverse_order = np.argsort(leaf_order)
	all_dates = all_dates[leaf_order]

	if regions is None: 
		region_name="global"
	elif region_name is None:
		region_name = ",".join(regions)
	# set the frequency of the root node to 1, the logit frequency to a large value
	tree.seed_node.pivots = get_pivots(time_interval[0], time_interval[1])
	tree.seed_node.virus_count[region_name] = np.histogram(all_dates, bins = tree.seed_node.pivots)
	tree.seed_node.freq[region_name] = np.ones_like(tree.seed_node.pivots)
	tree.seed_node.logit_freq[region_name] = 10*np.ones_like(tree.seed_node.pivots)

	# start estimating frequencies of subclades recursively
	estimate_sub_frequencies(tree.seed_node, all_dates, reverse_order, threshold = threshold, region_name = region_name)


def estimate_genotype_frequency(tree, gt, time_interval=None, regions = None, relevant_pos = None):
	'''
	estimate the frequency of a particular genotype specified 
	gt   --		[(position, amino acid), ....]
	'''
	all_dates = []
	observations = []
	total_leaf_count = 0
	for node in tree.leaf_iter():
		total_leaf_count+=1
		if isinstance(gt, basestring):
			if relevant_pos is None:
				is_gt = gt==node.aa_seq
			else:
				is_gt = gt==reduce_genotype(node.aa_seq, relevant_pos)
		else:
			is_gt = all([node.aa_seq[pos]==aa for pos, aa in gt])

		if time_interval is not None:
			good_time = (node.num_date>= time_interval[0]) and (node.num_date<time_interval[1])
		else:
			good_time=True
		if regions is not None:
			good_region = node.region in regions
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
	fe = frequency_estimator(zip(tps, obs), pivots=pivots, 
	               stiffness=flu_stiffness*float(len(observations))/total_leaf_count, 
                   logit=True, verbose = 0)
	fe.learn()
	return fe.frequency_estimate, (tps,obs)


def determine_clade_frequencies(tree, regions=None, plot=False):
	'''
	loop over different clades and determine their frequencies
	returns a dictionary with clades:frequencies
	'''
	import matplotlib.pyplot as plt
	xpol_pivots = get_extrapolation_pivots(time_interval[1], dt=0.5)
	clade_frequencies = {"pivots":list(get_pivots(time_interval[0], time_interval[1])),
						 "xpol_pivots":list(xpol_pivots)}

	for ci, (clade_name, clade_gt) in enumerate(clade_designations.iteritems()):
		print "estimating frequency of clade", clade_name, clade_gt
		freq, (tps, obs) = estimate_genotype_frequency(tree, [(pos+15, aa) for pos, aa in clade_gt], time_interval, regions)
		clade_frequencies[clade_name] = list(np.round(logit_inv(freq.y),3))
		if plot:
			grid_tps = np.linspace(time_interval[0], time_interval[1], 100)
			plt.plot(grid_tps, logit_inv(freq(grid_tps)), label=clade_name, lw=2, c=cols[ci%len(cols)])
			if debug:
				r_avg = running_average(obs, window_size)
				plt.plot(tps, r_avg, c=cols[ci%len(cols)])
	return clade_frequencies


def determine_mutation_frequencies(tree, regions=None, threshold=50, plot=False):
	'''
	determine the abundance of all single nucleotide variants and estimate the 
	frequency trajectory of the top 10, plot those optionally
	'''
	import matplotlib.pyplot as plt
	from collections import defaultdict
	from itertools import izip

	mut_counts = defaultdict(int)
	relevant_pos = defaultdict(int)
	ref_seq = tree.seed_node.aa_seq
	total_leaf_count = 0
	for node in tree.leaf_iter():
		if (node.num_date>= time_interval[0]) and (node.num_date<time_interval[1]):
			total_leaf_count+=1
			for pos, (a,b) in enumerate(izip(ref_seq, node.aa_seq)):
				mut_counts[(pos, b)]+=1


	mutation_frequencies = {"pivots":list(get_pivots(time_interval[0], time_interval[1]))}
	for mi, mut in enumerate(sorted(mut_counts.keys())):
		count = mut_counts[mut]
		if count>threshold or count<total_leaf_count-threshold:
			print "estimating freq of ", mut, "total count:", count
			freq, (tps, obs) = estimate_genotype_frequency(tree, [mut], time_interval, regions)
			mutation_frequencies[str(mut[0]-15)+mut[1]] = list(np.round(logit_inv(freq.y),3))
#			xpol_freq = fix_freq(extrapolation(interp1d(freq.x, logit_inv(freq.y)), xpol_pivots), pc)
#			mutation_frequencies['xpol_'+str(mut[0]-15)+mut[1]] = list(np.round(xpol_freq,3))
			if plot:
				grid_tps = np.linspace(time_interval[0], time_interval[1], 100)
				plt.plot(grid_tps, logit_inv(freq(grid_tps)), label=str(mut[0]-15)+mut[1], lw=2, c=cols[mi%len(cols)])
				if debug:
					r_avg = running_average(obs, window_size)
					plt.plot(tps, r_avg, c=cols[mi%len(cols)])

	return mutation_frequencies

def reduce_genotype(gt, pos):
	rgt = np.zeros(len(gt), dtype='S1')
	rgt[:] = '.'
	for p in pos:
		rgt[p] = gt[p]
	return "".join(rgt)

def add_genotype_at_pos(tree, positions):
	for node in tree.postorder_node_iter():
		node.gt = "".join([node.aa_seq[pos] for pos in positions])


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




def all_mutations(tree, region_list, plot=False):
	import matplotlib.pyplot as plt
	mutation_frequencies = {}
	for region_label, regions in region_list:
		print "--- "+"determining mutation frequencies in "+region_label+ " "  + time.strftime("%H:%M:%S") + " ---"
		if plot:
			plt.figure("mutations in "+region_label, figsize = (12,7))
			if regions is not None: plt.title("Region: "+", ".join(regions))
		mutation_frequencies[region_label] = determine_mutation_frequencies(tree, regions, plot=plot, threshold = 5)
		if plot:
			plt.legend()
			ticloc = np.arange(time_interval[0], int(time_interval[1])+1,1)
			plt.xticks(ticloc, map(str, ticloc))
			plt.xlim([time_interval[0], time_interval[1]+1])
			plt.ylim([-0.05, 1.05])
			plt.grid()
			plt.savefig('data/mutation_frequencies'+region_label+'.pdf')

	relevant_pos = []
	for mut, freq in mutation_frequencies["global"].iteritems():
		if "pivot" not in mut:
			if np.max(freq)-np.min(freq)>0.1:
				pos = int(mut.split('_')[-1][:-1])+15
				relevant_pos.append(pos)
	relevant_pos = sorted(set(relevant_pos))

	return mutation_frequencies, relevant_pos

def all_genotypes(tree, region_list, relevant_pos):
	from collections import defaultdict
	total_leaf_count = len(tree.leaf_nodes())
	gt_counts = defaultdict(int)
	gt_frequencies = {}
	for node in tree.leaf_iter():
		gt_counts[reduce_genotype(node.aa_seq, relevant_pos)]+=1
	for region_label, regions in region_list:
		tc_in_region, sum_gts = 0,0
		gt_frequencies[region_label]={"pivots":list(get_pivots(time_interval[0], time_interval[1]))}
		print "--- "+"determining genotype frequencies "+region_label+ " "  + time.strftime("%H:%M:%S") + " ---"
		for gt, c in gt_counts.iteritems():
			print gt, c
			tmp_freq , (tps, obs) = estimate_genotype_frequency(tree, gt, regions=regions, relevant_pos = relevant_pos)
			gt_frequencies[region_label][gt] = list(logit_inv(tmp_freq.y))
			sum_gts+=np.sum(obs)
			tc_in_region = len(tps)
		print "region: ", sum_gts, "out of", tc_in_region, "(",total_leaf_count," in total)"
	return gt_frequencies


def all_clades(tree, region_list, plot=False):
	clade_frequencies = {}
	import matplotlib.pyplot as plt
	for region_label, regions in region_list:
		print "--- "+"determining clade frequencies "+region_label+ " "  + time.strftime("%H:%M:%S") + " ---"
		if plot:
			plt.figure("region "+region_label, figsize = (12,7))
			if regions is not None: plt.title("Region: "+", ".join(regions))
		clade_frequencies[region_label] = determine_clade_frequencies(tree, regions=regions, plot=plot)
		if plot:
			plt.legend()
			ticloc = np.arange(time_interval[0], int(time_interval[1])+1,1)
			plt.xticks(ticloc, map(str, ticloc))
			plt.xlim([time_interval[0], time_interval[1]+1])
			plt.ylim([-0.05, 1.05])
			plt.grid()
			plt.savefig('data/clade_frequencies_'+region_label+'.pdf')
	return clade_frequencies

def main(tree_fname = 'data/tree_refine.json'):
	# load tree
	from io_util import read_json
	plot = debug
	tree =  json_to_dendropy(read_json(tree_fname))
	region_list = [("global", None), ("NA", ["NorthAmerica"]), ("EU", ["Europe"]), 
			("AS", ["China", "SoutheastAsia", "JapanKorea"]), ("OC", ["Oceania"]) ]

	out_fname = 'data/genotype_frequencies.json'
	gt_frequencies = {}

	gt_frequencies["mutations"], relevant_pos = all_mutations(tree, region_list, plot)
	write_json(gt_frequencies, out_fname, indent=None)

	gt_frequencies["genotypes"] = all_genotypes(tree, region_list, relevant_pos)
	write_json(gt_frequencies, out_fname, indent=None)

	gt_frequencies["clades"] = all_clades(tree, region_list, plot)

	# round frequencies
	for gt_type in gt_frequencies:
		for reg in region_list:
			for gt in gt_frequencies[gt_type][reg[0]]:
				tmp = gt_frequencies[gt_type][reg[0]][gt]
				gt_frequencies[gt_type][reg[0]][gt] = [round(x,3) for x in tmp]

	write_json(gt_frequencies, out_fname, indent=None)

	tree_out_fname = 'data/tree_frequencies.json'
	for region_label, regions in region_list:
		print "--- "+"adding frequencies to tree "+region_label+ " "  + time.strftime("%H:%M:%S") + " ---"
		estimate_tree_frequencies(tree, threshold = 10, regions=regions, region_name=region_label)
	write_json(dendropy_to_json(tree.seed_node), tree_out_fname, indent=None)
	return tree_out_fname

if __name__=="__main__":
	#test()
	main()


