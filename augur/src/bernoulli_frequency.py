# estimates clade frequencies using SMC

from scipy.interpolate import interp1d
import time
from io_util import *
from tree_util import *
from date_util import *
import numpy as np

pc=1e-3
class frequency_estimator(object):

	def __init__(self, observations, npivots = 10, stiffness = 2000.0):
		self.tps = np.array([x[0] for x in observations])
		self.obs = np.array([x[1]>0 for x in observations])
		self.stiffness = stiffness

		# make sure they are sorted
		tmp = np.argsort(self.tps)
		self.tps = self.tps[tmp]
		self.obs = self.obs[tmp]

		self.pivot_tps = np.linspace(self.tps[0], self.tps[-1], npivots)
		self.pivot_freq = np.mean(self.obs)*np.ones(npivots)

	def stiffLH(self, pivots):
		return -0.5*self.stiffness*np.sum(np.diff(pivots)**2/np.diff(self.pivot_tps))

	def logLH(self, pivots):
		freq = interp1d(self.pivot_tps, pivots)
		estfreq = freq(self.tps)
		LH = self.stiffLH(pivots) + np.sum(np.log(np.maximum(estfreq[self.obs],pc))) + np.sum(np.log(np.maximum(1-estfreq[~self.obs], pc)))
		print LH
		return -LH +np.sum(pivots<0)+np.sum(pivots>1)

	def learn(self):
		from scipy.optimize import fmin as minimizer
		self.pivot_freq = minimizer(self.logLH, self.pivot_freq)

def main():
	print "--- Frequencies at " + time.strftime("%H:%M:%S") + " ---"
	from scipy.interpolate import UnivariateSpline
	import matplotlib.pyplot as plt
	tree_fname = 'data/tree_refine_10y_50v.json'
	tree =  json_to_dendropy(read_json(tree_fname))
	dates = []
	for node in tree.postorder_node_iter():
		if node.is_leaf():
			node.date = datetime.datetime(*map(int, node.date.split('-')))
		else:
			node.date = min([c.date for c in node.child_nodes()])
		dates.append(node.date)
	dates.sort()
	ordinal_dates = [d.toordinal() for d in dates]
	dt=10
	y,x = np.histogram(ordinal_dates, bins = np.arange(min(ordinal_dates), max(ordinal_dates)+dt, dt))
	bc =  0.5*(x[:-1]+x[1:])
	sampling_intensity = UnivariateSpline(bc,y, w=1.0/np.sqrt(y+5)) 
	plt.plot(bc, y)
	plt.plot(bc, sampling_intensity(bc))

#	for node in tree.postorder_node_iter():
#		print [c.date.toordinal() for c in node.child_nodes()]

	return tree, dates, sampling_intensity

if __name__ == "__main__":
	fe = frequency_estimator(zip(np.arange(100), np.random.rand(100)<np.exp(-np.arange(100.0)/50)))
	fe.learn()
