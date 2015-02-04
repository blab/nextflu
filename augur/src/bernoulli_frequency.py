# estimates clade frequencies using SMC

from scipy.interpolate import interp1d
import time
from io_util import *
from tree_util import *
from date_util import *
import numpy as np

pc=1e-2
class frequency_estimator(object):

	def __init__(self, observations, npivots = 10, stiffness = 20.0):
		self.tps = np.array([x[0] for x in observations])
		self.obs = np.array([x[1]>0 for x in observations])
		self.stiffness = stiffness
		self.interolation_type = 'cubic'
		# make sure they are sorted
		tmp = np.argsort(self.tps)
		self.tps = self.tps[tmp]
		self.obs = self.obs[tmp]

		self.pivot_tps = np.linspace(self.tps[0], self.tps[-1], npivots)
		self.pivot_freq = np.mean(self.obs)*np.ones(npivots)

	def stiffLH(self, pivots):
		return -0.25*self.stiffness*np.sum(np.diff(pivots)**2/np.diff(self.pivot_tps)/
											(pivots[:-1]+pc)*(1-pivots[:-1]+pc))


	def logLH(self, pivots):
		freq = interp1d(self.pivot_tps, pivots, kind=self.interolation_type)
		estfreq = freq(self.tps)
		LH = self.stiffLH(pivots) + np.sum(np.log(np.maximum(estfreq[self.obs],pc))) + np.sum(np.log(np.maximum(1-estfreq[~self.obs], pc)))
		#LH =  np.sum(np.log(np.maximum(estfreq[self.obs],pc))) + np.sum(np.log(np.maximum(1-estfreq[~self.obs], pc)))
		return -LH/len(self.obs)+100000*(np.sum(pivots<0)+np.sum(pivots>1))

	def learn(self):
		from scipy.optimize import fmin as minimizer
		self.pivot_freq = minimizer(self.logLH, self.pivot_freq)
		self.frequency_estimate = interp1d(self.pivot_tps, self.pivot_freq, kind=self.interolation_type)


if __name__ == "__main__":
	import matplotlib.pyplot as plt
	tps = np.sort(100 * np.random.uniform(size=100))
	freq = [0.1]
	stiffness=1000
	s=-0.02
	for dt in np.diff(tps):
		freq.append(freq[-1]*np.exp(-s*dt)+np.sqrt(2*freq[-1]*(1-freq[-1])*dt/stiffness)*np.random.normal())
	obs = np.random.uniform(size=tps.shape)<freq
	fe = frequency_estimator(zip(tps, obs), npivots=10, stiffness=stiffness)
	fe.learn()
	plt.figure()
	plt.plot(tps, freq, 'o', label = 'actual frequency')
	plt.plot(fe.tps, fe.frequency_estimate(fe.tps), '-', label='interpolation')
	plt.plot(tps, (2*obs-1)*0.05, 'o')
	plt.plot(tps[obs], 0.05*np.ones(np.sum(obs)), 'o', c='r', label = 'observations')
	plt.plot(tps[~obs], -0.05*np.ones(np.sum(1-obs)), 'o', c='g')
	plt.plot(tps, np.zeros_like(tps), 'k')
	ws=20
	plt.plot(fe.tps[ws/2:-ws/2+1], np.convolve(np.ones(ws, dtype=float)/ws, obs, mode='valid'), 'r', label = 'running avg')
	plt.legend(loc=2)