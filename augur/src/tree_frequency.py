# estimates clade frequencies using SMC

import time
from io_util import *
from tree_util import *
from date_util import *
import numpy as np

START_DATE = "2011-01-01"
DAY_STEP = 10

def observations_for_clade(node, dates):
	"""Takes a node and a list of dates and returns an observation list"""
	node_dates = get_dates(node)
	obs = []
	for d in dates:
		if d in node_dates:
			obs.append(1)
			node_dates.remove(d)
		else:
			obs.append(0)
	return obs

def obs_likelihood(obs, freq):
	"""Bernoulli probability of observing binary presence / absence given circulating frequency"""
	"""obs is 0 or 1 binary indicator for absence / presence"""
	return (1-freq)**(1-obs) * freq**obs

class Particle(object):
	"""An individual SMC particle.  Stores value and weight."""
	min_value = 0.00001
	max_value = 0.99999

	def __init__(self, obs):
		self.value = np.random.random_sample()	# between [0,1]
		self.reweight(obs)

	def __str__(self):
		return "{value: " + str(self.value) + ", weight: " + str(self.weight) + "}"

	def simulate_step(self, dt, sigma):
		"""Single time step forward"""
		w = np.random.normal()
		x0 = self.value
		x1 = x0 + x0 * (1-x0) * sigma * np.sqrt(dt) * w
		if x1 < self.min_value:
			x1 = self.min_value
		if x1 > self.max_value:
			x1 = self.max_value
		self.value = x1

	def simulate(self, t, dt, sigma):
		"""Simulate forward t years with timesteps dt"""
		"""Implements Euler-Maruyama method"""
		steps = int(np.ceil(t / dt))
		if steps > 0:
			dt = t / steps
			sqrtdt_sigma = np.sqrt(dt) * sigma
			for i in range(steps):
				x0 = self.value
				x1 = x0 + x0 * (1-x0) * sqrtdt_sigma * np.random.normal()
				if x1 < self.min_value:
					x1 = self.min_value
				if x1 > self.max_value:
					x1 = self.max_value
				self.value = x1

	def reweight(self, obs):
		self.weight = obs_likelihood(obs, self.value) + 0.000000001

	def likelihood(self, obs):
		return obs_likelihood(obs, self.value)

class Filter(object):
	"""An SMC particle filter, comprising n particles."""
	"""Takes a list of dates and a list of (0,1) observations."""
	"""Time in encoded in continuous units of years."""
	particle_count = 500
	timestep = 0.01
	sigma = 5

	def __init__(self, dates, observations):
		self.num_dates = map(string_to_numerical_date, dates)
		self.observations = observations
		self.end_num_date = numerical_date(datetime.date.today())
		first_obs = self.observations[0]
		self.particles = [Particle(first_obs) for i in range(Filter.particle_count)]
		self.means = []

	def __str__(self):
		string = "[\n"
		for particle in self.particles:
			string += "    " + str(particle) + "\n"
		string += "]\n"
		return string

	def values(self):
		return [p.value for p in self.particles]

	def weights(self):
		return [p.weight for p in self.particles]

	def set_values(self, new_values):
		for (p,v) in zip(self.particles, new_values):
			p.value = v

	def set_weights(self, new_weights):
		for (p,w) in zip(self.particles, new_weights):
			p.weight = w

	def simulate(self, t):
		"""Simulate forward t years"""
		for p in self.particles:
			p.simulate(t, self.timestep, self.sigma)

	def reweight(self, obs):
		"""Reweight particles according to observation"""
		for p in self.particles:
			p.reweight(obs)

	def resample(self):
		"""Resample particles according to weights"""
		total = sum(self.weights())
		nweights = [float(w) / total for w in self.weights()]
		new_values = np.random.choice(self.values(), size=self.particle_count, replace=True, p=nweights)
		self.set_values(new_values)
		self.set_weights([1] * self.particle_count)

	def update(self, t, obs):
		self.simulate(t)
		self.reweight(obs)
		self.resample()

	def run(self):
		steps = [j-i for i, j in zip(self.num_dates[:-1], self.num_dates[1:])]
		for (s, obs, t) in zip(steps, self.observations[1:], self.num_dates[1:]):
			self.update(s, obs)
			self.means.append(self.mean())
		final_step = self.end_num_date - self.num_dates[-1]
		self.simulate(final_step)

	def mean(self):
		return np.mean(self.values())

	def likelihood(self, obs):
		l = []
		for p in self.particles:
			l.append(p.likelihood(obs))
		return np.mean(l)

	def run_with_likelihood(self):
		steps = [j-i for i, j in zip(self.num_dates[:-1], self.num_dates[1:])]
		ll = 0
		for (s, obs) in zip(steps, self.observations[1:]):
			self.update(s, obs)
			ll += np.log(self.likelihood(obs))
		final_step = self.end_num_date - self.num_dates[-1]
		self.simulate(final_step)
		return ll

def estimate_likelihood():
	tree = read_json('data/tree.json')
	dates = get_dates(tree)
	nodes = [n for n in all_descendants(tree)]
	ll = 0
	for node in nodes[1:100]:
		observations = observations_for_clade(node, dates)
		filter = Filter(dates, observations)
		ll += filter.run_with_likelihood()
		node['frequency'] = filter.mean()
		print str(node['clade']) + ": " + str(node['frequency'])

	print "log likelihood: " + str(ll)

def set_node_frequency(node, dates):
	observations = observations_for_clade(node, dates)
	filter = Filter(dates, observations)
	filter.run()
	node['frequency'] = round(filter.mean(), 5)
	means = []
	for (d, m) in zip(dates, filter.means):
		timepoint = {}
		timepoint['date'] = d
		timepoint['frequency'] = round(m, 5)
		means.append(timepoint)
	means = means[1:-1:10]
	node['frequencies'] = means

def main():
	print "--- Frequencies at " + time.strftime("%H:%M:%S") + " ---"

	tree = read_json('data/tree_clean.json')
	dates = get_dates(tree)

	nodes = [n for n in all_descendants(tree)]

	start = time.clock()
	for node in nodes:
		set_node_frequency(node, dates)
		print str(node['clade']) + ": " + str(node['frequency'])
		print str(node['clade']) + ": " + str(node['frequencies'])
	print "time: " + str(time.clock() - start)

	write_json(tree, "data/tree_freq.json")

if __name__ == "__main__":
	main()
