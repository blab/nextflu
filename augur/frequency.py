# estimates clade frequencies using SMC

import time
from io_util import *
from tree_util import *
from date_util import *
import numpy as np

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

class Particle:
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
		steps = t / dt
		array = np.linspace(0, t, num=steps)
		for step in array:
			self.simulate_step(step, sigma)

	def reweight(self, obs):
		self.weight = obs_likelihood(obs, self.value) + 0.000000001

class Filter:
	"""An SMC particle filter, comprising n particles."""
	"""Takes a list of dates and a list of (0,1) observations."""
	"""Time in encoded in continuous units of years."""
	particle_count = 100
	timestep = 0.001
	sigma = 1

	def __init__(self, dates, observations):
		self.num_dates = map(string_to_numerical_date, dates)
		self.observations = observations
		self.end_num_date = numerical_date(datetime.date.today())
		first_obs = self.observations[0]
		self.particles = [Particle(first_obs) for i in range(Filter.particle_count)]

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
		start = time.clock()
		self.simulate(t)
		print "simulate: " + str(time.clock() - start)
		start = time.clock()
		self.reweight(obs)
		print "reweight: " + str(time.clock() - start)
		start = time.clock()
		self.resample()
		print "resample: " + str(time.clock() - start)

	def run(self):
		times = [j-i for i, j in zip(self.num_dates[:-1], self.num_dates[1:])]
		for (t, obs) in zip(times, self.observations[1:]):
			self.update(t, obs)
		final_step = self.end_num_date - self.num_dates[-1]
		self.simulate(final_step)

	def mean(self):
		return np.mean(self.values())

def main():
	print "--- Frequencies at " + time.strftime("%H:%M:%S") + " ---"

	tree = read_json('tree.json')
	dates = get_dates(tree)

	for node in all_descendants(tree):
		observations = observations_for_clade(node, dates)
		filter = Filter(dates, observations)
		filter.run()
		node['frequency'] = filter.mean()
		print str(node['clade']) + ": " + str(node['frequency'])

if __name__ == "__main__":
    main()
