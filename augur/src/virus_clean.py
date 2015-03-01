# clean sequences after alignment, criteria based on sequences
# make inline with canonical ordering (no extra gaps)

import os, datetime, time, re
from itertools import izip
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from scipy import stats
import numpy as np
from io_util import *

class virus_clean(object):
	"""docstring for virus_clean"""
	def __init__(self):
		pass

	def remove_insertions(self):
		outgroup_ok = np.array(self.sequence_lookup[self.outgroup['strain']])!='-'
		for seq in self.viruses:
			seq.seq = Seq("".join(np.array(seq.seq)[outgroup_ok]).upper())

	def clean_gaps(self):
		self.viruses = filter(lambda x: '-' in x.seq, self.viruses)

	def clean_ambiguous(self):
		for v in self.viruses:
			v.seq = Seq(re.sub(r'[BDEFHIJKLMNOPQRSUVWXYZ]', '-',str(v.seq)))

	def unique_date(self):
		from date_util import numerical_date
		og = self.sequence_lookup[self.outgroup['strain']]
		og.num_date = numerical_date(og.date)
		for ii, v in enumerate(self.viruses):
			v.num_date = numerical_date(v.date) + 1e-7*(ii+1)

	def times_from_outgroup(self):
		self.unique_date()
		outgroup_date = self.sequence_lookup[self.outgroup['strain']].num_date
		return np.array([x.num_date-outgroup_date for x in self.viruses])

	def distance_from_outgroup(self):
		from seq_util import hamming_distance
		outgroup_seq = self.sequence_lookup[self.outgroup['strain']].seq
		return np.array([hamming_distance(x.seq, outgroup_seq) for x in self.viruses])

	def clean_distances(self, n_std = 5):
		"""Remove viruses that don't follow a loose clock """
		times = self.times_from_outgroup()
		distances = self.distance_from_outgroup()
		slope, intercept, r_value, p_value, std_err = stats.linregress(times, distances)
		residuals = slope*times - distances
		r_sd = residuals.std()
		if self.verbose:
			print "\tslope: " + str(slope)
			print "\tr: " + str(r_value)
			print "\tresiduals sd: " + str(r_sd)
		new_viruses = []
		for (v,r) in izip(self.viruses,residuals):		# filter viruses more than 5 sds up or down
			if np.abs(r)<n_std * r_sd or v.id == self.outgroup["strain"]:
				new_viruses.append(v)
			else:
				if self.verbose>1:
					print "\t\tresidual:", r, "\nremoved ",v.strain
		self.viruses = MultipleSeqAlignment(new_viruses)

	def clean_generic(self):
		print "Number of viruses before cleaning:",len(self.viruses)
		self.remove_insertions()
		self.clean_ambiguous()
		self.clean_distances()
		self.viruses.sort(key=lambda x:x.num_date)
		print "Number of viruses after outlier filtering:",len(self.viruses)
