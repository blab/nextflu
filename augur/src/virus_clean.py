# clean sequences after alignment, criteria based on sequences
# make inline with canonical ordering (no extra gaps)

import os, datetime, time, re
from itertools import izip
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from scipy import stats
import numpy as np

class virus_clean(object):
	"""docstring for virus_clean"""
	def __init__(self,n_iqd  = 5, **kwargs):
		'''
		parameters
		n_std	-- number of interquartile distances accepted in molecular clock filter 
		'''
		self.n_iqd = n_iqd

	def remove_insertions(self):
		'''
		remove all columns from the alignment in which the outgroup is gapped
		'''
		outgroup_ok = np.array(self.sequence_lookup[self.outgroup['strain']])!='-'
		for seq in self.viruses:
			seq.seq = Seq("".join(np.array(seq.seq)[outgroup_ok]).upper())

	def clean_gaps(self):
		'''
		remove viruses with gaps -- not part of the standard pipeline
		'''
		self.viruses = filter(lambda x: '-' in x.seq, self.viruses)

	def clean_ambiguous(self):
		'''
		substitute all ambiguous characters with '-', 
		ancestral inference will interpret this as missing data
		'''
		for v in self.viruses:
			v.seq = Seq(re.sub(r'[BDEFHIJKLMNOPQRSUVWXYZ]', '-',str(v.seq)))

	def unique_date(self):
		'''
		add a unique numerical date to each leaf. uniqueness is achieved adding a small number
		'''
		from date_util import numerical_date
		og = self.sequence_lookup[self.outgroup['strain']]
		if hasattr(og, 'date'):
			try:
				og.num_date = numerical_date(og.date)
			except:
				print "cannot parse date"
				og.num_date="undefined";
		for ii, v in enumerate(self.viruses):
			if hasattr(v, 'date'):
				try:
					v.num_date = numerical_date(v.date, self.date_format['fields']) + 1e-7*(ii+1)
				except:
					print "cannot parse date"
					v.num_date="undefined";

	def times_from_outgroup(self):
		outgroup_date = self.sequence_lookup[self.outgroup['strain']].num_date
		return np.array([x.num_date-outgroup_date for x in self.viruses  if x.strain])

	def distance_from_outgroup(self):
		from seq_util import hamming_distance
		outgroup_seq = self.sequence_lookup[self.outgroup['strain']].seq
		return np.array([hamming_distance(x.seq, outgroup_seq) for x in self.viruses if x.strain])

	def clean_distances(self):
		"""Remove viruses that don't follow a loose clock """
		times = self.times_from_outgroup()
		distances = self.distance_from_outgroup()
		slope, intercept, r_value, p_value, std_err = stats.linregress(times, distances)
		residuals = slope*times + intercept - distances
		r_iqd = stats.scoreatpercentile(residuals,75) - stats.scoreatpercentile(residuals,25)
		if self.verbose:
			print "\tslope: " + str(slope)
			print "\tr: " + str(r_value)
			print "\tresiduals iqd: " + str(r_iqd)
		new_viruses = []
		for (v,r) in izip(self.viruses,residuals):
			# filter viruses more than n_std standard devitations up or down
			if np.abs(r)<self.n_iqd * r_iqd or v.id == self.outgroup["strain"]:
				new_viruses.append(v)
			else:
				if self.verbose>1:
					print "\t\tresidual:", r, "\nremoved ",v.strain
		self.viruses = MultipleSeqAlignment(new_viruses)

	def clean_generic(self):
		print "Number of viruses before cleaning:",len(self.viruses)
		self.unique_date()
		self.remove_insertions()
		self.clean_ambiguous()
		self.clean_distances()
		self.viruses.sort(key=lambda x:x.num_date)
		print "Number of viruses after outlier filtering:",len(self.viruses)
