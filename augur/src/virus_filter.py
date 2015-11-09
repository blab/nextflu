# filter viruses after ingest, criteria based on metadata
#  - viruses equal to or longer than 1701 bases
#  - viruses with exact dates
#  - viruses that are not egg-passaged
#  - a single sequence per virus strain, taken as first sequence in list

import os, re, time, datetime, csv, sys, gzip
from collections import defaultdict
from Bio import SeqIO
import numpy as np

def myopen(fname, mode='r'):
	if fname[-2:] == 'gz':
		return gzip.open(fname, mode)
	else:
		return open(fname, mode)

def fix_name(name):
	tmp_name = name.replace(' ', '').replace('\'','').replace('(','').replace(')','').replace('H3N2','').replace('Human','').replace('human','').replace('//','/')
	fields = tmp_name.split('/')
	if len(fields[-1])==2:
		try:
			y = int(fields[-1])
			if y>16:
				y=1900+y
			else:
				y=2000+y
			return '/'.join(fields[:-1])+'/'+str(y)
		except:
			return tmp_name
	else:
		return tmp_name

class virus_filter(object):

	def __init__(self, alignment_file='', fasta_fields=None, date_spec='full', **kwargs):
		'''
		parameters:
		alignment_file   -- a FASTA sequence file with all viruses included
		fasta_fields     -- match between FASTA header index and virus dict field
		date_spec        -- if 'full', dates with day are required, if 'year', only year is accepted
		'''
		if fasta_fields is None:
			self.fasta_fields = {0:'strain', 1:'date' }
		else:
			self.fasta_fields = fasta_fields
		self.alignment_file = alignment_file
		self.viruses = self.parse_fasta(self.alignment_file)
		self.strain_lookup = {}
		self.outgroup = None
		self.date_spec = date_spec
		
	def parse_fasta(self, fasta):
		"""Parse FASTA file with default header formating"""
		viruses = []
		try:
			handle = myopen(fasta, 'r')
		except IOError:
			print fasta, "not found"
		else:
			for record in SeqIO.parse(handle, "fasta"):
				words = map(lambda x:x.strip(),record.description.replace(">","").split('|'))
				v = {key: words[ii] if ii<len(words) else "" for ii, key in self.fasta_fields.iteritems()}
				v['seq']= str(record.seq)
				viruses.append(v)
			handle.close()
		return viruses		
		
	def filter(self):
		self.filter_generic()			

	def filter_generic(self, prepend_strains = None):
		'''
		filter viruses by length and accurate date, sort, add additioanl strains such
		as vaccine strains that are preferentially retained and prune to unique strains
		'''
		print len(self.viruses), "initial viruses"
		if hasattr(self, 'min_length'):
			self.filter_length(self.min_length)
			print len(self.viruses), "after filtering by length >=", self.min_length
		self.filter_noncanoncial_nucleotides()
		print len(self.viruses), "after filtering bad nucleotides"

		self.filter_date()
		print len(self.viruses), "after filtering for precise dates"
		self.sort_length()
		if prepend_strains is not None:
			self.viruses = prepend_strains + self.viruses
			print len(self.viruses), "after adding custom strains"
		self.filter_unique()
		print len(self.viruses), "after filtering for unique strains"

	def sort_length(self):
		'''	
		Sort by length, but randomize viruses of a given length
		'''	
		from random import shuffle
		shuffle(self.viruses)
		self.viruses.sort(key = lambda v: len(v['seq']), reverse = True)

	def filter_unique(self):
		'''		
		Keep only the first isolate of a strain
		'''	
		filtered_viruses = []
		for v in self.viruses:
			label = v['strain'].upper() 
			if not label in self.strain_lookup:
				filtered_viruses.append(v)
				self.strain_lookup[label]=v
		self.viruses=filtered_viruses

	def filter_length(self, min_length):
		self.viruses = filter(lambda v: len(v['seq']) >= min_length, self.viruses)

	def filter_noncanoncial_nucleotides(self, max_bad_pos=3):
		self.viruses = filter(lambda v: sum(v['seq'].count(nuc) for nuc in 'ACGTacgt') >= len(v['seq'])-max_bad_pos, self.viruses)

	def filter_date(self):
		if self.date_spec=='full':
			self.viruses = filter(lambda v: re.match(self.date_format['reg'], v['date']) is not None, self.viruses)
		elif self.date_spec=='year':
			self.viruses = filter(lambda v: re.match(r'\d\d\d\d', v['date']) is not None, self.viruses)
			for v in self.viruses:
				if re.match(r'\d\d\d\d-\d\d-\d\d', v['date']) is None:
					v['date'] = v['date'][:4]+'-'+format(np.random.randint(12)+1, '02d')+'-01' 

	def subsample(self, viruses_per_month, prioritize = None, all_priority=False, region_specific = True):
		'''
		Subsample x viruses per month
		Take from beginning of list - this will prefer longer sequences
		Take viruses 1 per region in a cycle to get geographic diversity
		But pad with additional viruses from populous regions if necessary
		'''
		if prioritize is None:
			prioritize=[]
		else:
			prioritize = [v.upper() for v in prioritize]
		if region_specific:
			select_func = self.select_viruses
		else:
			select_func = self.select_viruses_global

		priority_viruses = self.viruses_by_date_region([v for v in self.viruses if v['strain'].upper() in prioritize]) 
		other_viruses = self.viruses_by_date_region([v for v in self.viruses if v['strain'].upper() not in prioritize]) 

		filtered_viruses = []
		first_year = int(np.floor(self.time_interval[0]))
		first_month = int((self.time_interval[0]-first_year)*12)
		regions = list(set([v['region'] for v in self.viruses]))

		print "Filtering between " + str(first_month) + "/" + str(first_year), "and today"
		print "Selecting " + str(viruses_per_month), "viruses per month"
		y = first_year
		for m in range(first_month,13):
			filtered_viruses.extend(select_func(priority_viruses,other_viruses, 
												y, m, viruses_per_month, regions, all_priority=all_priority))
		for y in range(first_year+1,int(np.floor(self.time_interval[1]))+1):
			for m in range(1,13):
				filtered_viruses.extend(select_func(priority_viruses,other_viruses, 
												y, m, viruses_per_month, regions, all_priority=all_priority))
				if y+float(m)/12.0>self.time_interval[1]:
					break
		if self.outgroup is not None:
			filtered_viruses.append(self.outgroup)
			print len(filtered_viruses), "with outgroup"
		self.viruses = filtered_viruses

	def viruses_by_date_region(self, tmp_viruses):
		'''
		make dictionary of lists of viruses belonging to a certain date and region
		'''
		from collections import defaultdict
		virus_tuples = defaultdict(list)
		for v in tmp_viruses:
			try:
				vdate = datetime.datetime.strptime(v['date'], self.date_format['fields']).date()
			except:
				print "incomplete date!", v['strain'], v['date'], "adjusting to July 1st"
				v['date']+='-07-01'	
				vdate = datetime.datetime.strptime(v['date'], '%Y-%m-%d').date()
			virus_tuples[(vdate.year, vdate.month, v['region'])].append(v)

		return virus_tuples

	def select_viruses(self, priority_viruses,other_viruses, y, m, viruses_per_month, regions, all_priority = False):
		'''
		select viruses_per_month strains as evenly as possible from all regions
		'''
		from itertools import izip_longest
		from random import sample,shuffle
		select_set = []
		for vset in [priority_viruses, other_viruses]:
			select_set.append([])
			for representative in izip_longest(*[vset[(y,m,r)] for r in regions], fillvalue = None):
				tmp = [v for v in representative if v is not None]
				shuffle(tmp)
				select_set[-1].extend(tmp)
		if self.verbose>1:
			print "\tfound",len(select_set[-1]), 'in year',y,'month',m
		if all_priority:
			n_other = max(0,viruses_per_month-len(select_set[0]))
			return select_set[0] + select_set[1][:n_other]
		else:
			tmp = select_set[0] + select_set[1]
			return tmp[:viruses_per_month]

	def select_viruses_global(self, priority_viruses,other_viruses, y, m, viruses_per_month, regions, all_priority = False):
		'''
		select viruses_per_month strains as evenly as possible from all regions
		'''
		from random import sample
		priority_viruses_flat = []
		for r in regions: priority_viruses_flat.extend(priority_viruses[(y,m,r)])
		other_viruses_flat = []
		for r in regions: other_viruses_flat.extend(other_viruses[(y,m,r)])

		if self.verbose>1:
			print "\t\tfound",len(priority_viruses_flat)+len(other_viruses_flat), 'in year',y,'month',m
		n_other = max(0,viruses_per_month-len(priority_viruses_flat))
		return sample(priority_viruses_flat, len(priority_viruses_flat) if all_priority else min(len(priority_viruses_flat), viruses_per_month))\
				+ sample(other_viruses_flat, min(n_other, len(other_viruses_flat)))


class flu_filter(virus_filter):

	def __init__(self, alignment_file='', fasta_fields=None, **kwargs):	
		virus_filter.__init__(self, alignment_file = alignment_file, fasta_fields = fasta_fields, **kwargs)
		self.add_gisaid_metadata()
		self.fix_strain_names()
		self.vaccine_strains=[]

	def filter(self):
		self.filter_strain_names()
		print len(self.viruses), "with proper strain names"
		self.filter_passage()
		print len(self.viruses), "without egg passage"
		self.filter_generic(prepend_strains = self.vaccine_strains)	
		self.filter_geo(prune=False)
		print len(self.viruses), "with geographic information"
		
	def add_gisaid_metadata(self):
		for v in self.viruses:
			v['db']="GISAID"

	def filter_strain_names(self):
		self.viruses = filter(lambda v: re.match(r'^[AB]/', v['strain']) != None, self.viruses)

	def fix_strain_names(self):
		for v in self.viruses:
			v['strain'] = fix_name(v['strain'])

	def filter_passage(self):
		self.viruses = filter(lambda v: re.match(r'^E\d+', v.get('passage',''), re.I) == None, self.viruses)
		self.viruses = filter(lambda v: re.match(r'^Egg', v.get('passage',''), re.I) == None, self.viruses)

	def filter_geo(self, prune = True):
		"""Label viruses with geographic location based on strain name"""
		"""Location is to the level of country of administrative division when available"""
		reader = csv.DictReader(open("source-data/geo_synonyms.tsv"), delimiter='\t')		# list of dicts
		label_to_country = {}
		for line in reader:
			label_to_country[line['label'].lower()] = line['country']
		for v in self.viruses:
			if "country" not in v:
				v['country'] = 'Unknown'
				try:
					label = re.match(r'^[AB]/([^/]+)/', v['strain']).group(1).lower()						# check first for whole geo match
					if label in label_to_country:
						v['country'] = label_to_country[label]
					else:
						label = re.match(r'^[AB]/([^\-^\/]+)[\-\/]', v['strain']).group(1).lower()			# check for partial geo match
					if label in label_to_country:
						v['country'] = label_to_country[label]
					else:
						label = re.match(r'^[AB]/([A-Z][a-z]+)[A-Z0-9]', v['strain']).group(1).lower()			# check for partial geo match
					if label in label_to_country:
						v['country'] = label_to_country[label]							
					if v['country'] == 'Unknown':
						print "couldn't parse location for", v['strain']
				except:
					print "couldn't parse location for", v['strain']

		reader = csv.DictReader(open("source-data/geo_regions.tsv"), delimiter='\t')		# list of dicts
		country_to_region = {}
		for line in reader:
			country_to_region[line['country']] = line['region']
		for v in self.viruses:
			v['region'] = 'Unknown'
			if v['country'] in country_to_region:
				v['region'] = country_to_region[v['country']]
			if v['country'] != 'Unknown' and v['region'] == 'Unknown':
				print "couldn't parse region for", v['strain'], "country:", v["country"]		
		
		if prune:
			self.viruses = filter(lambda v: v['region'] != 'Unknown', self.viruses)
