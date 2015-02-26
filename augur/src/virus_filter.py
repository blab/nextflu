# filter viruses after ingest, criteria based on metadata
#  - viruses equal to or longer than 1701 bases
#  - viruses with exact dates
#  - viruses that are not egg-passaged
#  - a single sequence per virus strain, taken as first sequence in list
# outputs to virus_filter.json

import os, re, time, datetime, csv, sys
from io_util import *
from collections import defaultdict
sys.path.append('../source-data')

class virus_filter(object):

	def __init__(self,viruses=None):
		if viruses is None: viruses=[]
		self.viruses = viruses
		self.strain_lookup = {}
		self.outgroup = None

	def filter_generic(self, min_length=None, date_spec = 'full', prepend_strains = None):
		'''
		filter viruses by length and accurate date, sort, add additioanl strains such
		as vaccine strains that are preferentially retained and prune to unique strains
		'''
		print len(self.viruses), "initial viruses"
		if min_length is not None:
			self.filter_length(min_length)
			print len(self.viruses), "after filtering by length >=", min_length

		self.filter_date(date_spec)
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
			label = v['strain'].lower() 
			if not label in self.strain_lookup:
				filtered_viruses.append(v)
				self.strain_lookup[label]=v
		self.viruses=filtered_viruses

	def filter_length(self, min_length):
		self.viruses = filter(lambda v: len(v['seq']) >= min_length, self.viruses)

	def filter_date(self, date_spec):
		if date_spec=='full':
			self.viruses = filter(lambda v: re.match(r'\d\d\d\d-\d\d-\d\d', v['date']) != None, self.viruses)
		elif date_spec=='year':
			self.viruses = filter(lambda v: re.match(r'\d\d\d\d', v['date']) != None, self.viruses)

	def subsample(self, years_back, viruses_per_month, prioritize = None, all_priority=False, region_specific = True):
		'''
		Subsample x viruses per month
		Take from beginning of list - this will prefer longer sequences
		Take viruses 1 per region in a cycle to get geographic diversity
		But pad with additional viruses from populous regions if necessary
		'''
		if prioritize is None:
			prioritize=[]
		else:
			prioritize = [v.lower() for v in prioritize]
		if region_specific:
			select_func = self.select_viruses
		else:
			select_func = self.select_viruses_global

		priority_viruses = self.viruses_by_date_region([v for v in self.viruses if v['strain'].lower() in prioritize]) 
		other_viruses = self.viruses_by_date_region([v for v in self.viruses if v['strain'].lower() not in prioritize]) 

		filtered_viruses = []
		first_year = datetime.datetime.today().year - years_back
		first_month = datetime.datetime.today().month
		regions = list(set([v['region'] for v in self.viruses]))

		print "Filtering between " + str(first_month) + "/" + str(first_year), "and today"
		print "Selecting " + str(viruses_per_month), "viruses per month"
		y = first_year
		for m in range(first_month,13):
			filtered_viruses.extend(select_func(priority_viruses,other_viruses, 
												y, m, viruses_per_month, regions, all_priority=all_priority))
		for y in range(first_year+1,datetime.datetime.today().year+1):
			for m in range(1,13):
				filtered_viruses.extend(select_func(priority_viruses,other_viruses, 
												y, m, viruses_per_month, regions, all_priority=all_priority))
		if self.outgroup is not None:
			filtered_viruses.append(self.outgroup)
			print len(filtered_viruses), "with outgroup"
		self.virus_subsample = filtered_viruses

	def viruses_by_date_region(self, tmp_viruses):
		'''
		make dictionary of lists of viruses belonging to a certain date and region
		'''
		from collections import defaultdict
		virus_tuples = defaultdict(list)
		for v in tmp_viruses:
			vdate = datetime.datetime.strptime(v['date'], '%Y-%m-%d').date()
			virus_tuples[(vdate.year, vdate.month, v['region'])].append(v)

		return virus_tuples

	def select_viruses(self, priority_viruses,other_viruses, y, m, viruses_per_month, regions, all_priority = False):
		'''
		select viruses_per_month strains as evenly as possible from all regions
		'''
		from itertools import izip_longest
		select_set = []
		for vset in [priority_viruses, other_viruses]:
			select_set.append([])
			for representative in izip_longest(*[vset[(y,m,r)] for r in regions], fillvalue = None):
				select_set[-1].extend([v for v in representative if v is not None])
			print "found",len(select_set[-1]), 'in year',y,'month',m
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
		priority_viruses_flat = sum(priority_viruses[(y,m,r)] for r in regions)
		other_viruses_flat = sum(other_viruses[(y,m,r)] for r in regions)
		n_other = max(0,viruses_per_month-len(priority_viruses))
		return sample(priority_viruses[:viruses_per_month], min(len(priority_viruses), viruses_per_month)\
				+ sample(other_viruses, min(n_other, len(other_viruses)))


class flu_filter(virus_filter):

	def __init__(self,fasta_fname, fasta_header=None):
		if fasta_header is None:
			self.fasta_header = {0:'strain', 1:'accession', 3:'passage', 5:'date' }
		else:
			self.fasta_header = fasta_header
		viruses = self.parse_gisaid(fasta_fname)
		virus_filter.__init__(self, viruses)
		self.fix_strain_names()
		self.vaccine_strains=[]

	def parse_gisaid(self, fasta):
		"""Parse FASTA file from GISAID with default header formating"""
		viruses = []
		try:
			handle = open(fasta, 'r')
		except IOError:
			print fasta, "not found"
		else:
			for record in SeqIO.parse(handle, "fasta"):
				words = record.description.replace(">","").replace(" ","").split('|')
				v = {key:words[ii] for ii, key in self.fasta_header.iteritems()}
				v['db']="GISAID"
				v['seq']=str(record.seq).upper()
				if 'passage' not in v: v['passage']=''
				viruses.append(v)
			handle.close()
		return viruses

	def filter(self):
		self.filter_generic(prepend_strains = self.vaccine_strains)	
		self.filter_strain_names()
		print len(self.viruses), "with proper strain names"
		self.filter_passage()
		print len(self.viruses), "without egg passage"
		self.filter_geo()
		print len(self.viruses), "with geographic information"

	def filter_strain_names(self):
		self.viruses = filter(lambda v: re.match(r'^A/', v['strain']) != None, self.viruses)

	def fix_strain_names(self):
		for v in self.viruses:
			v['strain'] = v['strain'].replace('\'','').replace('(','').replace(')','').replace('H3N2','').replace('Human','').replace('human','').replace('//','/')

	def filter_passage(self):
		self.viruses = filter(lambda v: re.match(r'^E\d+', v.get('passage',''), re.I) == None, self.viruses)
		self.viruses = filter(lambda v: re.match(r'^Egg', v.get('passage',''), re.I) == None, self.viruses)

	def filter_geo(self):
		"""Label viruses with geographic location based on strain name"""
		"""Location is to the level of country of administrative division when available"""
		reader = csv.DictReader(open("../source-data/geo_synonyms.tsv"), delimiter='\t')		# list of dicts
		label_to_country = {}
		for line in reader:
			label_to_country[line['label'].lower()] = line['country']
		for v in self.viruses:
			v['country'] = 'Unknown'
			try:
				label = re.match(r'^A/([^/]+)/', v['strain']).group(1).lower()	# check first for whole geo match
				if label in label_to_country:
					v['country'] = label_to_country[label]
				else:
						label = re.match(r'^A/([^\-^\/]+)[\-\/]', v['strain']).group(1).lower()		# check for partial geo match
						if label in label_to_country:
							v['country'] = label_to_country[label]
			except:
				print "couldn't parse", v['strain']

		reader = csv.DictReader(open("../source-data/geo_regions.tsv"), delimiter='\t')		# list of dicts
		country_to_region = {}
		for line in reader:
			country_to_region[line['country']] = line['region']
		for v in self.viruses:
			v['region'] = 'Unknown'
			if v['country'] in country_to_region:
				v['region'] = country_to_region[v['country']]
		
		self.viruses = filter(lambda v: v['region'] != 'Unknown', self.viruses)

