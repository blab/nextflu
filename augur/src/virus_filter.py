# filter viruses after ingest, criteria based on metadata
#  - viruses equal to or longer than 1701 bases
#  - viruses with exact dates
#  - viruses that are not egg-passaged
#  - a single sequence per virus strain, taken as first sequence in list

import os, re, time, datetime, csv, sys
from collections import defaultdict
from Bio import SeqIO
import numpy as np

def get_geo(label):
	from geopy.geocoders import GoogleV3
	geoobj = GoogleV3()
	loc = geoobj.geocode(label)
	if loc is not None:
		country, location = None, None
		for c in loc.raw['address_components']:
			if 'country' in c['types']: 
				country = c['long_name']
			if 'administrative_area_level_1' in c['types']: 
				location = c['long_name']
		if location is None:
			location=country
		return location, country
	else:
		return None, None

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
			handle = open(fasta, 'r')
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

		self.filter_date()
		print len(self.viruses), "after filtering for precise dates"
		self.sort_length()
		if prepend_strains is not None:
			self.viruses = prepend_strains + self.viruses
			print len(self.viruses), "after adding custom strains"
		self.parse_and_filter_label()
		print len(self.viruses), "after filtering for label"
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

	def parse_and_filter_label(self):
		tmp_viruses = []
		for v in self.viruses:
			strain_info, strain_ok = self.parse_strain_name(v['strain'])
			if strain_ok:
				for key, val in strain_info.iteritems():
					v[key]=val
				tmp_viruses.append(v)
		self.viruses=tmp_viruses

#	def parse_strain_name(self, strain_name):
#		'''
#		meant to me over written by subclass
#		'''
#		return {}, True

	def filter_length(self, min_length):
		self.viruses = filter(lambda v: len(v['seq']) >= min_length, self.viruses)

	def filter_date(self):
		if self.date_spec=='full':
			self.viruses = filter(lambda v: re.match(self.date_format['reg'], v['date']) != None, self.viruses)
		elif self.date_spec=='year':
			self.viruses = filter(lambda v: re.match(r'\d\d\d\d', v['date']) != None, self.viruses)

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
			prioritize = [v.lower() for v in prioritize]
		if region_specific:
			select_func = self.select_viruses
		else:
			select_func = self.select_viruses_global

		priority_viruses = self.viruses_by_date_region([v for v in self.viruses if v['strain'].lower() in prioritize]) 
		other_viruses = self.viruses_by_date_region([v for v in self.viruses if v['strain'].lower() not in prioritize]) 

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

	def __init__(self, alignment_file='', fasta_fields=None, strict_geo=True, 
				strict_host=False, **kwargs):	
		virus_filter.__init__(self, alignment_file = alignment_file, fasta_fields = fasta_fields, **kwargs)
		self.strict_host = strict_host
		self.strict_geo = strict_geo
		self.add_gisaid_metadata()
		self.fix_strain_names()
		self.vaccine_strains=[]
		self.load_strain_name_parsing_info()
		self.fix_geo=False

	def load_strain_name_parsing_info(self):
		from csv import DictReader
		self.label_to_country = {}
		for line in DictReader(open("source-data/geo_synonyms.tsv"), delimiter='\t'):
			self.label_to_country[line['label'].lower()] = line['country']

		self.country_to_region = {}
		for line in DictReader(open("source-data/geo_regions.tsv"), delimiter='\t'):
			self.country_to_region[line['country']] = line['region']

		self.label_to_animal = {}
		for line in DictReader(open("source-data/host_synonyms.tsv"), delimiter='\t'):
			self.label_to_animal[line['label'].lower()] = line['animal']

		self.animal_to_group = {}
		for line in DictReader(open("source-data/host_animal_groups.tsv"), delimiter='\t'):
			self.animal_to_group[line['animal']] = line['group']

	def parse_strain_name(self, strain_name):
		fields = map(lambda x:x.strip().lower(), strain_name.split('/'))
		strain_info = {}
		def add_geo(place):
			country = self.label_to_country[place]
			strain_info['country'] = country
			if country in self.country_to_region:
				strain_info['region'] = self.country_to_region[country]
			else:
				strain_info['region'] = 'Unknown'
		def add_host(host):
			animal = self.label_to_animal[host]
			strain_info['host'] = animal
			if animal in self.animal_to_group:
				strain_info['group'] = self.animal_to_group[animal]
			else:
				strain_info['group'] = 'Unknown'

		if fields[0].upper() not in ['A', 'B']:
			return strain_info, False
		if fields[2] in self.label_to_country:
			add_geo(fields[2])
			if fields[1] in self.label_to_animal:
				add_host(fields[1])
			else:
				strain_info['host'] = 'Unknown'
				strain_info['group'] = 'Unknown'
		elif fields[1] in self.label_to_animal and fields[1]!='turkey':
				add_host(fields[1])
				if self.fix_geo:
					country = self.determine_country(fields[2])
					if country is not None:
						self.label_to_country[fields[2]]=country
						add_geo(fields[2])
					else:
						strain_info['country'] = country
						strain_info['region'] = 'Unknown'
				else:
					strain_info['country'] = country
					strain_info['region'] = 'Unknown'
		elif fields[1] in self.label_to_country:
			add_geo(fields[1])
			add_host("human")
		else:
			if self.fix_geo:
				country = self.determine_country(fields[1])
				if country is not None:
					self.label_to_country[fields[1]]=country
					add_geo(fields[1])
				else:
					strain_info['country'] = country
					strain_info['region'] = 'Unknown'
			else:
				strain_info['country'] = country
				strain_info['region'] = 'Unknown'
			strain_info['host'] = 'Unknown'
			strain_info['group'] = 'Unknown'
		return strain_info, not(((strain_info['country']=='Unknown') and self.strict_geo) or\
								((strain_info['host']=='Unknown') and self.strict_host))


	def determine_country(self, label):
		location, country = get_geo(label)
		if country is not None:
			with open('source-data/geo_synonyms.tsv', 'a') as outf:
				try:
					print "adding:",label+'\t'+country+'\t'+location
					outf.write(label+'\t'+country+'\t'+location+'\n')
				except:
					print "Can't write to file:",label+'\t'+country+'\t'+location
		return country


	def filter(self):
		self.filter_generic(prepend_strains = self.vaccine_strains)	
		self.filter_passage()
		print len(self.viruses), "without egg passage"
		
	def add_gisaid_metadata(self):
		for v in self.viruses:
			v['db']="GISAID"

	def fix_strain_names(self):
		for v in self.viruses:
			v['strain'] = v['strain'].replace(' ', '').replace('\'','').replace('(','').replace(')','').replace('H3N2','').replace('Human','').replace('human','').replace('//','/')

	def filter_passage(self):
		self.viruses = filter(lambda v: re.match(r'^E\d+', v.get('passage',''), re.I) == None, self.viruses)
		self.viruses = filter(lambda v: re.match(r'^Egg', v.get('passage',''), re.I) == None, self.viruses)

