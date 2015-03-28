import time, argparse,re,os
from virus_filter import flu_filter
from virus_clean import virus_clean
from tree_refine import tree_refine
from tree_titer import HI_tree
from fitness_model import fitness_model
from process import process
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Align import MultipleSeqAlignment
import numpy as np
from itertools import izip

epitope_mask = np.fromstring("0000000000000000000000000000000000000000000011111011011001010011000100000001001011110011100110101000001100000100000001000110101011111101011010111110001010011111000101011011111111010010001111101110111001010001110011111111000000111110000000101010101110000000000011100100000001011011100000000000001001011000110111111000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000", dtype='S1')
	
virus_config = {
	# data source and sequence parsing/cleaning/processing
	'virus':'H1N1pdm',
	'alignment_file':'data/H1N1pdm_gisaid_epiflu_sequence.fasta',
	'fasta_fields':{0:'strain', 1:'accession', 3:'passage', 5:'date' },
	'outgroup':'A/Swine/Indiana/P12439/00',
	'force_include':'source-data/H1N1pdm_HI_strains.txt',
	'force_include_all':True,
	'date_spec':'year',
	'max_global':True,   # sample as evenly as possible from different geographic regions 
	'cds':[30,None], # define the HA1 start i n 0 numbering #CHECK
	'n_iqd':3,     # standard deviations from clock

	# frequency estimation parameters
	'aggregate_regions': [  ("global", None), ("NA", ["NorthAmerica"]), ("EU", ["Europe"]), 
							("AS", ["China", "SoutheastAsia", "JapanKorea"]), ("OC", ["Oceania"]) ],
	'frequency_stiffness':10.0,
	'time_interval':(2009.0, 2015.3),
	'pivots_per_year':3.0,
	'min_freq':0.10,
	# define relevant clades in canonical HA1 numbering (+1)
	'clade_designations': {},
	'verbose':2, 
	'tol':1e-4, #tolerance for frequency optimization
	'pc':1e-3, #pseudocount for frequencies 
	'extra_pivots': 6,  # number of pivot point for or after the last observations of a mutations
	'inertia':0.7,		# fraction of frequency change carry over in the stiffness term
	'auspice_frequency_fname':'../auspice/data/H1N1pdm_frequencies.json',
	'auspice_sequences_fname':'../auspice/data/H1N1pdm_sequences.json',
	'auspice_tree_fname':'../auspice/data/H1N1pdm_tree.json',
	'HI_fname':'source-data/H1N1pdm_HI_titers.txt',
}


class H1N1pdm_filter(flu_filter):
	def __init__(self,min_length = 987, **kwargs):
		'''
		parameters
		min_length  -- minimal length for a sequence to be acceptable
		'''
		flu_filter.__init__(self, **kwargs)
		self.min_length = min_length
		self.vaccine_strains =[
			]
		self.outgroup = {
			'strain': 'A/Swine/Indiana/P12439/00',
			'db': 'IRD',
			'accession': 'AF455680',
			'date': '2002-03-14',
			'country': 'USA',
			'region': 'NorthAmerica',
			'seq': 'ATGAAGGCAATACTAGTAGTCCTGCTATATACATTTACAACCGCAAATGCAGACACATTATGTATAGGTTATCATGCGAACAATTCAACTGACACTGTAGACACAGTACTAGAAAAGAATGTAACAGTAACACACTCTGTTAACCTTCTAGAAGACAGGCATAACGGGAAACTATGTAAACTAAGAGGGGTAGCCCCATTGCATTTGGGTAAATGTAACATTGCTGGCTGGCTCCTGGGAAATCCAGAGTGTGAATCACTCTTCACAGCAAGCTCATGGTCCTACATTGTGGAAACATCTAGTTCAGATAATGGGACGTGTTACCCAGGAGATTTCATCAATTATGAAGAGCTAAGAGAGCAATTGAGCTCAGTGTCATCATTTGAAAGATTTGAGATATTCCCCAAGACAAGTTCATGGCCCAATCATGACACGAACAGAGGTGTGACGGCAGCATGTCCTTATGCTGGAGCAAAAAGCTTCTACAGAAATTTAATATGGCTGGTCAAAAAAGAAAATTCATACCCAAAGCTCAGCAAATCCTATATTAACAATAAGGGGAAGGAAGTCCTCGTGCTATGGGGCATTCACCATCCATCTACCAGTGCCGACCAACAAAGTCTCTACCAGAATGCAGATGCATATGTTTTTGTGGGGTCATCAAGATACAGCAAGAAGTTCAAGCCAGAAATAGCAGCCAGACCCAAGGTGAGGGACCAAGCAGGGAGAATAAACTATTACTGGACACTAGTAGAGCCTGGAGACAAAATAACATTCGAAGCAACTGGAAATCTAGTGGTACCGAGATATGCCTTCGCAATGGAAAGAAATTCTGGATCTGGTATTATCATTTCAGATACATCAGTCCACGATTGTAATACGACTTGTCAGACACCCAAGGGTGCTATAAACACCAGCCTCCCATTTCAGAATATACATCCAGTCACAATTGGAGAATGTCCAAAATATGTAAAAAGCACAAAATTGAGAATGGCCACAGGATTAAGGAATGTCCCGTCTATTCAATCTAGAGGCCTGTTTGGGGCCATTGCCGGCTTTATTGAGGGGGGATGGACAGGAATGATAGATGGATGGTACGGTTATCACCATCAAAATGAGCAGGGATCAGGATATGCAGCCGACCTGAAGAGCACACAGAATGCCATTGACGGGATCACTAACAAAGTAAATTCTGTTATTGAAAAGATGAACACACAATTCATAGCAGTAGGTAAAGAGTTCAACCACCTGGAAAAAAGAATAGAGAATTTAAATAAAAAGGTTGATGATGGTTTTCTGGATATTTGGACTTACAATGCCGAACTGTTGATTCTGTTGGAAAATGAAAGAACTTTGGATTACCACGATTCAAATGTGAAGAACTTATATGAAAAGGTAAGAAGCCAGCTAAAAAACAATGCCAGGGAAATTGGGAATGGCTGCTTTGAATTTTACCACAAATGTGATGACAAGTGCATGGAAAGCGTCAAAAATGGGACTTATGATTACCCAAAATACTCAGAGGAAGCAAAACTAAACAGAGAGGAAATAGATGGGGTAAAGCTGGAATCAACAAGGATTTACCAGATTTTGGCGATCTATTCAACTGTCGCCAGTTCATTGGTACTGGTAGTCTCCCTGGGGGCAATCAGTTTCTGGATGTGCTCTAATGGGTCTCTACAGTGTAGAATATGTATTTAA'
		}


class H1N1pdm_clean(virus_clean):
	def __init__(self,**kwargs):
		virus_clean.__init__(self, **kwargs)

	def clean_outbreaks(self):
		"""Remove duplicate strains, where the geographic location, date of sampling and sequence are identical"""
		virus_hashes = set()
		new_viruses = []
		for v in self.viruses:
			geo = re.search(r'A/([^/]+)/', v.strain).group(1)
			if geo:
				vhash = (geo, v.date, str(v.seq))
				if vhash not in virus_hashes:
					new_viruses.append(v)
					virus_hashes.add(vhash)

		self.viruses = MultipleSeqAlignment(new_viruses)
		return new_viruses

	def clean(self):
		self.clean_generic()
		self.clean_outbreaks()
		print "Number of viruses after outbreak filtering:",len(self.viruses)


class H1N1pdm_refine(tree_refine):
	def __init__(self, **kwargs):
		tree_refine.__init__(self, **kwargs)

	def refine(self):
		self.refine_generic()  # -> all nodes now have aa_seq, xvalue, yvalue, trunk, and basic virus properties
		self.add_H1N1pdm_attributes()

	def add_H1N1pdm_attributes(self):
		for v in self.viruses:
			if v.strain in self.node_lookup:
				node = self.node_lookup[v.strain]
				try:
					node.passage=v.passage
				except:
					pass

class H1N1pdm_process(process, H1N1pdm_filter, H1N1pdm_clean, H1N1pdm_refine, HI_tree, fitness_model):
	"""docstring for H1N1pdm_process, H1N1pdm_filter"""
	def __init__(self,verbose = 0, force_include = None, 
				force_include_all = False, max_global= True, **kwargs):
		self.force_include = force_include
		self.force_include_all = force_include_all
		self.max_global = max_global
		process.__init__(self, **kwargs)
		H1N1pdm_filter.__init__(self,**kwargs)
		H1N1pdm_clean.__init__(self,**kwargs)
		H1N1pdm_refine.__init__(self,**kwargs)
		HI_tree.__init__(self,**kwargs)
		fitness_model.__init__(self,**kwargs)
		self.verbose = verbose

	def run(self, steps, years_back=3, viruses_per_month=50, raxml_time_limit = 1.0, reg=2):
		if 'filter' in steps:
			print "--- Virus filtering at " + time.strftime("%H:%M:%S") + " ---"
			self.filter()
			if self.force_include is not None and os.path.isfile(self.force_include):
				with open(self.force_include) as infile:
					forced_strains = [line.strip().upper() for line in infile]
			else:
				forced_strains = []
			self.subsample(years_back, viruses_per_month, 
				prioritize=forced_strains, all_priority=self.force_include_all, 
				region_specific = self.max_global)
			self.dump()
		else:
			self.load()
		if 'align' in steps:
			self.align()   	# -> self.viruses is an alignment object
		if 'clean' in steps:
			print "--- Clean at " + time.strftime("%H:%M:%S") + " ---"
			self.clean()   # -> every node as a numerical date
			self.dump()
		if 'tree' in steps:
			print "--- Tree	 infer at " + time.strftime("%H:%M:%S") + " ---"
			self.infer_tree(raxml_time_limit)  # -> self has a tree
			self.dump()
		if 'ancestral' in steps:
			print "--- Infer ancestral sequences " + time.strftime("%H:%M:%S") + " ---"
			self.infer_ancestral()  # -> every node has a sequence
		if 'refine' in steps:
			print "--- Tree refine at " + time.strftime("%H:%M:%S") + " ---"
			self.refine()
			self.dump()
		if 'frequencies' in steps:
			print "--- Estimating frequencies at " + time.strftime("%H:%M:%S") + " ---"
			self.determine_variable_positions()
			self.estimate_frequencies(tasks = ['mutations', 'tree'])
			self.dump()
		if 'HI' in steps:
			print "--- Adding HI titers to the tree " + time.strftime("%H:%M:%S") + " ---"
			self.map_HI_to_tree(training_fraction=0.9, method = 'nnl1reg', lam_HI=reg, lam_avi=reg, lam_pot = reg)
			self.add_titers()
			self.dump()
		if 'export' in steps:
			self.temporal_regional_statistics()
			# exporting to json, including the H1N1pdm specific fields
			self.export_to_auspice(tree_fields = ['ep', 'ne', 'rb', 'dHI', 'cHI', 'HI_titers', 'serum', 
				'HI_info', 'avidity', 'potency', 'aa_muts'], annotations = ['3c3.a', '3c2.a'])

if __name__=="__main__":
	all_steps = ['filter', 'align', 'clean', 'tree', 'ancestral', 'refine', 'frequencies', 'HI', 'export']
	parser = argparse.ArgumentParser(description='Process virus sequences, build tree, and prepare of web visualization')
	parser.add_argument('-y', '--years_back', type = int, default=3, help='number of past years to sample sequences from')
	parser.add_argument('-v', '--viruses_per_month', type = int, default = 50, help='number of viruses sampled per month')
	parser.add_argument('-r', '--raxml_time_limit', type = float, default = 1.0, help='number of hours raxml is run')
	parser.add_argument('--prefix', type = str, default = 'data/H1N1pdm_', help='path+prefix of file dumps')
	parser.add_argument('--test', default = False, action="store_true",  help ="don't run the pipeline")
	parser.add_argument('--start', default = 'filter', type = str,  help ="start pipeline at virus selection")
	parser.add_argument('--stop', default = 'export', type=str,  help ="run to end")
	params = parser.parse_args()

	steps = all_steps[all_steps.index(params.start):(all_steps.index(params.stop)+1)]
	# add all arguments to virus_config (possibly overriding)
	virus_config.update(params.__dict__)
	# pass all these arguments to the processor: will be passed down as kwargs through all classes
	myH1N1pdm = H1N1pdm_process(**virus_config) 
	if params.test:
		myH1N1pdm.load()
	else:
		myH1N1pdm.run(steps, years_back=virus_config['years_back'], 
			viruses_per_month = virus_config['viruses_per_month'], 
			raxml_time_limit = virus_config['raxml_time_limit'])
