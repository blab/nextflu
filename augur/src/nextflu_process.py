import time, argparse,os,subprocess, shutil, glob, sys
sys.path.append('./src')
from Bio import SeqIO
from io_util import write_json, read_json, write_fasta, read_fasta
from tree_util import dendropy_to_json, json_to_dendropy, delimit_newick
import dendropy

class nextflu(object):
	def __init__(self):
		self.viruses = None
		self.tree = None
		self.frequencies = {}
		self.initial_virus_fname = 'data/virus_ingest.json'
		self.clean_virus_fname = 'data/virus_clean.json'
		self.intermediate_tree_fname = 'data/tree_refine.json'
		self.frequency_fname = 'data/frequencies.json'

	def load_from_file(self, tree_fname=None, virus_fname = None):
		if tree_fname is None: tree_fname = self.intermediate_tree_fname
		if os.path.isfile(tree_fname):
			self.tree = json_to_dendropy(read_json(tree_fname))
		if virus_fname is None: virus_fname = self.clean_virus_fname
		if os.path.isfile(virus_fname):
			self.viruses = read_json(virus_fname)
		if os.path.isfile(self.frequency_fname):
			self.frequencies = read_json(self.frequency_fname)

	def load_viruses(self, aln_fname = None, years_back=3, viruses_per_month=50):
		if config['virus']:
			from H3N2_filter import H3N2_filter as virus_filter
			fasta_fields = config['fasta_fields']
			if 'force_include' in config and os.path.isfile(config['force_include']):
				with open(config['force_include']) as force_include_file:
					force_include_strains = [line.strip() for line in force_include_file]
			else:
				force_include_strains = []
		else:
			from virus_filter import virus_filter as virus_filter
			fasta_fields = {0:'strain'}
		if aln_fname is None: aln_fname = config['alignment_file']

		my_filter = virus_filter(aln_fname, fasta_fields)
		my_filter.filter()
		my_filter.subsample(years_back, viruses_per_month, prioritize = force_include_strains, 
								all_priority = True, region_specific=config['max_global'])

		self.viruses = my_filter.virus_subsample
		write_json(self.viruses, self.initial_virus_fname)

	def clean_viruses(self):
		import virus_clean
		self.viruses = virus_clean.main(self.viruses)
		write_json(self.viruses, self.clean_virus_fname)

	def align(self):
		import virus_align
		self.viruses = virus_align.main(self.viruses)
		out_fname = 'data/virus_align.json'
		write_json(self.viruses, out_fname)

	def infer_tree(self, raxml_time_limit = 1.0):
		import tree_infer
		tree_fname = tree_infer.main(self.viruses, raxml_time_limit, config['outgroup'])
		delimit_newick(tree_fname, "temp.newick")
		self.tree = dendropy.Tree.get_from_path("temp.newick", "newick")
		os.remove('temp.newick')

	def infer_ancestral(self, virus_fname = None):
		import tree_ancestral
		self.tree = tree_ancestral.main(self.tree, self.viruses)

	def refine_tree(self):
		import tree_refine
		tree_refine.main(self.tree, self.viruses, config['outgroup'], config['cds'])
		write_json(dendropy_to_json(self.tree.seed_node), self.intermediate_tree_fname)

	def estimate_frequencies(self, tasks = ['mutations','genotypes', 'clades', 'tree']):
		import bernoulli_frequency as freq_est
		plot=False
		freq_est.flu_stiffness = config['frequency_stiffness']
		freq_est.time_interval = config['time_interval']
		freq_est.pivots_per_year = config['pivots_per_year']
		freq_est.relevant_pos_cutoff = 0.1

		if 'mutations' in tasks or 'genotypes' in tasks:
			self.frequencies['mutations'], relevant_pos = freq_est.all_mutations(self.tree, config['aggregate_regions'], 
														threshold = config['min_mutation_count'], plot=plot)
		if 'genotypes' in tasks:
			self.frequencies['genotypes'] = freq_est.all_genotypes(self.tree, config['aggregate_regions'], relevant_pos)
		if 'clades' in tasks:
			self.frequencies['clades'] = freq_est.all_clades(self.tree, config['clade_designations'], 
															config['aggregate_regions'], plot)
		if any(x in tasks for x in ['mutations','clades', 'genotypes']):
			write_json(self.frequencies, self.frequency_fname)

		if 'tree' in tasks:
			for region_label, regions in config['aggregate_regions']:
				print "--- "+"adding frequencies to tree "+region_label+ " "  + time.strftime("%H:%M:%S") + " ---"
				freq_est.estimate_tree_frequencies(self.tree, threshold = 10, regions=regions, region_name=region_label)

	def export_to_auspice(self):
		import streamline
		tree_json = dendropy_to_json(self.tree.seed_node)
		streamline.main(tree_json, self.frequencies)

	def run(self,years_back=3, viruses_per_month=50, raxml_time_limit = 1.0,  **kwargs):
		self.load_viruses(years_back=years_back, viruses_per_month=viruses_per_month)
		self.align()
		self.clean_viruses()
		self.infer_tree(raxml_time_limit = raxml_time_limit)
		self.infer_ancestral()
		self.refine_tree()
		self.estimate_frequencies()
		self.export_to_auspice()

if __name__=="__main__":
	parser = argparse.ArgumentParser(description='Process virus sequences, build tree, and prepare of web visualization')
	parser.add_argument('-y', '--years_back', type = int, default=3, help='number of past years to sample sequences from')
	parser.add_argument('-v', '--viruses_per_month', type = int, default = 50, help='number of viruses sampled per month')
	parser.add_argument('-r', '--raxml_time_limit', type = float, default = 1.0, help='number of hours raxml is run')
	parser.add_argument('--config', default = "nextflu_config.py" , type=str, help ="config file")
	parser.add_argument('--test', default = False, action="store_true",  help ="don't run the pipeline")
	parser.add_argument('--virus', default = False, action="store_true",  help ="only select viruses")
	parser.add_argument('--tree', default = False, action="store_true",  help ="only build tree")
	parser.add_argument('--frequencies', default = False, action="store_true",  help ="only estimate frequencies")
	params = parser.parse_args()

	execfile(params.config)
	print config

	my_nextflu = nextflu()
	my_nextflu.load_from_file()
	if params.virus:
		my_nextflu.load_viruses(years_back=params.years_back, viruses_per_month = params.viruses_per_month)
		my_nextflu.align()
		my_nextflu.clean_viruses()
	elif params.tree:
		my_nextflu.infer_tree(raxml_time_limit=params.raxml_time_limit)
		my_nextflu.infer_ancestral()
		my_nextflu.refine_tree()
	elif params.frequencies:
		my_nextflu.estimate_frequencies()
	elif not params.test:
		my_nextflu.run(**params.__dict__)
