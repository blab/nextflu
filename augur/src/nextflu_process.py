import time, argparse,os,subprocess, shutil, glob, sys
sys.path.append('./src')
from nextflu_config import config
from Bio import SeqIO
from io_util import write_json, read_json, write_fasta, read_fasta
from tree_util import dendropy_to_json, json_to_dendropy, delimit_newick
import dendropy

class nextflu(object):
	def __init__(self):
		self.viruses = None
		self.tree = None
		self.initial_virus_fname = 'data/virus_ingest.json'
		self.clean_virus_fname = 'data/virus_clean.json'

	def load_from_file(self, tree_fname=None, virus_fname = None):
		if tree_fname is None: tree_fname = 'data/tree_ancestral.json'
		if os.path.isfile(tree_fname):
			self.tree = dendropy_to_json(read_json(tree_fname))
		if virus_fname is None: virus_fname = 'data/virus_clean.json'
		if os.path.isfile(virus_fname):
			self.viruses = dendropy_to_json(read_json(virus_fname))

	def load_viruses(self, aln_fname = None, years_back=3, viruses_per_month=50):
		if config['virus']:
			from H3N2_filter import H3N2_filter as virus_filter
			fasta_fields = config['fasta_fields']
			force_include_strains = [seq.name for seq in SeqIO.parse('data/strains_with_HI.fasta', 'fasta')]
		else:
			from virus_filter import virus_filter as virus_filter
			fasta_fields = {0:'strain'}
		if aln_fname is None: aln_fname = config['alignment_file']

		my_filter = virus_filter(aln_fname, fasta_fields)
		my_filter.filter()
		my_filter.subsample(years_back, viruses_per_month, prioritize = force_include_strains, 
								all_priority = True, region_specific=False)

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
		tree_refine.main(self.tree, self.viruses)
		write_json(dendropy_to_json(self.tree.seed_node), 'data/tree_refine.json')

	def export_to_auspice(self):
		import streamline
		tree_json = dendropy_to_json(self.tree.seed_node)
		streamline.main(tree_json)

	def run(self,years_back=3, viruses_per_month=50, raxml_time_limit = 1.0):
		self.load_viruses(years_back=years_back, viruses_per_month=viruses_per_month)
		self.clean_viruses()
		self.align()
		self.infer_tree(raxml_time_limit = raxml_time_limit)
		self.infer_ancestral()
		self.refine_tree()
		self.export_to_auspice()

if __name__=="__main__":
	parser = argparse.ArgumentParser(description='Process virus sequences, build tree, and prepare of web visualization')
	parser.add_argument('-y', '--years_back', type = int, default=3, help='number of past years to sample sequences from')
	parser.add_argument('-v', '--viruses_per_month', type = int, default = 50, help='number of viruses sampled per month')
	parser.add_argument('-r', '--raxml_time_limit', type = float, default = 1.0, help='number of hours raxml is run')
	params = parser.parse_args()

	my_nextflu = nextflu()
	my_nextflu.run(**params.__dict__)
