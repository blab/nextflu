import time, argparse,os,subprocess, shutil, glob
from nextflu_config import config
from Bio import SeqIO
from io_util import write_json, read_json, write_fasta, read_fasta
from tree_util import dendropy_to_json, json_to_dendropy
import dendropy

class nextflu(object):
	def __init__(self):
		self.viruses = None
		self.tree = None
		self.initial_virus_fname = 'data/virus_ingest.json'
		self.clean_virus_fname = 'data/virus_clean.json'

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

	def clean_viruses(self, virus_fname=None):
		import virus_clean
		print "--- Clean at " + time.strftime("%H:%M:%S") + " ---"
		if virus_fname is None:
			if self.viruses is None:
				self.viruses = read_json(self.initial_virus_fname)
		else:
			self.viruses = read_json(virus_fname)

		print str(len(self.viruses)) + " initial self.viruses"
		# mask extraneous columns and ambiguous bases
		virus_clean.mask_from_outgroup(self.viruses)
		virus_clean.clean_ambiguous(self.viruses)

		# clean gapped sequences
	#	self.viruses = clean_gaps(self.viruses)
	#	print str(len(self.viruses)) + " with complete HA"

		# clean sequences by outbreak
		self.viruses = virus_clean.clean_outbreaks(self.viruses)
		print str(len(self.viruses)) + " with outbreak sequences removed"

		# clean reassortant sequences
		self.viruses = virus_clean.clean_reassortants(self.viruses)
		print str(len(self.viruses)) + " with triple reassortants removed"

		# clean sequences by distance
		self.viruses = virus_clean.clean_distances(self.viruses)
		print str(len(self.viruses)) + " with clock"

		write_json(self.viruses, self.clean_virus_fname)

	def align(self):
		import virus_align
		write_fasta(self.viruses, 'temp_in.fasta')
		os.system("mafft --nofft temp_in.fasta > temp_out.fasta")
		alignment = read_fasta('temp_out.fasta')
		virus_align.update_viruses(alignment, self.viruses)
		out_fname = 'data/virus_align.json'
		write_json(self.viruses, out_fname)
		virus_align.cleanup()

	def infer_tree(self, virus_fname = None, raxml_time_limit = 1.0):
		print "--- Tree infer at " + time.strftime("%H:%M:%S") + " ---"
		import tree_infer
		if virus_fname is not None:
			self.viruses = read_json(virus_fname)
		else:
			if self.viruses is None:
				self.viruses = read_json(self.clean_virus_fname)

		tree_infer.cleanup()
		write_fasta(self.viruses, 'temp.fasta')
		print "Building initial tree with FastTree"
		os.system("fasttree -gtr -nt -gamma -nosupport -mlacc 2 -slownni temp.fasta > initial_tree.newick")
		tree_infer.delimit_newick("initial_tree.newick", "temp.newick")
		self.tree = dendropy.Tree.get_from_path("temp.newick", "newick")
		self.tree.resolve_polytomies()
		self.tree.write_to_path("initial_tree.newick", "newick")

		if raxml_time_limit>0:
			print "RAxML tree optimization with time limit " + str(raxml_time_limit) + " hours"
			os.system("seqmagick convert temp.fasta temp.phyx")
			# using exec to be able to kill process
			end_time = time.time() + int(raxml_time_limit*3600)
			process = subprocess.Popen("exec raxml -f d -T 6 -j -s temp.phyx -n topology -c 25 -m GTRCAT -p 344312987 -t initial_tree.newick", shell=True)
			while (time.time() < end_time):
				if os.path.isfile('raxml_result.topology'):
					break
				time.sleep(10)
			process.terminate()

			checkpoint_files = [file for file in glob.glob("RAxML_checkpoint*")]
			if os.path.isfile('raxml_result.topology'):
				checkpoint_files.append('raxml_result.topology')
			if len(checkpoint_files) > 0:
				last_tree_file = checkpoint_files[-1]
				shutil.copy(last_tree_file, 'raxml_tree.newick')
			else:
				shutil.copy("initial_tree.newick", 'raxml_tree.newick')

			print "RAxML branch length optimization and rooting"
			os.system("raxml -f e -T 6 -s temp.phyx -n branches -c 25 -m GTRGAMMA -p 344312987 -t raxml_tree.newick -o " + config["outgroup"])

			tree_infer.delimit_newick("RAxML_result.branches", "temp.newick")
			self.tree = dendropy.Tree.get_from_path("temp.newick", "newick")
			self.tree.resolve_polytomies()
		tree_infer.cleanup()
		self.tree.write_to_path("data/tree_infer.newick", "newick")

	def infer_ancestral(self, tree_fname="data/tree_infer.newick", virus_fname = None):
		from tree_ancestral import ancestral_sequences
		from seq_util import json_to_Bio_alignment
		from tree_util import BioPhylo_to_json
		print "--- Ancestral inference at " + time.strftime("%H:%M:%S") + " ---"
		if virus_fname is not None:
			self.viruses = read_json(virus_fname)
		else:
			if self.viruses is None:
				self.viruses = read_json(self.clean_virus_fname)
		from Bio import Phylo
		aln = json_to_Bio_alignment(self.viruses)
		biotree = Phylo.read(tree_fname, 'newick')
		print "--- Set-up ancestral inference at " + time.strftime("%H:%M:%S") + " ---"
		anc_seq = ancestral_sequences(biotree, aln, seqtype='str')
		anc_seq.calc_ancestral_sequences()
		anc_seq.cleanup_tree()
		out_fname = "data/tree_ancestral.json"
		write_json(BioPhylo_to_json(anc_seq.T.root), out_fname)
		self.tree = json_to_dendropy(read_json(out_fname))

	def refine_tree(self):
		import tree_refine
		print "--- Tree refine at " + time.strftime("%H:%M:%S") + " ---"
		print "Remove outgroup"
		tree_refine.remove_outgroup(self.tree)
		print "Remove outlier branches"
		tree_refine.reduce(self.tree)
		print "Collapse internal nodes"
		tree_refine.collapse(self.tree)
		print "Ladderize tree"
		tree_refine.ladderize(self.tree)
		print "Append node attributes"
		tree_refine.add_virus_attributes(self.viruses, self.tree)
		tree_refine.add_node_attributes(self.tree)
		print "Translate nucleotide sequences"
		tree_refine.translate_all(self.tree)
		print "Enumerate leaves of ladderized tree and calculate unique numerical date"
		tree_refine.unique_date(self.tree)
		print "Define trunk"
		tree_refine.define_trunk(self.tree)
		out_fname = "data/self.tree_refine.json"
		write_json(dendropy_to_json(self.tree.seed_node), out_fname)
		return out_fname

	def streamline(self):
		from tree_util import all_descendants
		print "--- Streamline at " + time.strftime("%H:%M:%S") + " ---"
		# Move sequence data to separate file
		print "Writing sequences"	
		tree_json = dendropy_to_json(self.tree.seed_node)
		elems = []
		for node in all_descendants(tree_json):
			elem = {}
			if 'clade' in node:
				elem['clade'] = node['clade']
			if 'aa_seq' in node:
				elem['aa_seq'] = node['aa_seq']			
			elems.append(elem)
		write_json(elems, "../auspice/data/sequences.json", indent=None)

		# Streamline tree for auspice
		print "Writing streamlined tree"
		for node in all_descendants(tree_json):
			node.pop("seq", None)
			node.pop("aa_seq", None)
			node.pop("logit_freq", None)

		out_fname_tree = "../auspice/data/tree.json"
		write_json(tree_json, out_fname_tree, indent=None)
		try:
			read_json(out_fname_tree)
		except:
			print "Read failed, rewriting with indents"	
			write_json(self.tree, out_fname_tree, indent=1)
			
		# Include genotype frequencies
		shutil.copy2("data/genotype_frequencies.json", "../auspice/data/frequencies.json")


	def run(self,years_back=3, viruses_per_month=50, raxml_time_limit = 1.0):
		self.load_viruses(years_back=years_back, viruses_per_month=viruses_per_month)
		self.clean_viruses()
		self.align()
		self.infer_tree(raxml_time_limit = raxml_time_limit)
		self.infer_ancestral()
		self.refine_tree()
		self.streamline()

if __name__=="__main__":
	parser = argparse.ArgumentParser(description='Process virus sequences, build tree, and prepare of web visualization')
	parser.add_argument('-y', '--years_back', type = int, default=3, help='number of past years to sample sequences from')
	parser.add_argument('-v', '--viruses_per_month', type = int, default = 50, help='number of viruses sampled per month')
	parser.add_argument('-r', '--raxml_time_limit', type = float, default = 1.0, help='number of hours raxml is run')
	params = parser.parse_args()

	my_nextflu = nextflu()
	my_nextflu.run(**params.__dict__)
