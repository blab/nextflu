import sys, time, os, argparse,shutil,subprocess, glob
sys.path.append('src')
sys.setrecursionlimit(10000)  # needed since we are dealing with large trees
from Bio import SeqIO, AlignIO,Phylo
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import dendropy
from bernoulli_frequency import virus_frequencies
from tree_util import delimit_newick
import numpy as np

parser = argparse.ArgumentParser(description='Process virus sequences, build tree, and prepare of web visualization')
parser.add_argument('-y', '--years_back', type = float, default=3, help='number of past years to sample sequences from')
parser.add_argument('-v', '--viruses_per_month', type = int, default = 50, help='number of viruses sampled per month')
parser.add_argument('-r', '--raxml_time_limit', type = float, default = 1.0, help='number of hours raxml is run')
parser.add_argument('--interval', nargs = '+', type = float, default = None, help='interval from which to pull sequences')
parser.add_argument('--prefix', type = str, default = 'data/', help='path+prefix of file dumps')
parser.add_argument('--test', default = False, action="store_true",  help ="don't run the pipeline")
parser.add_argument('--start', default = 'filter', type = str,  help ="start pipeline at specified step")
parser.add_argument('--stop', default = 'export', type=str,  help ="run to end")
parser.add_argument('--skip', nargs='+', type = str,  help ="analysis steps to skip")	


virus_config = {
	'fasta_fields':{0:'strain', 1:'isolate_id', 3:'passage', 5:'date', 7:'lab', 8:"accession"},
	# frequency estimation parameters
	'aggregate_regions': [  ("global", None), ("NA", ["NorthAmerica"]), ("EU", ["Europe"]), 
							("AS", ["China", "SoutheastAsia", "JapanKorea"]), ("OC", ["Oceania"]) ],
	'frequency_stiffness':5.0,
	'verbose':2, 
	'tol':1e-4, #tolerance for frequency optimization
	'pc':1e-3, #pseudocount for frequencies 
	'extra_pivots': 12,  # number of pivot point for or after the last observations of a mutations
	'inertia':0.7,		# fraction of frequency change carry over in the stiffness term
	'n_iqd':3,     # standard deviations from clock
}



class process(virus_frequencies):
	"""generic template class for processing virus sequences into trees"""
	def __init__(self, prefix = 'data/', time_interval = (2012.0, 2015.0), 
	             auspice_prefix = '', 
				 min_mutation_frequency = 0.01, min_genotype_frequency = 0.1, **kwargs):
		self.tree_fname = prefix+'tree.pkl'
		self.virus_fname = prefix+'virus.pkl'
		self.frequency_fname = prefix+'frequencies.pkl'
		self.aa_seq_fname = prefix+'aa_seq.pkl'
		self.min_mutation_frequency = min_mutation_frequency
		self.min_genotype_frequency = min_genotype_frequency
		self.time_interval = tuple(time_interval)

		self.auspice_tree_fname = 		'../auspice/data/' + auspice_prefix + 'tree.json'
		self.auspice_sequences_fname = 	'../auspice/data/' + auspice_prefix + 'sequences.json'
		self.auspice_frequency_fname = 	'../auspice/data/' + auspice_prefix + 'frequencies.json'
		self.auspice_meta_fname = 		'../auspice/data/' + auspice_prefix + 'meta.json'
		self.nuc_alphabet = 'ACGT-N'
		self.aa_alphabet = 'ACDEFGHIKLMNPQRSTVWY*X'
		virus_frequencies.__init__(self, **kwargs)

	def dump(self):
		import cPickle
		if hasattr(self, 'tree'):
			with open(self.tree_fname, 'w') as outfile:
				cPickle.dump(self.tree, outfile)
		if hasattr(self, 'viruses'):
			with open(self.virus_fname, 'w') as outfile:
				cPickle.dump(self.viruses, outfile)
		if hasattr(self, 'frequencies'):
			with open(self.frequency_fname, 'w') as outfile:
				cPickle.dump(self.frequencies, outfile)
		if hasattr(self, 'aa_aln'):
			with open(self.aa_seq_fname, 'w') as outfile:
				cPickle.dump(self.aa_aln, outfile)

	def load(self):
		import cPickle
		if os.path.isfile(self.tree_fname):
			with open(self.tree_fname, 'r') as infile:
				self.tree = cPickle.load(infile)
				try:
					self.node_lookup = {l.strain:l for l in self.tree.leaf_iter()}
					self.node_lookup.update({node.strain.lower():node for node in self.tree.leaf_iter()})
				except:
					pass
		if os.path.isfile(self.virus_fname):
			with open(self.virus_fname, 'r') as infile:
				self.viruses = cPickle.load(infile)
				try:
					self.sequence_lookup = {v.strain:v for v in self.viruses}
				except:
					pass
		if os.path.isfile(self.frequency_fname):
			with open(self.frequency_fname, 'r') as infile:
				self.frequencies = cPickle.load(infile)
		if os.path.isfile(self.aa_seq_fname):
			with open(self.aa_seq_fname, 'r') as infile:
				self.aa_aln = cPickle.load(infile)

	def export_to_auspice(self, tree_fields = [], tree_pop_list = [], annotations = []):
		from tree_util import dendropy_to_json, all_descendants
		from io_util import write_json, read_json
		print "--- Streamline at " + time.strftime("%H:%M:%S") + " ---"
		# Move sequence data to separate file
		print "Writing sequences"
		elems = {}
		for node in self.tree:
			if hasattr(node, "clade") and hasattr(node, "aa_seq"):
				elems[node.clade] = node.aa_seq
		write_json(elems, self.auspice_sequences_fname, indent=None)

		print "Writing tree"
		self.tree_json = dendropy_to_json(self.tree.seed_node, tree_fields)
		for node in all_descendants(self.tree_json):
			for attr in tree_pop_list:
				if attr in node:
					node.pop(attr, None)
			if "freq" in node:
				for reg in node["freq"]:
					try:
						node["freq"][reg] = [round(x,3) for x in node["freq"][reg]]
					except:
						node["freq"][reg] = "undefined"				

		if hasattr(self,"clade_designations"):
			# find basal node of clade and assign clade x and y values based on this basal node
			clade_xval = {}
			clade_yval = {}
			for clade, gt in self.clade_designations.iteritems():
				if clade in annotations:
					print "Annotating clade", clade
					base_node = sorted((x for x in self.tree.postorder_node_iter() if all([x.aa_seq[pos-1]==aa for pos, aa in gt])), key=lambda x: x.xvalue)[0]
					clade_xval[clade] = base_node.xvalue
					clade_yval[clade] = base_node.yvalue
			# append clades, coordinates and genotype to meta
			self.tree_json["clade_annotations"] = [(clade, clade_xval[clade],clade_yval[clade], 
								"/".join([str(pos)+aa for pos, aa in gt]))
							for clade, gt in self.clade_designations.iteritems() if clade in annotations
							]
		write_json(self.tree_json, self.auspice_tree_fname, indent=None)
		try:
			read_json(self.auspice_tree_fname)
		except:
			print "Read failed, rewriting with indents"	
			write_json(self.tree_json, self.auspice_tree_fname, indent=1)
			
		# Include genotype frequencies
		if hasattr(self, 'frequencies'):
			write_json(self.frequencies, self.auspice_frequency_fname)

		# Write out metadata
		print "Writing out metadata"		
		meta = {}
		meta["updated"] = time.strftime("X%d %b %Y").replace('X0','X').replace('X','')
		try:
			from pygit2 import Repository, discover_repository
			current_working_directory = os.getcwd()
			repository_path = discover_repository(current_working_directory)
			repo = Repository(repository_path)
			commit_id = repo[repo.head.target].id
			meta["commit"] = str(commit_id)
		except ImportError:
			meta["commit"] = "unknown"
		
		if hasattr(self,"date_region_count"):
			meta["regions"] = self.regions
			meta["virus_stats"] = [ [str(y)+'-'+str(m)] + [self.date_region_count[(y,m)][reg] for reg in self.regions]
									for y,m in sorted(self.date_region_count.keys()) ]
		write_json(meta, self.auspice_meta_fname, indent=0)

	def align(self):
		'''
		aligns viruses using mafft. produces temporary files and deletes those at the end
		after this step, self.viruses is a BioPhython multiple alignment object
		'''
		SeqIO.write([SeqRecord(Seq(v['seq']), id=v['strain']) for v in self.viruses], "temp_in.fasta", "fasta")
		os.system("mafft --nofft temp_in.fasta > temp_out.fasta")
		aln = AlignIO.read('temp_out.fasta', 'fasta')
		for tmp_file in ['temp_in.fasta', 'temp_out.fasta']:
			try:
				os.remove(tmp_file)
			except OSError:
				pass

		self.sequence_lookup = {seq.id:seq for seq in aln}
		# add attributes to alignment
		for v in self.viruses:
			self.sequence_lookup[v['strain']].__dict__.update({k:val for k,val in v.iteritems() if k!='seq'})
		self.viruses = aln

	def infer_tree(self, raxml_time_limit):
		'''
		builds a tree from the alignment using fasttree and RAxML. raxml runs for 
		raxml_time_limit and is terminated thereafter. raxml_time_limit can be 0.
		'''
		def cleanup():
			for file in glob.glob("RAxML_*") + glob.glob("temp*") + ["raxml_tree.newick", "initial_tree.newick"]:
				try:
					os.remove(file)
				except OSError:
					pass

		cleanup()
		AlignIO.write(self.viruses, 'temp.fasta', 'fasta')

		print "Building initial tree with FastTree"
		os.system("fasttree -gtr -nt -gamma -nosupport -mlacc 2 -slownni temp.fasta > initial_tree.newick")
		self.tree = dendropy.Tree.get_from_string(delimit_newick('initial_tree.newick'),'newick', as_rooted=True)
		self.tree.resolve_polytomies()
		self.tree.write_to_path("initial_tree.newick", "newick")

		AlignIO.write(self.viruses,"temp.phyx", "phylip-relaxed")
		if raxml_time_limit>0:
			print "RAxML tree optimization with time limit " + str(raxml_time_limit) + " hours"
			# using exec to be able to kill process
			end_time = time.time() + int(raxml_time_limit*3600)
			process = subprocess.Popen("exec raxml -f d -T 6 -j -s temp.phyx -n topology -c 25 -m GTRCAT -p 344312987 -t initial_tree.newick", shell=True)
			while (time.time() < end_time):
				if os.path.isfile('RAxML_result.topology'):
					break
				time.sleep(10)
			process.terminate()

			checkpoint_files = [file for file in glob.glob("RAxML_checkpoint*")]
			if os.path.isfile('RAxML_result.topology'):
				checkpoint_files.append('RAxML_result.topology')
			if len(checkpoint_files) > 0:
				last_tree_file = checkpoint_files[-1]
				shutil.copy(last_tree_file, 'raxml_tree.newick')
			else:
				shutil.copy("initial_tree.newick", 'raxml_tree.newick')
		else:
			shutil.copy("initial_tree.newick", 'raxml_tree.newick')

		print "RAxML branch length optimization and rooting"
		os.system("raxml -f e -T 6 -s temp.phyx -n branches -c 25 -m GTRGAMMA -p 344312987 -t raxml_tree.newick -o " + self.outgroup['strain'])

		out_fname = "data/tree_infer.newick"
		os.rename('RAxML_result.branches', out_fname)
		Phylo.write(Phylo.read(out_fname, 'newick'),'temp.newick','newick')
		self.tree = self.tree = dendropy.Tree.get_from_string(delimit_newick(out_fname), 'newick', as_rooted=True)
		cleanup()

	def infer_ancestral(self):
		from tree_util import to_Biopython
		from tree_ancestral import ancestral_sequences
		anc_seq = ancestral_sequences(self.tree, self.viruses,seqtype='str')
		anc_seq.calc_ancestral_sequences()

	def temporal_regional_statistics(self):
		'''
		produces a dictionary with (year, month) keys, each entry of which is a
		a dictionary that contains the isolate count in each region observed
		stored as:

		self.date_region_count
		self.regions
		self.region_totals
		'''
		from collections import defaultdict, Counter
		self.date_region_count = defaultdict(lambda:defaultdict(int))
		regions = set()
		# count viruses in every month and every region
		for v in self.viruses:
			if v.strain != self.outgroup['strain']:
				year, month, day = map(int, v.date.split('-'))
				self.date_region_count[(year, month)][v.region]+=1
				regions.add(v.region)
		# add a sorted list of all regions to self and calculate region totals
		self.regions = sorted(regions)
		self.region_totals = {reg:sum(val[reg] for val in self.date_region_count.values()) for reg in self.regions}

	def determine_variable_positions(self):
		'''
		calculates nucleoties_frequencies and aa_frequencies at each position of the alignment
		also computes consensus sequences and position at which the major allele is at less than 1-min_mutation_frequency
		results are stored as
		self.nucleoties_frequencies
		self.aa_frequencies
		self.variable_nucleotides
		self.variable_aa
		'''
		aln_array = np.array(self.viruses)
		self.nucleoties_frequencies = np.zeros((len(self.nuc_alphabet),aln_array.shape[1]))
		for ni,nuc in enumerate(self.nuc_alphabet):
			self.nucleoties_frequencies[ni,:]=(aln_array==nuc).mean(axis=0)

		self.variable_nucleotides = np.where(np.max(self.nucleoties_frequencies,axis=0)<1.0-self.min_mutation_frequency)[0]
		self.consensus_nucleotides = "".join(np.fromstring(self.nuc_alphabet, 'S1')[np.argmax(self.nucleoties_frequencies,axis=0)])

		if hasattr(self, 'aa_aln'):
			aln_array = np.array(self.aa_aln)
			self.aa_frequencies = np.zeros((len(self.aa_alphabet),aln_array.shape[1]))
			for ai,aa in enumerate(self.aa_alphabet):
				self.aa_frequencies[ai,:]=(aln_array==aa).mean(axis=0)

			self.variable_aa = np.where(np.max(self.aa_frequencies,axis=0)<1.0-self.min_mutation_frequency)[0]
			self.consensus_aa = "".join(np.fromstring(self.aa_alphabet, 'S1')[np.argmax(self.aa_frequencies,axis=0)])

	def estimate_frequencies(self, tasks = ['mutations','genotypes', 'clades', 'tree']):
		if 'mutations' in tasks:
			self.all_mutation_frequencies(threshold = self.min_mutation_frequency) 
		if 'genotypes' in tasks:
			self.all_genotypes_frequencies(threshold = self.min_genotype_frequency) 
		if 'clades' in tasks:
			self.all_clade_frequencies() 
		if 'tree' in tasks:
			self.all_tree_frequencies() 
