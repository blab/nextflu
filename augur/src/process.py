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
from itertools import izip

parser = argparse.ArgumentParser(description='Process virus sequences, build tree, and prepare of web visualization')
parser.add_argument('-y', '--years_back', type = float, default=3, help='number of past years to sample sequences from')
parser.add_argument('-v', '--viruses_per_month', type = int, default = 50, help='number of viruses sampled per month')
parser.add_argument('-r', '--raxml_time_limit', type = float, default = 1.0, help='number of hours raxml is run')
parser.add_argument('--interval', nargs = '+', type = float, default = None, help='interval from which to pull sequences')
parser.add_argument('--path', type = str, default = 'data/', help='path of file dumps')
parser.add_argument('--prefix', type = str, default = '', help='prefix of file dumps including auspice')
parser.add_argument('--test', default = False, action="store_true",  help ="don't run the pipeline")
parser.add_argument('--start', default = 'filter', type = str,  help ="start pipeline at specified step")
parser.add_argument('--stop', default = 'export', type=str,  help ="run to end")
parser.add_argument('--skip', nargs='+', type = str,  help ="analysis steps to skip")
parser.add_argument('--ATG', action="store_true", default=False, help ="include full HA sequence starting at ATG")
parser.add_argument('--resolution', type = str,  help ="label for the resolution")


virus_config = {
	'date_format':{'fields':'%Y-%m-%d', 'reg':r'\d\d\d\d-\d\d-\d\d'},
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

def shift_cds(shift, vc, epi_mask, rbs):
	vc['cds'] = (vc['cds'][0]+shift,vc['cds'][1])
	aashift = shift//3
	vc['clade_designations'] = {cl:[(pos-aashift, aa) for pos, aa in gt]
								for cl, gt in vc['clade_designations'].iteritems()}
	return vc, epi_mask[aashift:], [pos-aashift for pos in rbs]

class process(virus_frequencies):
	"""generic template class for processing virus sequences into trees"""
	def __init__(self, path = 'data/', prefix = 'virus', time_interval = (2012.0, 2015.0),
	             run_dir = None, virus = None, resolution = None, date_format={'fields':'%Y-%m-%d', 'reg':r'\d\d\d\d-\d\d-\d\d'},
				 min_mutation_frequency = 0.01, min_genotype_frequency = 0.1, nthreads=1, **kwargs):
		self.path = path
		self.virus_type=virus
		self.resolution = resolution
		self.prefix = prefix
		self.nthreads=nthreads
		if resolution is not None:
			self.resolution_prefix = resolution+'_'
		else:
			self.resolution_prefix = ''
		self.date_format = date_format
		self.min_mutation_frequency = min_mutation_frequency
		self.min_genotype_frequency = min_genotype_frequency
		self.time_interval = tuple(time_interval)
		self.kwargs = kwargs
		self.tree_fname = 		self.path + self.prefix + self.resolution_prefix + 'tree.pkl'
		self.virus_fname = 		self.path + self.prefix + self.resolution_prefix + 'virus.pkl'
		self.frequency_fname = 	self.path + self.prefix + self.resolution_prefix + 'frequencies.pkl'
		self.aa_seq_fname = 	self.path + self.prefix + self.resolution_prefix + 'aa_seq.pkl'
		if run_dir is None:
			import random
			self.run_dir = '_'.join(['temp', time.strftime('%Y%m%d-%H%M%S',time.gmtime()), str(random.randint(0,1000000))])
		else:
			self.run_dir = run_dir
		self.run_dir = self.run_dir.rstrip('/')+'/'
		self.auspice_tree_fname = 		'../auspice/data/' + self.prefix + self.resolution_prefix + 'tree.json'
		self.auspice_sequences_fname = 	'../auspice/data/' + self.prefix + self.resolution_prefix + 'sequences.json'
		self.auspice_frequency_fname = 	'../auspice/data/' + self.prefix + self.resolution_prefix + 'frequencies.json'
		self.auspice_meta_fname = 		'../auspice/data/' + self.prefix + self.resolution_prefix + 'meta.json'
		self.nuc_alphabet = 'ACGT-N'
		self.aa_alphabet = 'ACDEFGHIKLMNPQRSTVWY*X'
		virus_frequencies.__init__(self, **kwargs)

	def make_run_dir(self):
		if not os.path.isdir(self.run_dir):
			try:
				os.makedirs(self.run_dir)
			except OSError as e:
				print "Cannot create run_dir",e

	def remove_run_dir(self):
		if os.path.isdir(self.run_dir):
			import shutil
			shutil.rmtree(self.run_dir)

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

	def export_to_auspice(self, tree_fields = [], tree_pop_list = [], annotations = [], seq='aa'):
		from tree_util import dendropy_to_json, all_descendants
		from io_util import write_json, read_json
		print "--- Streamline at " + time.strftime("%H:%M:%S") + " ---"
		# Move sequence data to separate file
		print "Writing sequences"
		elems = {}
		for node in self.tree:
			if hasattr(node, "clade") and hasattr(node, "seq"):
				if seq == 'nuc':
					elems[node.clade] = {pos:state for pos, (state, ancstate) in
								enumerate(izip(node.seq, self.tree.seed_node.seq)) if state!=ancstate}
				else:
					elems[node.clade] = {pos:state for pos, (state, ancstate) in
								enumerate(izip(node.aa_seq, self.tree.seed_node.aa_seq)) if state!=ancstate}

		if seq == 'nuc':
			elems['root'] = self.tree.seed_node.seq
		else:
			elems['root'] = self.tree.seed_node.aa_seq
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
			clade_present = {}
			clade_xval = {}
			clade_yval = {}
			for clade, gt in self.clade_designations.iteritems():
				if clade in annotations:
					print "Annotating clade", clade
					tmp_nodes = sorted((node for node in self.tree.postorder_node_iter()
						if not node.is_leaf() and all([node.aa_seq[pos-1]==aa for pos, aa in gt])),
						key=lambda node: node.xvalue)
					if len(tmp_nodes):
						clade_present[clade] = True
						base_node = tmp_nodes[0]
						clade_xval[clade] = base_node.xvalue
						clade_yval[clade] = base_node.yvalue
					else:
						clade_present[clade] = False
						print "clade",clade, gt, "not in tree"
			# append clades, coordinates and genotype to meta
			self.tree_json["clade_annotations"] = [(clade, clade_xval[clade],clade_yval[clade],
								"/".join([str(pos)+aa for pos, aa in gt]))
							for clade, gt in self.clade_designations.iteritems() if clade in annotations and clade_present[clade] == True
							]
		write_json(self.tree_json, self.auspice_tree_fname, indent=None)
		try:
			read_json(self.auspice_tree_fname)
		except:
			print "Read failed, rewriting with indents"
			write_json(self.tree_json, self.auspice_tree_fname, indent=1)

		# Include genotype frequencies
		if hasattr(self, 'frequencies'):
			if not hasattr(self, 'aa_entropy') and not hasattr(self, 'nuc_entropy'):
				self.determine_variable_positions()

			if seq=='aa' and hasattr(self, 'aa_entropy'):
				self.frequencies["entropy"] = [ [pos, S, muts] for pos,S,muts in
						izip(xrange(self.aa_entropy.shape[0]), self.aa_entropy,self.variable_aa_identities) ]
			elif seq=='nuc' and hasattr(self, 'nuc_entropy'):
				self.frequencies["entropy"] = [ [pos, S, muts] for pos,S,muts in
						izip(xrange(self.nuc_entropy.shape[0]), self.nuc_entropy,self.variable_nuc_identities) ]

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

	def generate_indexHTML(self):
		htmlpath = '../auspice/'
		if self.virus_type is not None:
			htmlpath+=self.virus_type+'/'
		if self.resolution is not None:
			htmlpath+=self.resolution+'/'

		if not os.path.isdir(htmlpath): os.makedirs(htmlpath)
		if "layout" in self.kwargs:
			tmp_layout=self.kwargs["layout"]
		else:
			tmp_layout="auspice"
		with open(htmlpath+'index.html','w') as out:
			out.write("---\ntitle: nextflu / "+self.virus_type+" / "+self.resolution_prefix.rstrip('_')
					  +"\nlayout: "+tmp_layout
					  +"\nvirus: "+self.virus_type+"\nresolution: "+self.resolution_prefix.rstrip('_')+"\n")
			if "html_vars"  in self.kwargs:
				for vname, val in self.kwargs["html_vars"].iteritems():
					out.write(vname+": "+ val+'\n')
			dt=self.time_interval[1]-self.time_interval[0]
			step = 0.5 if dt<4 else 1 if dt<7 else dt//5
			out.write('---\n\n')
			out.write('<script>\n')
			out.write('var file_prefix = "'+self.prefix+self.resolution_prefix+'";\n')
			out.write('var time_window = '+str(max(1, dt//3))+';\n')
			out.write('var time_ticks=['+', '.join(map(str, np.arange(np.ceil(self.time_interval[0]), np.ceil(self.time_interval[1]), step)))+'];\n')
			if "js_vars" in self.kwargs:
				for vname, val in self.kwargs['js_vars'].iteritems():
					if isinstance(val, basestring):
						out.write('var '+vname+' = "'+val+'";\n')
					else:
						out.write('var '+vname+' = '+str(val)+';\n')
			out.write('{%include '+self.virus_type+'_meta.js %}\n')
			out.write('</script>\n\n')

	def align(self, fast=False):
		'''
		aligns viruses using mafft. produces temporary files and deletes those at the end
		after this step, self.viruses is a BioPhython multiple alignment object
		'''
		self.make_run_dir()
		os.chdir(self.run_dir)
		SeqIO.write([SeqRecord(Seq(v['seq']), id=v['strain']) for v in self.viruses], "temp_in.fasta", "fasta")
		if fast:
			os.system("mafft --anysymbol --thread " + str(self.nthreads) + " temp_in.fasta > temp_out.fasta")
		else:
			os.system("mafft --anysymbol --thread " + str(self.nthreads) + " temp_in.fasta > temp_out.fasta")
		aln = AlignIO.read('temp_out.fasta', 'fasta')
		self.sequence_lookup = {seq.id:seq for seq in aln}
		# add attributes to alignment
		for v in self.viruses:
			self.sequence_lookup[v['strain']].__dict__.update({k:val for k,val in v.iteritems() if k!='seq'})
		self.viruses = aln
		os.chdir('..')
		self.remove_run_dir()

	def infer_tree(self, raxml_time_limit):
		'''
		builds a tree from the alignment using fasttree and RAxML. raxml runs for
		raxml_time_limit and is terminated thereafter. raxml_time_limit can be 0.
		'''
		self.make_run_dir()
		os.chdir(self.run_dir)
		AlignIO.write(self.viruses, 'temp.fasta', 'fasta')

		print "Building initial tree with FastTree"
		os.system("fasttree -gtr -nt -gamma -nosupport temp.fasta > initial_tree.newick")
		self.tree = dendropy.Tree.get_from_string(delimit_newick('initial_tree.newick'),'newick', as_rooted=True)
		self.tree.resolve_polytomies()
		self.tree.write_to_path("initial_tree.newick", "newick")

		AlignIO.write(self.viruses,"temp.phyx", "phylip-relaxed")
		if raxml_time_limit>0:
			print "RAxML tree optimization with time limit " + str(raxml_time_limit) + " hours"
			# using exec to be able to kill process
			end_time = time.time() + int(raxml_time_limit*3600)
			process = subprocess.Popen("exec raxmlHPC -f d -T "+str(self.nthreads) +  " -j -s temp.phyx -n topology -c 25 -m GTRCAT -p 344312987 -t initial_tree.newick", shell=True)
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
		if raxml_time_limit>0:
			os.system("raxmlHPC -f e -T "+str(self.nthreads) +  " -s temp.phyx -n branches -c 25 -m GTRGAMMA -p 344312987 -t raxml_tree.newick -o " + self.outgroup['strain'])
			raxml_rooted=True
		else:
			shutil.copy('raxml_tree.newick', 'RAxML_result.branches')
			raxml_rooted=False

		out_fname = "tree_infer.newick"
		shutil.copy('RAxML_result.branches', out_fname)
		if not raxml_rooted:
			with open(out_fname) as ofile:
				tstr = "".join([x.strip() for x in ofile])
			if tstr.startswith('[&R]'):
				with open(out_fname,'w') as ofile:
					ofile.write(tstr[4:]+'\n')

			T = Phylo.read(out_fname, 'newick')
			try:
				outgroup_clade = [c for x in T.get_terminals() if c.strain == self.outgroup['strain']][0]
			except:
				print("Can't find outgroup in tree -- midpoint_rooting")
				self.midpoint_rooting = True
		else:
			T = Phylo.read(out_fname, 'newick')

		Phylo.write(T,'temp.newick','newick')
		self.tree = dendropy.Tree.get_from_string(delimit_newick(out_fname), 'newick', rooting="force-rooted")
		os.chdir('..')
		self.remove_run_dir()
		if self.midpoint_rooting:
			self.tree.reroot_at_midpoint()


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
		aln_array = np.array(self.nuc_aln)
		self.nuc_frequencies = np.zeros((len(self.nuc_alphabet),aln_array.shape[1]))
		for ni,nuc in enumerate(self.nuc_alphabet):
			self.nuc_frequencies[ni,:]=(aln_array==nuc).mean(axis=0)

		self.variable_nuc = np.where(np.max(self.nuc_frequencies,axis=0)<1.0-self.min_mutation_frequency)[0]
		self.consensus_nuc = "".join(np.fromstring(self.nuc_alphabet, 'S1')[np.argmax(self.nuc_frequencies,axis=0)])
		self.nuc_entropy = -np.sum(self.nuc_frequencies*np.log(np.maximum(1e-10,self.nuc_frequencies)), axis=0)
		self.variable_nuc_identities = [ [self.nuc_alphabet[ii] for ii in np.where(self.nuc_frequencies[:,pos])[0]]
											for pos in xrange(self.nuc_frequencies.shape[1])]

		if hasattr(self, 'aa_aln'):
			aln_array = np.array(self.aa_aln)
			self.aa_frequencies = np.zeros((len(self.aa_alphabet),aln_array.shape[1]))
			for ai,aa in enumerate(self.aa_alphabet):
				self.aa_frequencies[ai,:]=(aln_array==aa).mean(axis=0)

			self.variable_aa = np.where(np.max(self.aa_frequencies,axis=0)<1.0-self.min_mutation_frequency)[0]
			self.consensus_aa = "".join(np.fromstring(self.aa_alphabet, 'S1')[np.argmax(self.aa_frequencies,axis=0)])
			self.aa_entropy = -np.sum(self.aa_frequencies*np.log(np.maximum(1e-10,self.aa_frequencies)), axis=0)
			self.variable_aa_identities = [ [self.aa_alphabet[ii] for ii in np.where(self.aa_frequencies[:,pos])[0]]
											for pos in xrange(self.aa_frequencies.shape[1])]


	def estimate_frequencies(self, tasks = ['mutations','genotypes', 'clades', 'tree']):
		if 'mutations' in tasks:
			self.all_mutation_frequencies(threshold = self.min_mutation_frequency)
		if 'nuc_mutations' in tasks:
			self.all_mutation_frequencies(threshold = self.min_mutation_frequency, nuc=True)
		if 'genotypes' in tasks:
			self.all_genotypes_frequencies(threshold = self.min_genotype_frequency)
		if 'clades' in tasks:
			self.all_clade_frequencies()
		if 'nuc_clades' in tasks:
			self.all_clade_frequencies(nuc=True)
		if 'tree' in tasks:
			self.all_tree_frequencies()
