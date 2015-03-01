import time, os, argparse,shutil,subprocess, glob
from Bio import SeqIO, AlignIO,Phylo
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import dendropy
from tree_util import delimit_newick
from StringIO import StringIO
from itertools import izip

class process(object):
	"""generic template class for processing virus sequences into trees"""
	def __init__(self, tree_fname = 'data/tree.pkl', virus_fname = 'data/virus.pkl', 
				frequency_fname = 'data/frequency.pkl',**kwargs):
		self.tree_fname = tree_fname
		self.virus_fname = virus_fname
		self.frequency_fname = frequency_fname

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

	def load(self):
		import cPickle
		if os.path.isfile(self.tree_fname):
			with open(self.tree_fname, 'r') as infile:
				self.tree = cPickle.load(infile)
		if os.path.isfile(self.virus_fname):
			with open(self.virus_fname, 'r') as infile:
				self.viruses = cPickle.load(infile)
		if os.path.isfile(self.frequency_fname):
			with open(self.frequency_fname, 'r') as infile:
				self.frequencies = cPickle.load(infile)

	def align(self):
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
		anc_seq = ancestral_sequences(to_Biopython(self.tree), self.viruses,seqtype='str')
		anc_seq.calc_ancestral_sequences()
		# copy the inferred sequences into the  biopython tree
		for node, anc_node in izip(self.tree.postorder_internal_node_iter(), anc_seq.T.get_nonterminals(order='postorder')):
			node.seq = anc_node.seq
		for node, anc_node in izip(self.tree.leaf_iter(), anc_seq.T.get_terminals()):
			node.seq = anc_node.seq


