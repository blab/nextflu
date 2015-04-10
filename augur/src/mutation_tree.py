import time, re, os, argparse
from tree_refine import tree_refine
from virus_clean import virus_clean
from process import process, virus_config
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Align import MultipleSeqAlignment
import numpy as np
from itertools import izip

std_outgroup_file = 'source-data/outgroups.fasta'
virus_config.update({
	# data source and sequence parsing/cleaning/processing
	'cds':[0,None], # define the HA start i n 0 numbering
	'auspice_prefix':'H1N1pdm_',
	'verbose':3
	})


class mutation_tree(process, tree_refine, virus_clean):
	"""docstring for mutation_tree"""
	def __init__(self, aln_fname, outgroup, formats = ['pdf','svg','png'], verbose = 0, **kwargs):
		self.verbose = verbose
		self.formats = formats
		process.__init__(self, **kwargs)
		tree_refine.__init__(self, **kwargs)
		virus_clean.__init__(self, **kwargs)
		if os.path.isfile(aln_fname):
			self.aln_fname = aln_fname
			try:
				self.viruses = [{'strain':seq.name, 'seq':str(seq.seq).upper(), 'desc':seq.description}
								for seq in SeqIO.parse(self.aln_fname, 'fasta') ]
			except:
				raise ValueError("Parsing of fasta file %s failed!" % self.aln_fname)
		else:
			raise ValueError("file %s not found" % aln_fname)		
		if os.path.isfile(outgroup):
			tmp = [{'strain':seq.name, 'seq':str(record.seq).upper(), 'desc':seq.description}
								for seq in SeqIO.parse(outgroup, 'fasta') ]			
			if len(tmp):
				self.outgroup = tmp[0]
				if len(tmp)>1:
					print "More than one sequence in ", outgroup, "taking first"
				if self.verbose:
					print "using outgroup found in file ", outgroup
		elif isinstance(outgroup, basestring):
			seq_names = [x['strain'] for x in self.viruses]
			if outgroup in seq_names:
				self.outgroup = self.viruses.pop(seq_names.index(outgroup))
				if self.verbose:
					print "using outgroup found in alignment", outgroup
			else:
				standard_outgroups = [{'strain':seq.name, 'seq':str(record.seq).upper(), 'desc':seq.description}
										for seq in SeqIO.parse(std_outgroup_file, 'fasta') ]
				outgroup_names = [x['strain'] for x in standard_outgroups]
				if outgroup in outgroup_names:
					self.outgroup = standard_outgroups.index(outgroup)
					if self.verbose:
						print "using standard outgroup", outgroup
				else:
					raise ValueError("outgroup %s not found" % outgroup)
					return
		self.viruses.append(self.outgroup)

	def refine(self):
		self.node_lookup = {node.taxon.label:node for node in self.tree.leaf_iter()}
		self.remove_outgroup()
		self.ladderize()
		self.collapse()
		self.add_nuc_mutations()
		if self.cds is not None:
			self.translate_all()
			self.add_aa_mutations()
		self.layout()

	def export(self):
		from Bio import Phylo
		import matplotlib.pyplot as plt
		plt.rcParams.update({
          	'text.fontsize': 16,
          	'font.size':10,
          	'line.width':5,
			'font.sans-serif': 'Helvetica',}
		)
		
		from tree_util import to_Biopython
		tmp_tree = to_Biopython(self.tree)
		tmp_tree.ladderize()
		fig = plt.figure('Tree', figsize = (15,len(self.viruses)/3))
		ax = plt.subplot('111')

		Phylo.draw(tmp_tree, axes=ax, show_confidence=False, 
			label_func = lambda x: x.name,
			branch_labels = lambda x: x.aa_muts if hasattr(x,'aa_muts') else x.nuc_muts
			)
		tl = np.diff(ax.get_xticks())[0]
		lengthbar = tl/2
		plt.plot( [0,lengthbar],[len(self.viruses),len(self.viruses)], lw=10, c='k')
		plt.text(lengthbar/2, len(self.viruses)-1, str(lengthbar),horizontalalignment='center',fontsize=16)
		ax.set_axis_off()

		for t in tmp_tree.find_clades():
			if t.name is None:
				t.name=''
			muts = t.aa_muts if hasattr(t,'aa_muts') else t.nuc_muts
			if len(t.name) and len(muts): t.name+='-'
			t.name+='_'.join(muts.split(','))
		for fmt in self.formats:
			plt.savefig(self.prefix+'tree.'+fmt)

		Phylo.write(tmp_tree, self.prefix+'tree.nwk', 'newick')


	def run(self, raxml_time_limit):
		self.align()
		self.remove_insertions()
		print "--- Tree	 infer at " + time.strftime("%H:%M:%S") + " ---"
		self.infer_tree(raxml_time_limit)
		print "--- Infer ancestral sequences " + time.strftime("%H:%M:%S") + " ---"
		self.infer_ancestral()  # -> every node has a sequence
		print "--- Tree refine at " + time.strftime("%H:%M:%S") + " ---"
		self.refine()



if __name__=="__main__":
	parser = argparse.ArgumentParser(description='Build a tree given a fasta file and annotate braches with mutations')
	parser.add_argument('--aln', required = True, type = str,  help ="fasta file with input sequences")
	parser.add_argument('--outgroup', required = True, type = str,  help ="outgroup to root the tree, strain label or fasta file")
	parser.add_argument('--cds', nargs = '+', type = int, default = None, help='part of the outgroup sequence that is to be translated')
	parser.add_argument('--prefix', type = str, default = 'data/', help='path+prefix of file dumps')
	params = parser.parse_args()

	if params.cds is None:
		virus_config['cds']=None
	else:
		if len(params.cds)==2:
			virus_config['cds']=params.cds
		elif len(params.cds)==1:			
			virus_config['cds']=(params.cds[0], None)
		else:
			raise ValueError("Expecting a cds of length 1 (start only) or 2, got "+str(params.cds))
			exit()

	virus_config["prefix"]=params.prefix

	muttree = mutation_tree(params.aln, params.outgroup, **virus_config)
	muttree.run(raxml_time_limit=0.1)
	muttree.export()