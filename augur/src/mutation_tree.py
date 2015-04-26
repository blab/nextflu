import time, re, os, argparse,shutil
from tree_refine import tree_refine
from virus_clean import virus_clean
from virus_filter import flu_filter
from collections import defaultdict
from process import process, virus_config
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Align import MultipleSeqAlignment
import numpy as np
from itertools import izip

std_outgroup_file = 'source-data/outgroups.fasta'
virus_config.update({
	# data source and sequence parsing/cleaning/processing
	'fasta_fields':{0:'strain', 1:'date', 2:'isolate_id', 3:'passage', 4:'subtype', 5:'ori_lab', 6:'sub_lab', 7:'submitter'},
	'cds':[0,None], # define the HA start i n 0 numbering
	'auspice_prefix':'H1N1pdm_',
	'verbose':3
	})


class mutation_tree(process, flu_filter, tree_refine, virus_clean):
	"""docstring for mutation_tree"""
	def __init__(self, aln_fname, outgroup, outdir = './', formats = ['pdf','svg','png'], verbose = 0, **kwargs):
		process.__init__(self, **kwargs)
		flu_filter.__init__(self, alignment_file = aln_fname, **kwargs)
		tree_refine.__init__(self, **kwargs)
		virus_clean.__init__(self, **kwargs)
		self.verbose = verbose
		self.formats = formats
		self.outdir = outdir.rstrip('/')+'/'
		self.auspice_tree_fname = 		self.outdir + 'tree.json'
		self.auspice_sequences_fname = 	self.outdir + 'sequences.json'
		self.auspice_frequencies_fname = None
		self.auspice_meta_fname = 		self.outdir + 'meta.json'

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
				standard_outgroups = [{'strain':seq.name, 'seq':str(seq.seq).upper(), 'desc':seq.description}
										for seq in SeqIO.parse(std_outgroup_file, 'fasta') ]
				outgroup_names = [x['strain'] for x in standard_outgroups]
				if outgroup in outgroup_names:
					self.outgroup = standard_outgroups[outgroup_names.index(outgroup)]
					if self.verbose:
						print "using standard outgroup", outgroup
				else:
					raise ValueError("outgroup %s not found" % outgroup)
					return
		self.viruses.append(self.outgroup)
		self.filter_geo(prune=False)
		self.make_strain_names_unique()

	def refine(self):
		self.node_lookup = {node.taxon.label:node for node in self.tree.leaf_iter()}
		self.unique_date()
		self.remove_outgroup()
		self.ladderize()
		self.collapse()
		self.add_nuc_mutations()
		self.add_node_attributes()
		if self.cds is not None:
			self.translate_all()
			self.add_aa_mutations()
		self.layout()
		for v in self.viruses:
			if v.strain in self.node_lookup:
				node = self.node_lookup[v.strain]
				for attr in ['strain', 'desc']:
					try:
						node.__setattr__(attr, v.__getattribute__(attr))
					except:
						pass

	def export(self):
		from bio_draw import muttree_draw
		def select_fontsize(n):
			if n<10:
				return 12
			elif n<50:
				return 10
			else:
				return 8


		from Bio import Phylo
		import matplotlib.pyplot as plt
		plt.rcParams.update({'font.size':select_fontsize(len(self.viruses))})
		plt.ioff()
		from tree_util import to_Biopython
		tmp_tree = to_Biopython(self.tree)
		tmp_tree.ladderize()
		fig = plt.figure('Tree', figsize = (15,2+len(self.viruses)/5))
		ax = plt.subplot('111')

		muttree_draw(tmp_tree, axes=ax, show_confidence=False, do_show=False,
			label_func = lambda x: x.name,
			branch_labels = lambda x: x.aa_muts if hasattr(x,'aa_muts') else x.nuc_muts
			)
		ax.invert_yaxis()
		tl = np.diff(ax.get_xticks())[0]
		lengthbar = tl/2
		plt.plot( [0,lengthbar],[len(self.viruses),len(self.viruses)], lw=10, c='k')
		plt.text(lengthbar/2, len(self.viruses)+0.1, str(lengthbar),horizontalalignment='center',fontsize=16)
		ax.set_axis_off()
		for fmt in self.formats:
			plt.savefig(self.outdir+'tree.'+fmt)

		for t in tmp_tree.find_clades():
			if t.name is None:
				t.name=''
			muts = t.aa_muts if hasattr(t,'aa_muts') else t.nuc_muts
			if len(t.name) and len(muts): t.name+='-'
			t.name+='_'.join(muts.split(','))

		Phylo.write(tmp_tree, self.outdir+'tree.nwk', 'newick')

		self.export_to_auspice(tree_fields = ['aa_muts','num_date']+self.fasta_fields.values())

	def make_strain_names_unique(self):
		strain_to_seq = defaultdict(list)
		for v in self.viruses:
			strain_to_seq[v['strain']].append(v)
		for strain, strain_list in strain_to_seq.iteritems():
			if len(strain_list)>1:
				for ii, virus in enumerate(strain_list):
					virus['strain']+='-'+str(ii+1)



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
	parser.add_argument('--out', type = str, default = 'output/', help='output directory')
	params = parser.parse_args()

	# check and parse cds
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

	# check and create output directory
	if not os.path.isdir(params.out):
		try:
			os.makedirs(params.out)
			os.makedirs(params.out+'/js')
		except OSError as e:
			print "Cannot create output directory",e
	virus_config["outdir"]=params.out

	muttree = mutation_tree(params.aln, params.outgroup, **virus_config)
	muttree.run(raxml_time_limit=0.1)
	muttree.export()

	shutil.copy2('../auspice/_site/js/muttree.js', muttree.outdir+'js/muttree.js')
	shutil.copy2('../auspice/_site/muttree/index.html', muttree.outdir+'muttree.html')
	shutil.copytree('../auspice/_site/css', muttree.outdir+'css')


#	os.system('firefox '+muttree.outdir+'muttree.html &')
