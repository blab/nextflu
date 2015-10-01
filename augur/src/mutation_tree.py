import time, re, os, argparse,shutil
from tree_refine import tree_refine
from virus_clean import virus_clean
from virus_filter import flu_filter
from date_util import numerical_date
from collections import defaultdict
from process import process, virus_config
from Bio import SeqIO, AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
import numpy as np
from itertools import izip

std_outgroup_file_blast = 'source-data/outgroups.fasta'
std_outgroup_file_nuc = 'source-data/vaccrefMix_HAanno_nuc_Sep2015.fa'
virus_config.update({
	# data source and sequence parsing/cleaning/processing
	'fasta_fields':{0:'strain', 1:'isolate_id', 2:'date',  3:'subtype', 4:'country', 5:'region', 7:'host', 6:'passage'},
	'cds':[0,None], # define the HA start i n 0 numbering
	'verbose':3
	})

def get_date(strain):
	try:
		year = int(strain.split('|')[2][:4])
	except:
		print("cannot parse year of ", strain)
		return 1900

	if year<18:
		year +=2000
	elif year<100:
		year+=1900
	return year

class mutation_tree(process, flu_filter, tree_refine, virus_clean):
	"""docstring for mutation_tree"""
	def __init__(self, aln_fname, outgroup, include_ref_strains = True, outdir = './', formats = ['pdf','png'], verbose = 0, **kwargs):
		process.__init__(self, **kwargs)
		flu_filter.__init__(self, alignment_file = aln_fname, **kwargs)
		tree_refine.__init__(self, **kwargs)
		virus_clean.__init__(self, **kwargs)
		self.midpoint_rooting = False
		self.include_ref_strains = include_ref_strains
		self.verbose = verbose
		self.formats = formats
		self.outdir = outdir.rstrip('/')+'/'
		self.auspice_tree_fname = 		self.outdir + 'tree.json'
		self.auspice_align_fname = 		self.outdir + 'aln.fasta'
		self.auspice_aa_align_fname = 		self.outdir + 'aa_aln.fasta'
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
		elif outgroup=='auto':
			print "automatically determine outgroup"
			self.auto_outgroup_blast()
		elif isinstance(outgroup, basestring):
			seq_names = [x['strain'] for x in self.viruses]
			if outgroup in seq_names:
				self.outgroup = self.viruses.pop(seq_names.index(outgroup))
				if self.verbose:
					print "using outgroup found in alignment", outgroup
			else:
				standard_outgroups = self.load_standard_outgroups()
				if outgroup in standard_outgroups:
					self.outgroup = standard_outgroups[outgroup]
					if self.verbose:
						print "using standard outgroup", outgroup
				else:
					raise ValueError("outgroup %s not found" % outgroup)
					return
		if "anno:" in self.outgroup['desc']:
			anno = [x for x in self.outgroup['desc'].split() if "anno:" in x][0]
			anno = (anno.split(':')[1]).split('_')
			tmp = [(anno[2*i], int(anno[2*i+1])) for i in range(len(anno)/2)]
			self.anno = sorted(tmp, key=lambda x:x[1])
			print("Using annotation",self.anno)
		else:
			self.anno = None
			print("No annotation found")
		#self.anno = sorted((('SP',0), ('HA1',16), ('HA2',329+16)), key=lambda x:x[1])
		self.viruses.append(self.outgroup)
		self.filter_geo(prune=False)
		self.make_strain_names_unique()

	def load_standard_outgroups(self):
		return {'|'.join(seq.description.split()[1].split('|')[:2]).replace(' ',''):{
					'seq':str(seq.seq).upper(), 'strain':seq.description.split()[1].split('|')[1].replace(' ',''), 'desc':seq.description,
						'date':get_date(seq.description)}
				for seq in SeqIO.parse(std_outgroup_file_nuc, 'fasta')}


	def auto_outgroup_blast(self):
		from random import sample
		from Bio.Blast.Applications import NcbiblastxCommandline
		from Bio.Blast import NCBIXML

		self.make_run_dir()
		nvir = 10
		tmp_dates = []
		for v in self.viruses:
			try:
				tmp_dates.append(numerical_date(v["date"]))
			except:
				print("Can't parse date for",v['strain'], v['date'])
		earliest_date = np.min(tmp_dates)
		all_strains = [v["strain"] for v in self.viruses]
		representatives = [SeqRecord(Seq(v['seq']), id=v['strain']) for v in sample(self.viruses, min(nvir, len(self.viruses)))]
		standard_outgroups = self.load_standard_outgroups()
		SeqIO.write(representatives, self.run_dir+'representatives.fasta', 'fasta')
		blast_out = self.run_dir+"outgroup_blast.xml"
		blast_cline = NcbiblastxCommandline(query=self.run_dir+"representatives.fasta", db=std_outgroup_file_blast, evalue=0.01,
		                                     outfmt=5, out=blast_out)
		stdout, stderr = blast_cline()
		with open(blast_out, 'r') as bfile:
			og_blast = NCBIXML.parse(bfile)
			by_og = defaultdict(list)
			for rep in og_blast:
				for hit in rep.alignments:
					for aln in hit.hsps:
						by_og[hit.hit_def].append((rep.query, aln.score, aln.score/aln.align_length, 1.0*aln.identities/aln.align_length))
		by_og = by_og.items()
		print by_og[1]
		# sort by number of hits, then mean score
		by_og.sort(key = lambda x:(len(x[1]), np.mean([y[-2] for y in x[1]])), reverse=True)
		for og, hits in by_og:
			if standard_outgroups[og]['date']<earliest_date-5 or np.mean([y[-1] for y in hits])<0.8:
				break

		if np.mean([y[-1] for y in hits])<0.8:
			self.midpoint_rooting = True

		for ref, hits in by_og:
			if np.max([y[-1] for y in hits])>0.9 and ref!=og:
				self.viruses.append(standard_outgroups[ref])
				print("including reference strain ",ref, [y[-1] for y in hits])
		self.outgroup = standard_outgroups[og]
		self.outgroup['strain']+='OG'
		self.cds = [0,len(self.outgroup['seq'])]
		print("chosen outgroup",self.outgroup['strain'])

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
			if self.anno is not None:
				divides = np.array([x[1] for x in self.anno])
				for node in self.tree.postorder_node_iter():
					node.alt_aa_muts = ""
					tmp = defaultdict(list)
					if len(node.aa_muts):
						for mut in node.aa_muts.split(','):
							anc,pos,der = mut[0], int(mut[1:-1]), mut[-1]
							ii = divides.searchsorted(pos)-1
							if ii>0:
								tmp[ii].append(anc+str(pos-divides[ii])+der)
						for ii, anno in enumerate(self.anno):
							if len(tmp[ii]):
								node.alt_aa_muts+=anno[0]+': '+','.join(tmp[ii])+" "

		self.layout()
		for v in self.viruses:
			if v.strain in self.node_lookup:
				node = self.node_lookup[v.strain]
				for attr in ['strain', 'desc']:
					try:
						node.__setattr__(attr, v.__getattribute__(attr))
					except:
						pass
		# make an amino acid aligment
		from Bio.Align import MultipleSeqAlignment
		from Bio.Seq import Seq
		from Bio.SeqRecord import SeqRecord
		if self.cds is not None:
			tmp_aaseqs = [SeqRecord(Seq(node.aa_seq), id=node.strain, annotations = {'num_date':node.num_date, 'region':node.region}) for node in self.tree.leaf_iter()]
			tmp_aaseqs.sort(key = lambda x:x.annotations['num_date'])
			self.aa_aln = MultipleSeqAlignment(tmp_aaseqs)
		tmp_nucseqs = [SeqRecord(Seq(node.seq), id=node.strain, annotations = {'num_date':node.num_date, 'region':node.region}) for node in self.tree.leaf_iter()]
		tmp_nucseqs.sort(key = lambda x:x.annotations['num_date'])
		self.nuc_aln = MultipleSeqAlignment(tmp_nucseqs)



	def export(self):
		from bio_draw import muttree_draw
		def select_fontsize(n):
			if n<10:
				return 12
			elif n<50:
				return 10
			else:
				return 8

		def branch_label_func(n):
			max_muts = 5
			if hasattr(n,'aa_muts'):
				if alt:
					muts = n.alt_aa_muts
				else:
					muts = n.aa_muts
			else:
				print(n,"has no amino acid mutations")
				try:
					muts = n.nuc_muts
				except:
					print(n,"has no nucleotide mutations")
					muts = []
			tmp = muts.split(',')
			if len(tmp)>max_muts:
				return ', '.join(tmp[:max_muts])+' + '+str(len(tmp)-max_muts)+' others'
			else:
				return ', '.join(tmp)

		from Bio import Phylo
		import matplotlib.pyplot as plt
		plt.rcParams.update({'font.size':select_fontsize(len(self.viruses))})
		plt.ioff()
		from tree_util import to_Biopython
		tmp_tree = to_Biopython(self.tree)
		tmp_tree.ladderize()

		for alt in [False] if self.anno is None else [False, True]:
			fig = plt.figure('Tree', figsize = (15,2+len(self.viruses)/5))
			ax = plt.subplot('111')
			muttree_draw(tmp_tree, axes=ax, show_confidence=False, do_show=False,
				label_func = lambda x: x.name,
				branch_labels = branch_label_func
				)
			ax.invert_yaxis()
			tl = np.diff(ax.get_xticks())[0]
			lengthbar = tl/2
			plt.plot( [0,lengthbar],[len(self.viruses),len(self.viruses)], lw=10, c='k')
			plt.text(lengthbar/2, len(self.viruses)+0.1, str(lengthbar),horizontalalignment='center',fontsize=16)
			ax.set_axis_off()
			for fmt in self.formats:
				if alt:
					plt.savefig(self.outdir+'tree_alt.'+fmt)
				else:
					plt.savefig(self.outdir+'tree.'+fmt)

			for t in tmp_tree.find_clades():
				t.label = t.name # save original name
				if hasattr(t,"strain"):
					t.name = t.strain
				else:
					t.name = ""
				if alt:
					muts = t.alt_aa_muts if hasattr(t,'alt_aa_muts') else t.nuc_muts
				else:
					muts = t.aa_muts if hasattr(t,'aa_muts') else t.nuc_muts
				if len(t.name) and len(muts): t.name+='-'
				t.name+='_'.join(muts.split(',')).replace(' ','')

			if alt:
				Phylo.write(tmp_tree, self.outdir+'tree_alt.nwk', 'newick')
			else:
				Phylo.write(tmp_tree, self.outdir+'tree.nwk', 'newick')
			for t in tmp_tree.find_clades(): # revert to original name
				t.name = t.label
			plt.close('Tree')
		for n in self.tree.leaf_iter():
			for field in self.fasta_fields.values():
				if (not hasattr(n, field)) or n.__dict__[field]=="":
					n.__dict__[field]="Unknown"
		for n in self.tree.postorder_internal_node_iter():
			for field in self.fasta_fields.values():
				n.__dict__[field]="Unknown"


		if self.cds is None:
			self.export_to_auspice(tree_fields = ['nuc_muts','num_date']+self.fasta_fields.values(), seq='nuc')
		else:
			self.export_to_auspice(tree_fields = ['aa_muts','alt_aa_muts','num_date']+self.fasta_fields.values())

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
		AlignIO.write(self.viruses, self.auspice_align_fname, 'fasta')
		self.remove_insertions()
		print "--- Tree	 infer at " + time.strftime("%H:%M:%S") + " ---"
		self.infer_tree(raxml_time_limit)
		print "--- Infer ancestral sequences " + time.strftime("%H:%M:%S") + " ---"
		self.infer_ancestral()  # -> every node has a sequence
		print "--- Tree refine at " + time.strftime("%H:%M:%S") + " ---"
		self.refine()
		if self.cds:
			aa_aln = MultipleSeqAlignment([])
			for node in self.tree.leaf_iter():
				aa_aln.append(SeqRecord(id=node.strain, seq=Seq(node.aa_seq)))
			AlignIO.write(aa_aln, self.auspice_aa_align_fname, 'fasta')


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
			os.makedirs(params.out+'/css')
		except OSError as e:
			print "Cannot create output directory",e
	virus_config["outdir"]=params.out

	muttree = mutation_tree(params.aln, params.outgroup, **virus_config)
	muttree.run(raxml_time_limit=0.1)
	muttree.export()


	shutil.copy2('../auspice/_site/js/muttree.js', muttree.outdir+'js/muttree.js')
	shutil.copy2('../auspice/_site/js/msa.min.js', muttree.outdir+'js/msa.min.js')
	shutil.copy2('../auspice/_site/muttree/index.html', muttree.outdir+'index.html')
	shutil.copy2('../auspice/_site/css/style.css', muttree.outdir+'css/style.css')

	with open(muttree.outdir+'/js/fields.js', 'w') as ofile:
		for field in ['passage', 'host', 'subtype','region']:
			try:
				tmp = sorted(set([x.__dict__[field] for x in muttree.tree.leaf_iter()]))
			except:
				tmp = ["Unknown"]
			if "Unknown" not in tmp: tmp.append("Unknown")
			ofile.write(field + 's = [' + ', '.join(map(lambda x:'"'+str(x)+'"',tmp))+']\n')

#	os.system('firefox '+muttree.outdir+'index.html &')
