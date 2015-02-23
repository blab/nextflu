# clean, reroot, ladderize newick tree
# output to tree.json
import numpy as np
from io_util import *
from tree_util import *
import time
from collections import defaultdict

min_titer = 10.0

def strain_name_fixing(name):
	return name.replace(' ','').lower()

def parse_HI_matrix(fname):
	import csv
	def titer_to_number(val):
		try:
			if '<' in val:
				return min_titer
			else:
				return float(val)
		except:
			print "Bad HI measurement:", val
			return np.nan

	with open(fname) as infile:
		csv_reader = csv.reader(infile)

		# parse sera
		row1 = csv_reader.next()
		row2 = csv_reader.next()
		print row1
		print row2
		ref_sera = [(e1+'/'+e2).replace(' ','') for e1,e2 in zip(row1, row2)[4:]]
		print ref_sera
		for row in csv_reader: # advance until the reference virus
			if row[0].startswith('REFERENCE'):
				break

		ref_strains = []
		ref_matrix = []
		for row in csv_reader: 
			if row[0].startswith('TEST'):
				break
			else: # load matrices until the test virus section starts
				ref_strains.append(strain_name_fixing(row[0]))
				ref_matrix.append(map(titer_to_number, row[4:]))

		test_strains = []
		test_matrix = []
		for row in csv_reader: # load test viruses until it is no longer an A/ flu  name
			if not row[0].startswith('A/'):
				break
			else:
				test_strains.append(strain_name_fixing(row[0]))
				test_matrix.append(map(titer_to_number, row[4:]))
		return ref_strains, np.array(ref_matrix), test_strains, np.matrix(test_matrix)


def normalize_HI_matrices(ref_matrix, test_matrix):
	self_avidity = np.diag(ref_matrix)
	norm_test_matrix = np.array(test_matrix[:,:len(self_avidity)])
	norm_ref_matrix = np.array(ref_matrix[:,:len(self_avidity)])
	for ai, av in enumerate(self_avidity):
		print av
		norm_test_matrix[:,ai] -= av
		norm_ref_matrix[:,ai]  -= av

	return -norm_ref_matrix, -norm_test_matrix

def read_tables():
	import glob
	from itertools import product
	flist = glob.glob('../../../flu_HI_data/tables/Flu*csv')
	all_names = set()
	all_measurements = defaultdict(list)
	HI_matrices = []
	for fname in flist:
		ref_names, ref_matrix, test_names, test_matrix = parse_HI_matrix(fname)
#		try:
#			ref_matrix, test_matrix = normalize_HI_matrices(ref_matrix, test_matrix)
#		except:
#			continue
		HI_matrices.append( (ref_names, ref_matrix[:,:len(ref_names)], test_names, test_matrix[:,:len(ref_names)] ))
		all_names.update(ref_names)	
		all_names.update(test_names)

		for test, ref in product(test_names, ref_names):
			all_measurements[(test, ref)].append(test_matrix[test_names.index(test), ref_names.index(ref)])
		for ref_test, ref in product(ref_names, ref_names):
			all_measurements[(ref_test, ref)].append(ref_matrix[ref_names.index(ref_test), ref_names.index(ref)])

		print fname, len(all_measurements), "measurements"

	print "total from tables:", len(all_measurements), "measurements"
	trevor_table = '../../../flu_HI_data/tables/trevor_elife_H3N2_HI_data.tsv'
	import csv
	with open(trevor_table) as infile:
		table_reader = csv.reader(infile, delimiter="\t")
		header = table_reader.next()
		for row in table_reader:
			try:
				all_names.add(strain_name_fixing(row[1]))
				all_measurements[(strain_name_fixing(row[1]), strain_name_fixing(row[4]))].append(float(row[6]))
			except:
				pass #print row

	print "grand total:", len(all_measurements), "measurements"
	return all_names, all_measurements, HI_matrices


def assign_tree_graph(tree, HI_distances):
	node_count = 0
	for node in tree.postorder_node_iter():
		node.node_index = node_count
		node_count+=1

	tree_graph = np.zeros((len(HI_distances), node_count))
	for di, ((v1, v2), dist) in enumerate(HI_distances):
		tmp_path = get_path(tree, v1,v2)
		for n in tmp_path:
			node[di, n.node_index]=1
	return tree_graph

def get_path_dendropy(tree, v1, v2):
	p1 = [v1]
	p2 = [v2]
	while p1[-1].parent != tree.seed_node:
		p1.append(p1[-1].parent_node)
	p1.append(tree.seed_node)
	p1.reverse()

	while p2[-1].parent != tree.seed_node:
		p2.append(p2[-1].parent_node)
	p2.append(tree.seed_node)
	p2.reverse()

	for pi, tmp_v1, tmp_v2 in enumerate(izip(p1,p2)):
		if tmp_v1!=tmp_v2:
			break
	path = p1[pi-1:] + p2[pi:].reverse()

	return path

def get_path_biopython(tree, v1, v2):
	#print "getting path for", v1.name, v2.name
	p1 = tree.get_path(v1)
	p2 = tree.get_path(v2)
	path = None
	for ii, (c1,c2) in enumerate(zip(p1, p2)):
		if c1!=c2:
			#print "found match", ii
			path = p1[ii:] + p2[ii:]
			break

	return path

def calc_tree_HI_distance(v1,v2, tree, lca):
	dist = 0
	tmp_lca = lca[(v1,v2)]
	for tmp_v in [v1, v2]:
		while tmp_v.parent_node != tmp_lca:
			dist+=tmp_v.dHI
			tmp_v1 = tmp_v.parent_node

	return dist

def fit_func(dHI, tree_graph, distances):
	return np.sum( (distances - np.dot(tree_graph, dHI))**2 )


def prune_tree(tree, taxa_with_HI):
	seq_strains = set([leaf.label for leaf in tree.leaf_iter()])
	print seq_strains.intersection(taxa_with_HI)
	tree.retain_taxa_with_labels(taxa_with_HI)


def main(tree_fname = 'data/tree_ancestral.json', HI_fname='data/HI_titers.txt'):

	print "--- Fitting HI titers at " + time.strftime("%H:%M:%S") + " ---"

	tree =  json_to_dendropy(read_json(tree_fname))
	HI_distances = load_HI_distances(HI_fname)

def get_strains_with_HI_and_sequence():
	names, measurements, HI_matrices = read_tables()
	from Bio import SeqIO
	good_strains = set()
	with open("data/strains_with_HI.fasta", 'w') as outfile, \
		 open("data/20150222_all_H3N2_HA1.fasta", 'r') as infile:
		for seq_rec in SeqIO.parse(infile, 'fasta'):
			reduced_name = strain_name_fixing(seq_rec.name)
			if reduced_name in names and (reduced_name not in good_strains):
				SeqIO.write(seq_rec, outfile,'fasta')
				good_strains.add(reduced_name)
				print seq_rec.name

def fix_HI_tree():
	from Bio import Phylo
	tree = Phylo.read('data/HI_tree.newick', 'newick')
	early_leafs = [leaf for leaf in tree.get_terminals() if '1991' in leaf.name]
	tree.root_with_outgroup(early_leafs[0])
	tree.ladderize()
	Phylo.write(tree, 'data/HI_tree.newick',  'newick')
	return tree

if __name__ == "__main__":
	from Bio import Phylo
#	get_strains_with_HI_and_sequence()
#	tree = fix_HI_tree()
	tree = Phylo.read('data/HI_tree.newick', 'newick')
	tree.root.branch_length = 0.001
	for node in tree.get_nonterminals():
		if node.branch_length<0.00015:
			node.branch_length=0
			tree.collapse(node)
	tree.ladderize()

	names, measurements, tables = read_tables()
	normalized_measurements = {}
	ref_names = set()
	for (test, ref), val in measurements.iteritems():
		if (ref,ref) in measurements:
			ref_names.add(ref)
			normalized_val = np.log2(measurements[(ref, ref)]).mean() - np.log2(val).mean()
			if normalized_val<=10:
				normalized_measurements[(test, ref)] = max(0,normalized_val)

	names_to_clades = {leaf.name.lower(): leaf for leaf in tree.get_terminals()}
	# assign indices to branches
	branch_count = 0
	for node in tree.get_nonterminals():
		for c in node.clades:
			if True: #not c.is_terminal():
				c.branch_index = branch_count
				branch_count+=1

	tree_graph = []
	HI_diff = []

	for mi, (pair, val) in enumerate(normalized_measurements.iteritems()):
		if not np.isnan(val):
			try:
				if pair[0] != pair[1] and pair[0] in names_to_clades and pair[1] in names_to_clades:
					path = get_path_biopython(tree, names_to_clades[pair[0]], names_to_clades[pair[1]])
					tmp = np.zeros(branch_count)
					tmp[[c.branch_index for c in path]] = 1
					tree_graph.append(tmp)
					HI_diff.append(val)
				else:
					pair, "not found"
			except:
				print pair, "ERROR"

	HI_diff =  np.array(HI_diff)
	tree_graph= np.array(tree_graph)

	from scipy.optimize import nnls
	w = nnls(tree_graph, HI_diff)
	tree.root.cHI = 0
	for node in tree.get_nonterminals():
		node.ref=False
		for c in node.clades:
			if True: #not c.is_terminal():
				c.constraints = tree_graph[:,c.branch_index].sum()
				if c.constraints>0:
					c.dHI = w[0][c.branch_index]
				else:
					c.dHI = 0
			else:
				c.dHI=0
			c.cHI = node.cHI + c.dHI
			c.ref=False
			c.branch_length = c.dHI

	for ref in ref_names:
		if ref in names_to_clades:
			names_to_clades[ref].ref = True

	color_BioTree_by_attribute(tree, "cHI", transform = lambda x:x)
	Phylo.draw(tree, label_func = lambda  x: x.name if x.ref else '', show_confidence= False)
