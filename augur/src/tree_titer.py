# clean, reroot, ladderize newick tree
# output to tree.json
import numpy as np
from io_util import *
from tree_util import *
import time

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
		return ref_strains, np.log2(np.array(ref_matrix)), test_strains, np.log2(np.matrix(test_matrix))


def normalize_HI_matrices(ref_matrix, test_matrix):
	norm_test_matrix = test_matrix - np.diag(ref_matrix)
	norm_ref_matrix = ref_matrix - np.diag(ref_matrix)
	return -norm_ref_matrix, -norm_test_matrix

def load_HI_distances(HI_fname):
	pass

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

def get_path(tree, v1, v2):
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

def read_tables():
	import glob
	from itertools import product
	flist = glob.glob('../../../flu_HI_data/tables/Flu*csv')
	all_names = set()
	all_measurements = set()
	HI_matrices = []
	for fname in flist:
		ref_names, ref_matrix, test_names, test_matrix = parse_HI_matrix(fname)
		HI_matrices.append( (ref_names, ref_matrix, test_names, test_matrix ))
		all_names.update(ref_names)	
		all_names.update(test_names)
		all_measurements.update([ (test, ref) for ref, test in product(test_names, ref_names)])

	return all_names, all_measurements, HI_matrices

def prune_tree(tree, taxa_with_HI):
	seq_strains = set([leaf.label for leaf in tree.leaf_iter()])
	print seq_strains.intersection(taxa_with_HI)
	tree.retain_taxa_with_labels(taxa_with_HI)


def main(tree_fname = 'data/tree_ancestral.json', HI_fname='data/HI_titers.txt'):

	print "--- Fitting HI titers at " + time.strftime("%H:%M:%S") + " ---"

	tree =  json_to_dendropy(read_json(tree_fname))
	HI_distances = load_HI_distances(HI_fname)



if __name__ == "__main__":
	main()
