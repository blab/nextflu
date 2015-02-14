# clean, reroot, ladderize newick tree
# output to tree.json

import os, re, time
import dendropy
from io_util import *
from seq_util import *
from date_util import *
from tree_util import *

OUTGROUP = 'A/Beijing/32/1992'

def delimit_newick(infile_name):
	with open(infile_name, 'r') as file:
		newick = file.read().replace('\n', '')
		newick = re.sub(r'(A/[^\:^,^)]+)', r"'\1'", newick)
	return newick

def crossref_import(branches_tree_file, states_tree_file, states_file):
	"""RAxML won't return a single NEWICK tree with both ancestral states and branch lengths"""
	"""This opens the necessary RAxL output files and outputs a single Dendropy tree"""
	label_to_seq = {}
	with open(states_file) as file:
		for line in file:
			(label, seq) = line.split()
			label_to_seq[label] = seq
	branches_tree = dendropy.Tree.get_from_string(delimit_newick(branches_tree_file), "newick", as_rooted=True)
	states_tree = dendropy.Tree.get_from_string(delimit_newick(states_tree_file), "newick", as_rooted=True)
	for (bn, sn) in zip(branches_tree.postorder_node_iter(), states_tree.postorder_node_iter()):
		if sn.label:
			bn.seq = label_to_seq[sn.label]
	return branches_tree

def get_yvalue(node):
	"""Return y location based on recursive mean of daughter locations"""
	if hasattr(node, 'yvalue'):
		return node.yvalue
	if node.child_nodes():
		mean = 0
		for ch in node.child_nodes():
			mean += get_yvalue(ch)
		return mean / float(len(node.child_nodes()))

def get_xvalue(node):
	"""Return x location based on total distance from root"""
	root = node.get_tree_root()
	return node.get_distance(root)

def remove_outgroup(tree):
	"""Reroot tree to outgroup"""
	outgroup_node = None
	for node in tree.postorder_node_iter():
		if (str(node.taxon) == OUTGROUP):
			outgroup_node = node
	if outgroup_node:
		tree.prune_subtree(outgroup_node)

def collapse(tree):
	"""Collapse short edges to polytomies"""
	for edge in tree.postorder_edge_iter():
		if edge.length < 0.00001 and edge.is_internal():
			edge.collapse()

def reduce(tree):
	"""Remove outlier tips"""
	for node in tree.postorder_node_iter():
		if node.edge_length > 0.01 and node.is_leaf():
			parent = node.parent_node
			parent.remove_child(node)

def ladderize(tree):
	"""Sorts child nodes in terms of the length of subtending branches each child node has"""
	node_desc_counts = {}
	for node in tree.postorder_node_iter():
		if len(node._child_nodes) == 0:
			node_desc_counts[node] = node.edge_length
		else:
			total = 0
			if node.edge_length > 0:
				total += node.edge_length
			for child in node._child_nodes:
				total += node_desc_counts[child]
			node_desc_counts[node] = total
			node._child_nodes.sort(key=lambda n: node_desc_counts[n], reverse=True)

def layout(tree):
	"""Set yvalue of tips by post-order traversal"""
	yvalue = 0
	distance_matrix = dendropy.treecalc.PatristicDistanceMatrix(tree)
	tips = [node for node in tree.leaf_iter()]
	tips[0].yvalue = yvalue
	for (a,b) in zip(tips[:-1], tips[1:]):
		d = distance_matrix(a.taxon, b.taxon)
	#	print str(a.taxon) + " to " + str(b.taxon) + ": " + str(d)
		if b.is_leaf():
			yvalue += d
			b.yvalue = yvalue

	for node in tree.postorder_node_iter():
		node.yvalue = get_yvalue(node)

def add_virus_attributes(viruses, tree):
	"""Add date and loc attributes to all tips in tree"""
	strain_to_date = {}
	strain_to_country = {}
	strain_to_region = {}	
	for v in viruses:
		strain_to_date[v['strain']] = v['date']
		strain_to_country[v['strain']] = v['country']
		strain_to_region[v['strain']] = v['region']				
	for node in tree.postorder_node_iter():
		strain = str(node.taxon).replace("'", '')
		if strain_to_date.has_key(strain):
			node.date = strain_to_date[strain]
		if strain_to_country.has_key(strain):
			node.country = strain_to_country[strain]
		if strain_to_region.has_key(strain):
			node.region = strain_to_region[strain]						

def add_node_attributes(tree):
	"""Add clade, xvalue, yvalue, mutation and trunk attributes to all nodes in tree"""
	clade = 0
	yvalue = 0
	for node in tree.postorder_node_iter():
		node.clade = clade
		clade += 1
		if node.is_leaf():
			node.yvalue = yvalue
			yvalue += 1
	for node in tree.postorder_node_iter():
		node.yvalue = get_yvalue(node)
		node.xvalue = node.distance_from_root()
	root = tree.seed_node
	for node in tree.postorder_node_iter():
		node.ep = epitope_distance(translate(node.seq), translate(root.seq))
		node.ne = nonepitope_distance(translate(node.seq), translate(root.seq))
		node.rb = receptor_binding_distance(translate(node.seq), translate(root.seq))
	for node in tree.postorder_node_iter():
		node.trunk_count = 0
		node.trunk = False

def translate_all(tree):
	for node in tree.postorder_node_iter():
		node.aa_seq = translate(node.seq)

def unique_date(tree):
	leaf_count = 0
	for node in tree.postorder_node_iter():
		if node.is_leaf():
			# attach index to a leaf, to allow for array indexing later
			node.tip_index = leaf_count
			# modify date by a tiny amount << than a day to ensure uniqueness
			node.num_date = numerical_date(node.date) + 1e-7*node.tip_index
			leaf_count+=1
		else: # internal node preceed the oldest child
			node.num_date = min([c.num_date for c in node.child_nodes()])

def define_trunk(tree):
	"""Trace current lineages backward to define trunk"""

	# Find most recent tip
	dates = []
	for node in tree.postorder_node_iter():
		if node.is_leaf():
			dates.append(node.date)
	most_recent_date = string_to_date(sorted(dates)[-1])

	# Mark ancestry of recent tips
	number_recent = 0
	for node in tree.postorder_node_iter():
		if node.is_leaf():
			diff = year_difference(string_to_date(node.date), most_recent_date)
			if (diff < 1):
				number_recent += 1
				parent = node.parent_node
				while (parent != None):
					parent.trunk_count += 1
					parent = parent.parent_node

	# Mark trunk nodes
	for node in tree.postorder_node_iter():
		if node.trunk_count == number_recent:
			node.trunk = True;


def main(tree_fname = 'data/tree_ancestral.json', virus_fname='data/virus_clean.json'):

	print "--- Tree refine at " + time.strftime("%H:%M:%S") + " ---"

	viruses = read_json(virus_fname)
	tree =  json_to_dendropy(read_json(tree_fname))
	print "Remove outgroup"
	remove_outgroup(tree)
	print "Remove outlier branches"
	reduce(tree)
	print "Collapse internal nodes"
	collapse(tree)
	print "Ladderize tree"
	ladderize(tree)
	print "Append node attributes"
	add_virus_attributes(viruses, tree)
	add_node_attributes(tree)
	print "translate nucleotide sequences"
	translate_all(tree)
	print "enumerate leaves of ladderized tree and calculate unique numerical date"
	unique_date(tree)
	print "Define trunk"
	define_trunk(tree)
	out_fname = "data/tree_refine.json"
	write_json(dendropy_to_json(tree.seed_node), out_fname)
	return out_fname

if __name__ == "__main__":
	main()
