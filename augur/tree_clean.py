# clean, reroot, ladderize newick tree
# output to tree.json

import os, re, time
import dendropy
from io_util import *

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

def to_json(node):
	json = {}
	if hasattr(node, 'clade'):
		json['clade'] = node.clade
	if node.taxon:
		json['strain'] = str(node.taxon).replace("'", '')
	if hasattr(node, 'yvalue'):
		json['yvalue'] = round(node.yvalue, 5)
	if hasattr(node, 'xvalue'):
		json['xvalue'] = round(node.xvalue, 5)	
	if hasattr(node, 'date'):
		json['date'] = node.date
	if hasattr(node, 'seq'):
		json['seq'] = node.seq
	if node.child_nodes():
		json["children"] = []
		for ch in node.child_nodes():
			json["children"].append(to_json(ch))
	return json
	
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
	outgroup_node = tree.find_node_with_taxon_label(OUTGROUP)	
	if outgroup_node:
#		tree.to_outgroup_position(outgroup_node, update_splits=False)
		tree.prune_subtree(outgroup_node)
		
def collapse(tree):
	"""Collapse short edges to polytomies"""
	for edge in tree.postorder_edge_iter():
		if edge.length < 0.00001 and edge.is_internal():
			edge.collapse()
			
def reduce(tree):
	"""Remove outlier tips"""
	for node in tree.postorder_node_iter():
		if node.edge_length > 0.04 and node.is_leaf():
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

def add_node_attributes(tree):
	"""Add clade, xvalue and yvalue attributes to all nodes in tree"""
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
	"""Add date and seq attributes to all tips in tree"""
	strain_to_date = {}
	strain_to_seq = {}
	for v in viruses:
		strain_to_date[v['strain']] = v['date']
		strain_to_seq[v['strain']] = v['seq']
	for node in tree.postorder_node_iter():
		strain = str(node.taxon).replace("'", '')
		if strain_to_date.has_key(strain):
			node.date = strain_to_date[strain]
		if strain_to_seq.has_key(strain):
			node.seq = strain_to_seq[strain]
																		
def main():

	print "--- Tree clean at " + time.strftime("%H:%M:%S") + " ---"
		
	viruses = read_json('data/virus_clean.json')
	tree = crossref_import('data/tree_branches.newick', 'data/tree_states.newick', 'data/states.txt')
	print "Remove outgroup"
	remove_outgroup(tree)
	print "Remove outlier branches"	
	reduce(tree)
	print "Collapse internal nodes"		
	collapse(tree)	
	print "Ladderize tree"	
	ladderize(tree)
	add_node_attributes(tree)
	add_virus_attributes(viruses, tree)

	write_json(to_json(tree.seed_node), "data/tree_clean.json")
	
if __name__ == "__main__":
    main()