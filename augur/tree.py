# create phylogeny from alignment
# raxml will complain about identical sequences, but still includes them all in the resulting tree
# writes out newick.tree file

import os, time, seqmagick
from ete2 import Tree
from io_util import *

OUTGROUP = 'A/Beijing/32/1992'

def to_json(node):
	node.name = node.name.replace("'", '')
	json = {}
	if hasattr(node, 'clade'):
		json['clade'] = node.clade
	if hasattr(node, 'branch'):
		json['branch'] = round(node.branch, 5)
	if hasattr(node, 'layout'):
		json['layout'] = round(node.layout, 5)
	if hasattr(node, 'distance'):
		json['distance'] = round(node.distance, 5)
	if node.name != 'NoName':
		json['strain'] = node.name
	if hasattr(node, 'date'):
		json['date'] = node.date
	if hasattr(node, 'seq'):
		json['seq'] = node.seq
	if node.children:
		json["children"] = []
		for ch in node.children:
			json["children"].append(to_json(ch))
	return json
	
def get_layout(node):
	"""Return y location based on recursive mean of daughter locations"""	
	if hasattr(node, 'layout'):
		return node.layout	
	if node.children:
		mean = 0
		for ch in node.children:
			mean += get_layout(ch)
		return mean / float(len(node.children))
		
def get_distance(node):
	"""Return x location based on total distance from root"""
	root = node.get_tree_root()
	return node.get_distance(root)

def reroot(tree):
	"""Reroot tree to outgroup"""
	for node in tree:
		if node.name == OUTGROUP:
			outgroup = node
			break
	tree.set_outgroup(outgroup)

def add_node_attributes(tree):
	"""Add clade, distance and layout attributes to all nodes in tree"""
	clade = 0
	layout = 0
	for node in tree.traverse("postorder"):
		node.add_feature("clade", clade)
		clade += 1
		if node.is_leaf():
			node.add_feature("layout", layout)
			layout += 1

	for node in tree.traverse("postorder"):
		node.add_feature("layout", get_layout(node))
		node.add_feature("distance", get_distance(node))

def add_virus_attributes(viruses, tree):
	"""Add date and seq attributes to all tips in tree"""
	strain_to_date = {}
	strain_to_seq = {}
	for v in viruses:
		strain_to_date[v['strain']] = v['date']
		strain_to_seq[v['strain']] = v['seq']
	for node in tree.traverse("postorder"):
		if strain_to_date.has_key(node.name):
			node.add_feature("date", strain_to_date[node.name])
		if strain_to_seq.has_key(node.name):
			node.add_feature("seq", strain_to_seq[node.name])

def cleanup():
	try:
		os.remove('temp.fasta')
	except OSError:
		pass
	try:
		os.remove('temp.phyx')
	except OSError:
		pass
	try:
		os.remove('temp.phyx.reduced')
	except OSError:
		pass		
	try:
		os.remove('RAxML_info.out')
	except OSError:
		pass	
	try:
		os.remove('RAxML_log.out')
	except OSError:
		pass		
	try:
		os.remove('RAxML_parsimonyTree.out')
	except OSError:
		pass
	try:
		os.remove('RAxML_result.out')
	except OSError:
		pass
									
def main():

	print "--- Tree at " + time.strftime("%H:%M:%S") + " ---"
	
	cleanup()	
	
	viruses = read_json('virus_clean.json')
	write_fasta(viruses, 'temp.fasta')
	os.system("seqmagick convert temp.fasta temp.phyx")
	os.system("raxml -T 6 -s temp.phyx -n out -c 25 -f d -m GTRCAT -p 344312987")
	os.rename('RAxML_bestTree.out', 'tree.newick')

	tree = Tree("tree.newick", format=5)
	reroot(tree)
	tree.ladderize(direction=1)
	add_node_attributes(tree)
	add_virus_attributes(viruses, tree)
	
	write_json(to_json(tree), "tree.json")

	cleanup()
	
if __name__ == "__main__":
    main()