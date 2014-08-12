# create phylogeny from alignment
# raxml will complain about identical sequences, but still includes them all in the resulting tree
# writes out newick.tree file

import os, re, json, time, seqmagick
from ete2 import Tree
from share import *

RAXML = 'raxmlHPC-PTHREADS-AVX'

def to_json(node):
	node.name = node.name.replace("'", '').replace('NoName','')
	json = { 
		"name": node.name, 
		"branch": str(node.dist),
		"layout": str(node.layout),
		"distance": str(node.distance),		
	}
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
	os.system(RAXML + " -T 6 -s temp.phyx -n out -c 25 -f d -m GTRCAT -p 344312987")
	os.rename('RAxML_bestTree.out', 'tree.newick')

	ete_tree = Tree("tree.newick", format=5)
	reroot(ete_tree)
	ete_tree.ladderize(direction=1)	
	
	layout = 0
	for node in ete_tree.traverse("postorder"):
		if node.is_leaf():
			layout += 1
			node.add_feature("layout", layout)
					
	for node in ete_tree.traverse("postorder"):
		node.add_feature("layout", get_layout(node))
		node.add_feature("distance", get_distance(node))			
					
	write_json(to_json(ete_tree), "tree.json")

	cleanup()
	
if __name__ == "__main__":
    main()