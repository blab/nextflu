import dendropy, time
import numpy as np
from seq_util import epitope_distance
from seq_util import nonepitope_distance
from io_util import read_json
from io_util import write_json
from tree_util import json_to_dendropy
from tree_util import dendropy_to_json

def calc_epitope_distance(tree, attr='ep', ref = None):
	'''
	calculates the distance at epitope sites of any tree node  to ref
	tree   --   dendropy tree
	attr   --   the attribute name used to save the result
	'''
	if not hasattr(tree, "epitope_distance_assigned") or tree.epitope_distance_assigned==False:
		if ref == None:
			ref = tree.seed_node
		for node in tree.postorder_node_iter():
			node.__setattr__(attr, epitope_distance(node.seq, ref.seq))
		tree.epitope_distance_assigned=True

def calc_nonepitope_distance(tree, attr='ne', ref = None):
	'''
	calculates the distance at nonepitope sites of any tree node to ref
	tree   --   dendropy tree
	attr   --   the attribute name used to save the result
	'''
	if not hasattr(tree, "nonepitope_distance_assigned") or tree.nonepitope_distance_assigned==False:
		if ref == None:
			ref = tree.seed_node
		for node in tree.postorder_node_iter():
			node.__setattr__(attr, nonepitope_distance(node.seq, ref.seq))
		tree.nonepitope_distance_assigned=True

def calc_LBI(tree, attr = 'lbi', tau=0.0005, transform = lambda x:x):
	'''
	traverses the tree in postorder and preorder to calculate the
	up and downstream tree length exponentially weighted by distance.
	then adds them as LBI
	tree -- dendropy tree for whose node the LBI is being computed
	attr	 -- the attribute name used to store the result
	'''
	# traverse the tree in postorder (children first) to calculate msg to parents
	for node in tree.postorder_node_iter():
		node.down_polarizer = 0
		node.up_polarizer = 0
		for child in node.child_nodes():
			node.up_polarizer += child.up_polarizer
		bl =  node.edge_length/tau
		node.up_polarizer *= np.exp(-bl)
		if node.alive: node.up_polarizer += tau*(1-np.exp(-bl))

	# traverse the tree in preorder (parents first) to calculate msg to children
	for node in tree.preorder_internal_node_iter():
		for child1 in node.child_nodes():
			child1.down_polarizer = node.down_polarizer
			for child2 in node.child_nodes():
				if child1!=child2:
					child1.down_polarizer += child2.up_polarizer

			bl =  child1.edge_length/tau
			child1.down_polarizer *= np.exp(-bl)
			if child1.alive: child1.down_polarizer += tau*(1-np.exp(-bl))

	# go over all nodes and calculate the LBI (can be done in any order)
	for node in tree.postorder_node_iter():
		tmp_LBI = node.down_polarizer
		for child in node.child_nodes():
			tmp_LBI += child.up_polarizer
		node.__setattr__(attr, transform(tmp_LBI))


def main(tree_fname = 'data/tree_refine.json'):

	print "--- Testing predictor evaluations ---"
	tree =  json_to_dendropy(read_json(tree_fname))

	print "Calculating epitope distances"
	calc_epitope_distance(tree)

	print "Calculating nonepitope distances"
	calc_nonepitope_distance(tree)

	print "Calculating LBI"
#	calc_LBI(tree)

	print "Writing decorated tree"
	out_fname = "data/tree_predictors.json"
	write_json(dendropy_to_json(tree.seed_node), out_fname)
	return out_fname

if __name__=='__main__':
	main()
