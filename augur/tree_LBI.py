import dendropy, datetime
import numpy as np

def from_json(json):
	'''
	read a json dictionary and make a dendropy tree from it. 
	'''
	tree = dendropy.Tree()
	tree.get_from_string(';', 'newick')
	root = tree.leaf_nodes()[0]
	from_sub_json(json, root)
	root.edge_length=1.0
	return tree

def from_sub_json(json, node):
	'''
	recursively calls itself for all children of node and 
	builds up the tree. entries in json are added as node attributes
	'''
	for attr,val in json.iteritems():
		if attr=='children':
			for sub_json in val:
				child_node = dendropy.Node()
				from_sub_json(sub_json, child_node)
				node.add_child(child_node, edge_length = child_node.xvalue - node.xvalue)
		elif attr=='date':
			y,m,d = map(int, val.split('-')) 
			node.date = datetime.date(year=y, month=m, day=d).toordinal()/365.0
		else:
			try:
				node.__setattr__(attr, float(val))
			except:
				node.__setattr__(attr, val)
	if len(node.child_nodes())==0:
		node.taxon= True
	else:
		node.taxon = False

def get_average_T2(tree, dt):
	'''
	calculates the average pairwise distance between nodes from a time interval dt
	repeats this for a bunch of dt intervals, returns the average
	tree  -- dendropy tree
	dt	  -- time interval
	'''

	leaves = tree.leaf_nodes()
	leaves.sort(key=lambda x:x.date, reverse = True)
	upper_cutoff = leaves[0].date
	ndt = 4
	T2 = 0
	ntrees = 0
	for ii in xrange(ndt):
		leaf_set = set([leaf for leaf in leaves 
				if leaf.date<upper_cutoff-ii*dt 
				and leaf.date>=upper_cutoff-(ii+1)*dt])
		if len(leaf_set)>10:
			tmp_T2= avg_T2_of_leaf_set(tree, leaf_set) 
			T2 += tmp_T2
			ntrees += 1
	return T2/ntrees

def avg_T2_of_leaf_set(tree, leaf_set):
	'''
	walk through the nodes in postorder and for each edge increment the T2
	count by the edge length multiplied by the number of node pairs separated
	by this split
	tree	 -- dendropy tree containing the distances
	leaf_set -- subset of the leaves of tree whose pairwise distance is returned
	'''
	T2 = 0
	nleaves = len(leaf_set)
	for node in tree.postorder_node_iter():
		# get number of leaves in leaf_set that are descendants of node
		ndescendants = len(leaf_set.intersection(node.leaf_nodes()))
		# multiply edge_length by the number of pairs across the split
		T2 += node.edge_length*ndescendants*(nleaves-ndescendants)
	return 2*T2/nleaves/(nleaves-1)

def calc_LBI(tree, tau):   
	'''
	traverses the tree in postorder and preorder to calculate the
	up and downstream tree length exponentially weighted by distance.
	then adds them as LBI
	tree -- dendropy tree for whose node the LBI is being computed
	tau	 -- the memory time scale of tree length along branches of the tree
	'''
	import numpy as np
	# traverse the tree in postorder (children first) to calculate msg to parents
	for node in tree.postorder_node_iter():
		node.down_polarizer = 0
		node.up_polarizer = 0
		for child in node.child_nodes():
			node.up_polarizer += child.up_polarizer
		bl =  node.edge_length/tau
		node.up_polarizer *= np.exp(-bl)
		node.up_polarizer += tau*(1-np.exp(-bl))

	# traverse the tree in preorder (parents first) to calculate msg to children
	for node in tree.preorder_internal_node_iter():
		for child1 in node.child_nodes():
			child1.down_polarizer = node.down_polarizer
			for child2 in node.child_nodes():
				if child1!=child2:
					child.down_polarizer += child.up_polarizer
			
			bl =  child1.edge_length/tau
			child.down_polarizer *= np.exp(-bl)
			node.down_polarizer += tau*(1-np.exp(-bl))

	# go over all nodes and calculate the LBI (can be done in any order)
	for node in tree.postorder_node_iter():
		node.LBI = node.down_polarizer
		for child in node.child_nodes():
			node.LBI += child.up_polarizer
	#done

def test():
	from io_util import read_json 
	tree = from_json(read_json('tree.json'))
	print "calculate local branching index"
	T2 = get_average_T2(tree, 1)
	tau =  T2*2**-4
	print "avg pairwise distance:", T2
	print "memory time scale:", tau
	calc_LBI(tree, tau = tau)
	return tree


if __name__=='__main__':
	tree = test()


