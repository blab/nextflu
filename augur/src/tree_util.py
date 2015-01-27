import dendropy
from io_util import *

def tip_descendants(node):
	"""Take node, ie. dict, and return a flattened list of all tips descending from this node"""
	if 'children' in node:
		for child in node['children']:
			for desc in tip_descendants(child):
				yield desc
	else:
		yield node

def all_descendants(node):
	"""Take node, ie. dict, and return a flattened list of all nodes descending from this node"""
	yield node
	if 'children' in node:
		for child in node['children']:
			for desc in all_descendants(child):
				yield desc

def get_dates(node):
	"""Return ordered list of dates of descendants of a node"""
	return sorted([n['date'] for n in tip_descendants(node)])

def dendropy_to_json(node):
	json = {}
	if hasattr(node, 'clade'):
		json['clade'] = node.clade
	if node.taxon:
		json['strain'] = str(node.taxon).replace("'", '')
	if hasattr(node, 'xvalue'):
		json['xvalue'] = round(node.xvalue, 5)
	if hasattr(node, 'yvalue'):
		json['yvalue'] = round(node.yvalue, 5)
	if hasattr(node, 'ep'):
		json['ep'] = node.ep
	if hasattr(node, 'ne'):
		json['ne'] = node.ne
	if hasattr(node, 'rb'):
		json['rb'] = node.rb
	if hasattr(node, 'date'):
		json['date'] = node.date
	if hasattr(node, 'seq'):
		json['seq'] = node.seq
	if hasattr(node, 'LBI'):
		json['LBI'] = round(node.LBI,5)
	if node.child_nodes():
		json["children"] = []
		for ch in node.child_nodes():
			json["children"].append(dendropy_to_json(ch))
	return json

def BioPhylo_to_json(node):
	json = {}
	if hasattr(node, 'clade'):
		json['clade'] = node.clade
	if node.name:
		json['strain'] = str(node.name).replace("'", '')
	if hasattr(node, 'branch_length'):
		json['branch_length'] = round(node.branch_length, 5)
	if hasattr(node, 'xvalue'):
		json['xvalue'] = round(node.xvalue, 5)
	if hasattr(node, 'yvalue'):
		json['yvalue'] = round(node.yvalue, 5)
	if hasattr(node, 'ep'):
		json['ep'] = node.ep
	if hasattr(node, 'ne'):
		json['ne'] = node.ne
	if hasattr(node, 'rb'):
		json['rb'] = node.rb
	if hasattr(node, 'date'):
		json['date'] = node.date
	if hasattr(node, 'seq'):
		json['seq'] = str(node.seq)
	if hasattr(node, 'LBI'):
		json['LBI'] = round(node.LBI,5)
	if len(node.clades):
		json["children"] = []
		for ch in node.clades:
			json["children"].append(BioPhylo_to_json(ch))
	return json


def json_to_dendropy(json):
	'''
	read a json dictionary and make a dendropy tree from it.
	'''
	tree = dendropy.Tree()
	tree.get_from_string(';', 'newick')
	root = tree.leaf_nodes()[0]
	json_to_dendropy_sub(json, root)
	root.edge_length=0.0
	return tree

def json_to_dendropy_sub(json, node):
	'''
	recursively calls itself for all children of node and
	builds up the tree. entries in json are added as node attributes
	'''
	for attr,val in json.iteritems():
		if attr=='children':
			for sub_json in val:
				child_node = dendropy.Node()
				json_to_dendropy_sub(sub_json, child_node)
				if hasattr(child_node, 'xvalue'):
					node.add_child(child_node, edge_length = child_node.xvalue - node.xvalue)
				elif hasattr(child_node, 'branch_length'):
					node.add_child(child_node, edge_length = child_node.branch_length)
				else:
					node.add_child(child_node, edge_length = 1.0)
		else:
			try:
				node.__setattr__(attr, float(val))
			except:
				node.__setattr__(attr, val)
	if len(node.child_nodes())==0:
		node.taxon = json['strain']

def main():

	tree = read_json('tree.json')

#	print "Whole tree"
#	for tip in descendants(tree):
#		print tip['date']

#	node = tree['children'][0]

#	dates = get_dates(tree)
#	print dates

	for node in all_descendants(tree):
		dates = get_dates(node)
		print str(node['clade']) + ": " + str(len(dates))

if __name__ == "__main__":
	main()