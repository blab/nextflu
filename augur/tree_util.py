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
	if hasattr(node, 'distance_ep'):
		json['distance_ep'] = node.distance_ep
	if hasattr(node, 'distance_ne'):
		json['distance_ne'] = node.distance_ne
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
	
def json_to_dendropy(json):
    '''
    read a json dictionary and make a dendropy tree from it. 
    '''
    tree = dendropy.Tree()
    tree.get_from_string(';', 'newick')
    root = tree.leaf_nodes()[0]
    json_to_dendropy_sub(json, root)
    root.edge_length=0.01
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
                node.add_child(child_node, edge_length = child_node.xvalue - node.xvalue)
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