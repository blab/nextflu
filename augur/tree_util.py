# tree class
# each node is a dict in which the key 'children' is used to reference more dicts

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