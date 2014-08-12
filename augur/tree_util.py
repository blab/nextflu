# tree class
# each node is a dict in which the key 'children' is used to reference more dicts

import os, re, json

def descendants(node):
	"""Take node, ie. dict, and return a flattened list of all tips descending from this node"""
	if 'children' in node:
		for child in node['children']:
			for desc in descendants(child):
				yield desc
	else:
		yield node
			
def main():

	try:
		handle = open('tree.json', 'r')  
	except IOError:
		pass
	else:	
 		tree = json.load(handle)
 		handle.close()

	print "Whole tree"
	for tip in descendants(tree):
		print tip
		
	print "Internal node"		
	node = tree['children'][0]
	for tip in descendants(node):
		print tip
		
if __name__ == "__main__":
    main()