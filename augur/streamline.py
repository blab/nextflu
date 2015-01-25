import time, os
from io_util import *
from tree_util import *

def main():
	"""Prep for auspice"""
	
	print "--- Streamline at " + time.strftime("%H:%M:%S") + " ---"	
    
    # Streamline tree for auspice
	tree = read_json('data/tree_LBI.json')
	for node in all_descendants(tree):
		node.pop("seq", None)
		node.pop("clade", None)
	write_json(tree, "auspice/tree.json")		
    	
	# Write out metadata
	meta = {"updated": time.strftime("%d %b %Y")}
	write_json(meta, "auspice/meta.json")

if __name__ == "__main__":
    main()
