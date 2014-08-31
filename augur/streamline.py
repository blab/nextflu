# streamline tree for visualization in auspice

import time, os
from io_util import *
from tree_util import *

def main():
	print "--- Streamline at " + time.strftime("%H:%M:%S") + " ---"

	tree = read_json('tree.json')

	for node in all_descendants(tree):
		node.pop("seq", None)
		node.pop("clade", None)
	
	if not os.path.exists("data"):
		os.makedirs("data")	
	write_json(tree, "data/tree.json")		

if __name__ == "__main__":
    main()
