import time, os
from io_util import *
from tree_util import *

def main(in_fname=None):
	"""Prep tree for auspice, stripping sequence data"""

	print "--- Streamline at " + time.strftime("%H:%M:%S") + " ---"

	# Streamline tree for auspice
	print "Writing streamlined tree"
	if in_fname is None: in_fname='data/tree_refine.json'
	tree = read_json(in_fname)
	for node in all_descendants(tree):
		node.pop("seq", None)
		node.pop("clade", None)

	out_fname_tree = "../auspice/data/tree.json"
	write_json(tree, out_fname_tree, indent=0)

if __name__ == "__main__":
	main()
