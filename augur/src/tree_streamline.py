import time, os
from io_util import *
from tree_util import *

def main(in_fname=None):
	"""Prep tree for auspice, stripping sequence data"""

	print "--- Streamline at " + time.strftime("%H:%M:%S") + " ---"

	# Streamline tree for auspice
	print "Writing streamlined tree to tree_streamline.json"
	if in_fname is None: in_fname='data/tree_refine.json'
	tree = read_json(in_fname)
	for node in all_descendants(tree):
		node.pop("seq", None)
		node.pop("clade", None)

	out_fname_tree = "data/tree_streamline.json"
	write_json(tree, out_fname_tree, indent=0)

	# Write out metadata
	print "Writing out metadata to meta.json"
	meta = {"updated": time.strftime("%d %b %Y")}
	out_fname_meta = "data/meta.json"
	write_json(meta, out_fname_meta)
	return out_fname_tree, out_fname_meta

if __name__ == "__main__":
	main()
