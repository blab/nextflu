import time, os
from io_util import *
from tree_util import *

def main():
	"""Prep tree for auspice, stripping sequence data"""

	print "--- Streamline at " + time.strftime("%H:%M:%S") + " ---"

	# Streamline tree for auspice
	print "Writing streamlined tree to tree_streamline.json"
	tree = read_json('data/tree_LBI.json')
	for node in all_descendants(tree):
		node.pop("seq", None)
		node.pop("clade", None)
		node.pop("LBI", None)
	write_json(tree, "data/tree_streamline.json", indent=0)

	# Write out metadata
	print "Writing out metadata to meta.json"
	meta = {"updated": time.strftime("%d %b %Y")}
	write_json(meta, "data/meta.json")

if __name__ == "__main__":
	main()
