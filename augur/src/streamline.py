import time, os, shutil
from io_util import *
from tree_util import *

def main(in_fname='data/tree_frequencies.json'):
	"""Prep tree for auspice, stripping sequence data"""

	print "--- Streamline at " + time.strftime("%H:%M:%S") + " ---"

	# Streamline tree for auspice
	print "Writing streamlined tree"
	tree = read_json(in_fname)
	for node in all_descendants(tree):
		node.pop("seq", None)
		node.pop("clade", None)

	out_fname_tree = "../auspice/data/tree.json"
	write_json(tree, out_fname_tree, indent=None)
	try:
		read_json(out_fname_tree)
	except:
		print "Read failed, rewriting with indents"	
		write_json(tree, out_fname_tree, indent=1)
		
	# Include genotype frequencies
	shutil.copy2("data/genotype_frequencies.json", "../auspice/data/frequencies.json")

if __name__ == "__main__":
	main()
