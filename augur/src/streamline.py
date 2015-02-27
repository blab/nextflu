import time, os, shutil
from io_util import *
from tree_util import *

def main(tree_json):
	"""Prep tree for auspice, stripping sequence data"""

	print "--- Streamline at " + time.strftime("%H:%M:%S") + " ---"

	# Move sequence data to separate file
	print "Writing sequences"
	elems = []
	for node in all_descendants(tree_json):
		elem = {}
		if 'clade' in node:
			elem['clade'] = node['clade']
		if 'aa_seq' in node:
			elem['aa_seq'] = node['aa_seq']			
		elems.append(elem)
	write_json(elems, "../auspice/data/sequences.json", indent=None)

	# Streamline tree for auspice
	print "Writing streamlined tree"
	for node in all_descendants(tree_json):
		node.pop("seq", None)
		node.pop("aa_seq", None)
		node.pop("logit_freq", None)

	out_fname_tree = "../auspice/data/tree.json"
	write_json(tree_json, out_fname_tree, indent=None)
	try:
		read_json(out_fname_tree)
	except:
		print "Read failed, rewriting with indents"	
		write_json(tree_json, out_fname_tree, indent=1)
		
	# Include genotype frequencies
	shutil.copy2("data/genotype_frequencies.json", "../auspice/data/frequencies.json")

if __name__ == "__main__":
	main()
