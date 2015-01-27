import time, os
import virus_download, virus_filter, virus_align, virus_clean
import tree_infer, tree_refine, tree_LBI, tree_streamline
from io_util import *

def main():
	"""Run full pipeline"""

	print "--- Start processing at " + time.strftime("%H:%M:%S") + " ---"

	# Run pipeline
#	virus_download.main()		# Download from GISAID
	fname = 'data/gisaid_epiflu_sequence.fasta'
	fname = virus_filter.main(fname)		# Filter sequences
	fname = virus_align.main(fname)			# Align sequences
	fname = virus_clean.main(fname)			# Clean sequences
	tree_infer.main(fname)					# Make tree, creates raxml files
	# infer ancestral states using the cleaned viruses and the raxml tree
	tree_fname = tree_infer_ancestral(tree_fname = 'data/raxml_branches.newick', 
	                          		  virus_fname = fname)
	# Clean tree, reads viruses in fname + raxml files
	tree_fname = tree_refine.main(tree_fname=tree_fname, 
	                              	virus_fname = fname)
	tree_fname = tree_streamline.main(tree_fname)		# Streamline tree for auspice

	# Write out metadata
	print "Writing out metadata"
	meta = {"updated": time.strftime("%d %b %Y")}
	fname_meta = "../auspice/data/meta.json"
	write_json(meta, fname_meta)

if __name__ == "__main__":
	main()
