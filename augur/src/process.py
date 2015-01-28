import time, os
import virus_download, virus_filter, virus_align, virus_clean
import tree_infer, tree_ancestral, tree_refine, tree_streamline
from io_util import *

def main():
	"""Run full pipeline"""

	print "--- Start processing at " + time.strftime("%H:%M:%S") + " ---"

	# Run pipeline
#	virus_download.main()		# Download from GISAID
	virus_fname = 'data/gisaid_epiflu_sequence.fasta'

	# Filter sequences
	virus_fname = virus_filter.main(virus_fname)

	# Align sequences
	virus_fname = virus_align.main(virus_fname)

	# Clean sequences
	virus_fname = virus_clean.main(virus_fname)

	# Make tree, creates raxml files
	tree_fname = tree_infer.main(virus_fname)

	# infer ancestral states using the cleaned viruses and the raxml tree
	tree_fname = tree_ancestral.main(tree_fname=tree_fname, virus_fname=virus_fname)

	# Clean tree, reads viruses in fname + raxml files
	tree_fname = tree_refine.main(tree_fname=tree_fname, virus_fname = virus_fname)

	# Calculate virus fitness based on mutations tolerance from Thyagarajan and Bloom
	tree_fname = fitness_tolerance.main(tree_fname)

	# Streamline tree for auspice
	tree_fname = tree_streamline.main(tree_fname, tree=True)

	# Write out metadata
	print "Writing out metadata"
	meta = {"updated": time.strftime("%d %b %Y")}
	meta_fname = "../auspice/data/meta.json"
	write_json(meta, meta_fname)

if __name__ == "__main__":
	main()
