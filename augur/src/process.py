import time, os, sys
import virus_download, virus_filter, virus_align, virus_clean
import tree_infer, tree_ancestral, tree_refine, tree_streamline
from io_util import *

def main(years_back=3, viruses_per_month=50):
	"""Run full pipeline"""

	print "--- Start processing at " + time.strftime("%H:%M:%S") + " ---"

	# Run pipeline
#	virus_download.main()		# Download from GISAID
	virus_fname = 'data/gisaid_epiflu_sequence.fasta'

	# Filter sequences
	virus_fname = virus_filter.main(virus_fname, years_back=years_back, viruses_per_month=viruses_per_month)

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

	# Streamline tree for auspice
	tree_fname = tree_streamline.main(tree_fname)

	# Write out metadata
	print "Writing out metadata"
	meta = {"updated": time.strftime("%d %b %Y")}
	meta_fname = "../auspice/data/meta.json"
	write_json(meta, meta_fname)

if __name__ == "__main__":
	years_back = 3
	viruses_per_month = 50
	if (len(sys.argv) > 1):
		years_back = int(sys.argv[1])
	if (len(sys.argv) > 2):
		viruses_per_month = int(sys.argv[2])
	main(years_back=years_back, viruses_per_month=viruses_per_month)
