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
	fname = tree_refine.main(fname)			# Clean tree, reads viruses in fname + raxml files
	fname = tree_streamline.main(fname)		# Streamline tree for auspice

	# Write out metadata
	print "Writing out metadata"
	meta = {"updated": time.strftime("%d %b %Y")}
	fname_meta = "../auspice/data/meta.json"
	write_json(meta, fname_meta)

if __name__ == "__main__":
	main()
