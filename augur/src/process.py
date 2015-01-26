import time, os
import virus_download, virus_filter, virus_align, virus_clean
import fitness_epitope, fitness_nonepitope
import tree_infer, tree_refine, tree_LBI, tree_streamline

def main():
	"""Run full pipeline"""

	print "--- Start processing at " + time.strftime("%H:%M:%S") + " ---"

#	virus_download.main()		# Download from GISAID
	fname = 'data/gisaid_epiflu_sequence.fasta'
	fname = virus_filter.main(fname)		# Filter sequences
	fname = virus_align.main(fname)			# Align sequences
	fname = virus_clean.main(fname)			# Clean sequences
#	fname = fitness_epitope.main(fname)		# Calculate epitope fitness
#	fname = fitness_nonepitope.main(fname)	# Calculate non-epitope fitness
	tree_infer.main(fname)					# Make tree, creates raxml files
	fname = tree_refine.main(fname)			# Clean tree, reads viruses in fname + raxml files
	fname_tree, fname_meta = tree_streamline.main(fname)		# Streamline tree for auspice

if __name__ == "__main__":
	main()
