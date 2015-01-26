import time, os
import virus_download, virus_filter, virus_align, virus_clean
import fitness_epitope, fitness_nonepitope
import tree_infer, tree_refine, tree_LBI, tree_streamline

def main():
	"""Run full pipeline"""

	print "--- Start processing at " + time.strftime("%H:%M:%S") + " ---"

#	virus_download.main()		# Download from GISAID
	virus_filter.main()			# Filter sequences
	virus_align.main()			# Align sequences
	virus_clean.main()			# Clean sequences
#	fitness_epitope.main()		# Calculate epitope fitness
#	fitness_nonepitope.main()	# Calculate non-epitope fitness
	tree_infer.main()			# Make tree
	tree_refine.main()			# Clean tree
	tree_LBI.main()				# Calculate LBI across tree
	tree_streamline.main()		# Streamline tree for auspice

if __name__ == "__main__":
	main()
