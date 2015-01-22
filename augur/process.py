import time, os
import virus_filter, virus_align, virus_clean, virus_reduce
import fitness_epitope
import tree_infer, tree_clean

def main():
	"""Run full pipeline"""
	
	print "--- Start processing at " + time.strftime("%H:%M:%S") + " ---"	
    
	virus_filter.main()			# Filter sequences
	virus_align.main()			# Align sequences
	virus_clean.main()			# Clean sequences
	fitness_epitope.main()		# Calculate epitope fitness
	virus_reduce.main()			# Reduce sequences	
	tree_infer.main()			# Make tree
	tree_clean.main()			# Clean tree	

if __name__ == "__main__":
    main()
