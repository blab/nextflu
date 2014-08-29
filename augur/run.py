import ingest, filter, align, clean, tree, streamline

def main():
	"""Run full pipeline"""
	ingest.main()		# Ingest sequences
	filter.main()		# Filter sequences
	align.main()		# Align sequences
	clean.main()		# Clean sequences	
	tree.main()			# Make tree
	streamline.main()	# Streamline tree			

if __name__ == "__main__":
    main()
