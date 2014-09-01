import schedule, time, os
import virus_ingest, virus_filter, virus_align, virus_clean
import tree_infer, tree_clean, tree_streamline

def pipeline():
	"""Run full pipeline"""
	virus_ingest.main()			# Ingest sequences
	virus_filter.main()			# Filter sequences
	virus_align.main()			# Align sequences
	virus_clean.main()			# Clean sequences	
	tree_infer.main()			# Make tree
	tree_clean.main()			# Clean tree	
	tree_streamline.main()		# Streamline tree

def main():
	"""Run every day"""

	if not os.path.exists("log"):
		os.makedirs("log")
    
	if not os.path.exists("data"):
		os.makedirs("data")    

	time.sleep(10)
	pipeline()
	schedule.every().day.do(pipeline)
	while True:
		schedule.run_pending()
		time.sleep(1)	

if __name__ == "__main__":
    main()
