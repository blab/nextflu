import schedule, time, os
import ingest, filter, align, clean, tree, streamline, upload

def pipeline():
	"""Run full pipeline"""
	ingest.main()		# Ingest sequences
	filter.main()		# Filter sequences
	align.main()		# Align sequences
	clean.main()		# Clean sequences	
	tree.main()			# Make tree
	streamline.main()	# Streamline tree

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
