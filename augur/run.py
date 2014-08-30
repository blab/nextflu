import schedule, time
import ingest, filter, align, clean, tree, streamline, upload

def pipeline():
	"""Run full pipeline"""
#	ingest.main()		# Ingest sequences
	filter.main()		# Filter sequences
	align.main()		# Align sequences
	clean.main()		# Clean sequences	
	tree.main()			# Make tree
	streamline.main()	# Streamline tree
	upload.main()		# Upload tree

def main():
	"""Run every day"""

	pipeline()
	schedule.every().minute.do(pipeline)
	while True:
		schedule.run_pending()
		time.sleep(1)	

if __name__ == "__main__":
    main()
