import time, os
from io_util import *
import tree_auspice

def main():
	"""Prep for auspice"""
	
	print "--- Start processing at " + time.strftime("%H:%M:%S") + " ---"	
    
    # Streamline tree for auspice
	tree_auspice.main()			
	
	# Write out metadata
	meta = {"updated": time.strftime("%d %b %Y")}
	write_json(meta, "auspice/meta.json")

if __name__ == "__main__":
    main()
