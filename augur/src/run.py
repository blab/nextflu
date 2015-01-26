import time, os, threading, argparse
import ingest, process, sync

def main():
	"""Ingest, process, sync"""

	if not os.path.exists("data"):
		os.makedirs("data")

	parser = argparse.ArgumentParser(description='Ingest virus sequences')
	parser.add_argument('--headless', action='store_true', help='Run firefox in headless state (requires xvfb and x11vnc)', required=False)
	parser.add_argument('--clock', action='store_true', help='Run ntpdate to fix date', required=False)
	args = vars(parser.parse_args())
	headless = args['headless']
	clock = args['clock']

	ingest.main([headless])
	process.main()
	sync.main([clock])

if __name__ == "__main__":
	main()
