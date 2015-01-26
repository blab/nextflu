import os, time, sys
import subprocess
import virus_download
from io_util import *

def str2bool(obj):
	if isinstance(obj, basestring):
		return obj.lower() in ("yes", "true", "t", "1")
	if isinstance(obj, bool):
		return obj
	else:
		return False

def main(argv):
	"""Download sequence data"""
	"""Pass True/False for headless"""

	print "--- Start ingest at " + time.strftime("%H:%M:%S") + " ---"

	processes = []

	headless = False
	if len(argv) > 0:
		headless = str2bool(argv[0])
	if headless:
		processes.append(subprocess.Popen("exec Xvfb :99 -shmem -screen 0 1366x768x16", shell=True))
		processes.append(subprocess.Popen("exec x11vnc -display :99 -N -forever", shell=True))
		time.sleep(10)

	virus_download.main()

	for process in processes:
		process.terminate()

if __name__ == "__main__":
	main(sys.argv[1:])
