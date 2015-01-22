import os, schedule, time, sys

def config():
	try:
		handle = open(os.path.expanduser("~") + "/.s3cfg", 'w')
	except IOError:
		pass
	else:
		handle.write("[default]\n")
		handle.write("access_key = " + os.environ['S3_KEY'] + "\n")
		handle.write("secret_key = " + os.environ['S3_SECRET'] + "\n")		
		handle.close()	
	os.system("s3cmd mb s3://" + os.environ['S3_BUCKET'])			

def data():
	print "Syncing data with S3"
	os.system("s3cmd sync --acl-public auspice/ s3://" + os.environ['S3_BUCKET'] + "/auspice/")	

def str2bool(obj):
	if isinstance(obj, basestring):
		return obj.lower() in ("yes", "true", "t", "1")
	if isinstance(obj, bool):
		return obj
	else:
		return False

def main(argv):
	"""Sync data with Amazon S3"""
	
	clock = False
	if len(argv) > 0:
		clock = str2bool(argv[0])
	if clock:
		print "Setting clock with ntpdate"
		os.system("ntpdate ntp.ubuntu.com")
	
	config()
	data()
	schedule.every().minute.do(data)	
	while True:
		schedule.run_pending()
		time.sleep(1)		
	
if __name__ == "__main__":
    main(sys.argv[1:])
