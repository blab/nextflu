import os

def main():
	"""Upload to Amazon S3"""
	
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
	os.system("s3cmd sync --acl-public data/ s3://" + os.environ['S3_BUCKET'])	
	
if __name__ == "__main__":
    main()