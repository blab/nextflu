import os

def main():
	"""Upload to Amazon S3"""
	os.system("s3_website cfg apply --headless")	
	os.system("s3_website push")		

if __name__ == "__main__":
    main()