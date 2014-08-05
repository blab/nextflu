import os, time
from selenium import webdriver
from Bio import SeqIO
import rethinkdb as r
from rethinkdb.errors import RqlRuntimeError, RqlDriverError

RDB_HOST =  os.environ.get('RDB_HOST') or 'localhost'
RDB_PORT = os.environ.get('RDB_PORT') or 28015
RDB_DB = 'augur'

def download_gisaid(start_year, end_year):

	# start fresh
	try:
		os.remove('gisaid_epiflu_sequence.fasta')
	except OSError:
		pass # okay

	# start firefox driver
	profile = webdriver.FirefoxProfile()
	profile.set_preference("browser.download.folderList",2)
	profile.set_preference("browser.download.manager.showWhenStarting",False)
	profile.set_preference("browser.download.dir", os.getcwd())
	profile.set_preference("browser.helperApps.neverAsk.saveToDisk", "application/octet-stream")
	driver = webdriver.Firefox(firefox_profile=profile)

	# open GISAID 
	driver.get('http://platform.gisaid.org/epi3/')
	assert 'GISAID' in driver.title

	# login
	username = driver.find_element_by_name('login')
	username.send_keys(os.environ['GISAID_USER'])
	password = driver.find_element_by_name('password')
	password.send_keys(os.environ['GISAID_PASS'])
	driver.execute_script("return doLogin();")

	# navigate to EpiFlu
	time.sleep(10)
	driver.execute_script("return sys.call('c_n7h14c_4l','Go',new Object({'page':'epi3'}));")

	# fill in form and submit
	time.sleep(10)
	type = driver.find_element_by_id('ce_n7h14c_82_select')
	type.send_keys("A")
	time.sleep(10)
	ha = driver.find_element_by_id('ce_n7h14c_83_select')
	ha.send_keys("3")
	time.sleep(10)
	na = driver.find_element_by_id('ce_n7h14c_84_select')
	na.send_keys("2")
	time.sleep(10)
	species = driver.find_element_by_id('ce_n7h14c_86_select')
	species.send_keys("H")
	time.sleep(10)
	start_date = driver.find_element_by_id('ce_n7h14c_8a_input')
	start_date.send_keys(str(start_year)+"-01-01")
	time.sleep(10)
	end_date = driver.find_element_by_id('ce_n7h14c_8b_input')
	end_date.send_keys(str(end_year)+"-12-31")
	time.sleep(10)	
	button = driver.find_element_by_xpath("//button[@accesskey='g']")
	button.click()

	# select sequences
	time.sleep(10)
	checkbox = driver.find_element_by_xpath("//span[@class='yui-dt-label']/input[@type='checkbox']")
	checkbox.click()
	time.sleep(30)
	button = driver.find_element_by_xpath("//div[@id='ce_n7h14c_9f']//button")
	button.click()

	# set download options
	time.sleep(10)
	driver.switch_to_frame(driver.find_element_by_tag_name("iframe"))
	checkbox = driver.find_element_by_xpath("//input[@name='ce_n7h14c_9k_name' and @value='dna']")
	checkbox.click()
	time.sleep(10)
	checkbox = driver.find_element_by_xpath("//input[@name='ce_n7h14c_9m_name' and @value='HA']")
	checkbox.click()

	# download
	time.sleep(10)
	button = driver.find_element_by_xpath("//div[@id='ce_n7h14c_9y']//button")
	button.click()
	
	# wait for download to complete
	time.sleep(60)		
	if os.path.isfile('gisaid_epiflu_sequence.fasta'):
		file_size = os.stat('gisaid_epiflu_sequence.fasta').st_size
		time.sleep(10)	
		while (file_size < os.stat('gisaid_epiflu_sequence.fasta').st_size):
			file_size = os.stat('gisaid_epiflu_sequence.fasta').st_size
			time.sleep(10)	

def add_gisaid(start_year, end_year):
	
	# Each thread or multiprocessing PID should be given its own connection object.	
	try:
		connection = r.connect(host=RDB_HOST, port=RDB_PORT, db=RDB_DB)
	except RqlDriverError:
		abort(503, "No database connection could be established.")
	
	try:
		handle = open('gisaid_epiflu_sequence.fasta', "rU")
	except IOError:
		print "gisaid_epiflu_sequence.fasta not found"

	for record in SeqIO.parse(handle, "fasta"):

		# parse fasta
		words = record.description.replace(">","").replace(" ","").split('|')
		strain = words[0]
		accession = words[1]
		date = words[5]
		seq = str(record.seq).upper()
		entry = {
			"id": accession,
			"db": "gisaid",
			"strain": strain,		
			"date": date,
			"seq": seq				
		}

		# push to db
		r.table('seq').insert(entry).run(connection)
			
	handle.close()	
	
	try:
		connection.close()
	except RqlDriverError:
		abort(503, "Database connection broken.")					
		
def main():

#	download_gisaid(2010, 2020)
	add_gisaid(2010, 2020)
		
	
if __name__ == "__main__":
    main()
