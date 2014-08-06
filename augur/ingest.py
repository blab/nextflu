# loads sequences from GISAID
# outputs to v_ingest.json

import os, time, json
from selenium import webdriver
from Bio import SeqIO

GISAID_FASTA = 'gisaid_epiflu_sequence.fasta'

def download_gisaid(start_year, end_year):

	# start fresh
	try:
		os.remove(GISAID_FASTA)
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
	if not os.path.isfile(GISAID_FASTA):
		time.sleep(30)
	if not os.path.isfile(GISAID_FASTA):
		time.sleep(30)
	if not os.path.isfile(GISAID_FASTA):
		time.sleep(30)
	if not os.path.isfile(GISAID_FASTA):
		time.sleep(30)			
	if os.path.isfile(GISAID_FASTA):						
		while (os.stat(GISAID_FASTA).st_size == 0):
			time.sleep(5)
	
	# close driver
	driver.quit()
	
def parse_gisaid():	
	viruses = []
	try:
		handle = open(GISAID_FASTA, 'r')
	except IOError:
		print GISAID_FASTA + " not found"
	else:
		for record in SeqIO.parse(handle, "fasta"):
			words = record.description.replace(">","").replace(" ","").split('|')
			strain = words[0]
			accession = words[1]
			passage = words[3]
			date = words[5]
			nt = str(record.seq).upper()
			v = {
				"strain": strain,	
				"date": date,	
				"accession": accession,
				"db": "gisaid",
				"nt": nt				
			}
			if passage != "":
				v['passage'] = passage
			viruses.append(v)
		handle.close()	
		
	# start fresh
	try:
		os.remove(GISAID_FASTA)
	except OSError:
		pass # okay		
		
	return viruses
	
def download_and_parse_gisaid(start_year, end_year):

	download_gisaid(start_year, end_year)	# leaves GISAID_FASTA in dir
	return parse_gisaid()
		
def main():

	print "--- Virus ingest ---"

	prior_length = 0
	try:
		handle = open('v_ingest.json', 'r')  
	except IOError:
		pass
	else:	
  		viruses = json.load(handle)
  		prior_length = len(viruses)
  		handle.close()

	viruses = []
	
	print "Downloading 1990 to 2008 viruses"	
	viruses.extend(download_and_parse_gisaid(1990, 2008))
	
	print "Downloading 2009 to 2013 viruses"	
	viruses.extend(download_and_parse_gisaid(2009, 2013))
	
	print "Downloading 2014+ viruses"	
	viruses.extend(download_and_parse_gisaid(2014, 2020))		
	current_length = len(viruses)
	
	if (current_length > prior_length):
		print "Writing new v_ingest.json with " + str(current_length) + " viruses"
		try:
			handle = open('v_ingest.json', 'w') 
		except IOError:
			pass
		else:				
  			json.dump(viruses, handle, indent=2)
  			handle.close()
  	else:
  		print "Keeping old v_ingest.json with " + str(prior_length) + " viruses"
  		
if __name__ == "__main__":
    main()
