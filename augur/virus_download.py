# loads sequences from GISAID
# outputs to virus_ingest.json

import os, time, json, sys
from selenium import webdriver
from Bio import SeqIO
from io_util import *

def download_gisaid(start_year, end_year):

	VAR = "ncszee"
	GISAID_FASTA = 'gisaid_epiflu_sequence.fasta'	

	# start fresh
	try:
		os.remove(GISAID_FASTA)
	except OSError:
		pass

	# start firefox driver
	print "Start webdriver"
	profile = webdriver.FirefoxProfile()
	profile.set_preference("browser.download.folderList",2)
	profile.set_preference("browser.download.manager.showWhenStarting",False)
	profile.set_preference("browser.download.dir", os.getcwd())
	profile.set_preference("browser.helperApps.neverAsk.saveToDisk", "application/octet-stream")
	driver = webdriver.Firefox(firefox_profile=profile)

	# open GISAID
	print "Open GISAID"
	driver.get('http://platform.gisaid.org/epi3/')
	assert 'GISAID' in driver.title

	# login
	print "Login"
	username = driver.find_element_by_name('login')
	username.send_keys(os.environ['GISAID_USER'])
	password = driver.find_element_by_name('password')
	password.send_keys(os.environ['GISAID_PASS'])
	driver.execute_script("return doLogin();")

	# navigate to EpiFlu
	print "Navigate to EpiFlu"
	time.sleep(10)
	driver.execute_script("return sys.call('c_" + VAR + "_4l','Go',new Object({'page':'epi3'}));")

	# fill in form and submit
	print "Fill in form and submit"
	time.sleep(10)
	type = driver.find_element_by_id('ce_nbp4cn_7x_select')
	type.send_keys("A")
	time.sleep(10)
	ha = driver.find_element_by_id('ce_nbp4cn_7y_select')
	ha.send_keys("3")
	time.sleep(10)
	na = driver.find_element_by_id('ce_nbp4cn_7z_select')
	na.send_keys("2")
	time.sleep(10)
	species = driver.find_element_by_id('ce_nbp4cn_81_select')
	species.send_keys("H")
	time.sleep(10)
	start_date = driver.find_element_by_id('ce_nbp4cn_85_input')
	start_date.send_keys(str(start_year)+"-01-01")
	time.sleep(10)
	end_date = driver.find_element_by_id('ce_nbp4cn_86_input')
	end_date.send_keys(str(end_year)+"-12-31")
	time.sleep(10)
	button = driver.find_element_by_xpath("//button[@accesskey='g']")
	button.click()

	# select sequences
	print "Select sequences"
	time.sleep(10)
	checkbox = driver.find_element_by_xpath("//span[@class='yui-dt-label']/input[@type='checkbox']")
	checkbox.click()
	time.sleep(30)
	button = driver.find_element_by_xpath("//div[@id='ce_nbp4cn_9a']//button")
	button.click()

	# set download options
	print "Set download options"
	time.sleep(10)
	driver.switch_to_frame(driver.find_element_by_tag_name("iframe"))
	checkbox = driver.find_element_by_xpath("//input[@name='ce_nbp4cn_9f_name' and @value='dna']")
	checkbox.click()
	time.sleep(10)
	checkbox = driver.find_element_by_xpath("//input[@name='ce_nbp4cn_9h_name' and @value='HA']")
	checkbox.click()

	# download
	print "Download"
	time.sleep(10)
	button = driver.find_element_by_xpath("//div[@id='ce_nbp4cn_9t']//button")
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
	
def main():

	print "--- Ingest at " + time.strftime("%H:%M:%S") + " ---"

	download_gisaid(2008, 2020)
	os.rename('gisaid_epiflu_sequence.fasta', 'data/gisaid_epiflu_sequence.fasta')

if __name__ == "__main__":
	main()
