import os, time
from selenium import webdriver

driver = webdriver.Firefox()
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
button = driver.find_element_by_xpath("//button[@accesskey='g']")
button.click()

# download
time.sleep(10)
checkbox = driver.find_element_by_xpath("//span[@class='yui-dt-label']/input[@type='checkbox']")
checkbox.click()
time.sleep(30)
button = driver.find_element_by_xpath("//div[@id='ce_n7h14c_9f']/button")
button.click() # this is failing because of the 20k isolate limit
