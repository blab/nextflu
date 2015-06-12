import csv
from collections import defaultdict
import json
from datetime import datetime

fname = '../MERS-Cases/data/cases.csv'

case_numbers = defaultdict(lambda : defaultdict(int))

with open(fname) as table:
	mers_table = csv.DictReader(table)
	for line in mers_table:
		case_numbers[line['country']][line['reported']] +=1

countries = ['other', 'Korea/China','Saudi Arabia']
reduced_casenumbers = {c:defaultdict(int) for c in countries}


for c, cases in case_numbers.iteritems():
	if c in ['China', 'South Korea']:
		roughc='Korea/China'
	elif c in ['KSA']:
		roughc='Saudi Arabia'
	else:
		roughc='other'

	for rep_date, val in cases.iteritems():
		if rep_date!='':
			y,w,d = datetime.strptime(rep_date, '%Y-%m-%d').isocalendar()
			reduced_casenumbers[roughc][(y,w)]+=val

for roughc in countries:
	valid_dates = sorted(reduced_casenumbers[roughc].keys())
	print valid_dates
	for y in [2012, 2013, 2014,2015]:
		for w in range(1,53):
			if (y,w-1) in valid_dates or (y,w+1) in valid_dates:
				if (y,w-1) != valid_dates[-1]:
					reduced_casenumbers[roughc][(y,w)]+=0


processed_casenumbers = {}
for c, cases in reduced_casenumbers.iteritems():
	tmpcases = sorted(cases.items())
	processed_casenumbers.update({'x'+c:[round(x[0][0]+x[0][1]*7/365.25,2) for x in tmpcases],
								  c:[x[1] for x in tmpcases]})


with open('../auspice/data/mers_case_numbers.json', 'w') as cn:
    json.dump(processed_casenumbers, cn)
