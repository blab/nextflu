import csv
from collections import defaultdict
import json

countries = ['Guinea', 'Liberia', 'SierraLeone']
datestr = '2015-06-03'
case_numbers = {}

for country in countries:
    with open('data/'+datestr+'_'+country+'.csv') as table: 
        who_table = csv.DictReader(table)
        cases = defaultdict(int)
        for line in who_table:
            if line['Country'].replace(' ', '')==country and line['Epi week']!='':
                try:
                    n=int(line['Numeric'].split('.')[0]) 
                except:
                    print "not a number:",line['Numeric'] 
                    continue
                else:
                    tmp = line['Epi week'].split()[-1][1:-1].split('-')
                    date = int(tmp[0])+float(tmp[1][1:])*7/365.25
                    cases[date]+=n
                tmpcases = sorted(cases.items())
                case_numbers.update({'x'+country: [round(val[0],2) for val in tmpcases],
                                     country: [val[1] for val in tmpcases]})

with open('../auspice/data/case_numbers.json', 'w') as cn:
    json.dump(case_numbers, cn)
