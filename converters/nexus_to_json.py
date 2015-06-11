import dendropy
import json
import time
import argparse

from datetime import datetime as dt
from dateutil.parser import parse

def toYearFraction(date):
    def sinceEpoch(date): # returns seconds since epoch
        return time.mktime(date.timetuple())
    s = sinceEpoch

    year = date.year
    startOfThisYear = dt(year=year, month=1, day=1)
    startOfNextYear = dt(year=year+1, month=1, day=1)

    yearElapsed = s(date) - s(startOfThisYear)
    yearDuration = s(startOfNextYear) - s(startOfThisYear)
    fraction = yearElapsed/yearDuration

    return date.year + fraction

def parse_taxon_fields(taxon, node):
    fields = taxon.__str__().split("|")
    
    node["strain"] = fields[2]
    node["lab"] = fields[3]
    node["country"] = fields[4]
    node["region"] = fields[5]
    
    date = parse(fields[8])
    node["date"] = date.strftime('%Y-%m-%d')
    node["datevalue"] = toYearFraction(date)

def process_node(node, x, y):    
    y0 = 0.0
    
    node1 = {}
    if node.is_leaf():
        parse_taxon_fields(node.taxon, node1)
        node1['name'] = node.taxon.__str__()
        y0 = y[0]
        y[0] += 1
    else:
        node1['clade'] = 0;
        node1['children'] = [];
        for child in node.child_nodes():
            x1 = x + child.edge_length
            child1, y1 = process_node(child, x1, y)
            y0 += y1;
            node1['children'].append(child1)
            
        y0 /= len(node.child_nodes())
        
    node1['xvalue'] = x
    node1['yvalue'] = y0
    
    return node1, y0

def main():

	#using argparse incase switches need to be added
	parser = argparse.ArgumentParser(description='NEXUS to JSON converter')
	parser.add_argument('args', nargs=argparse.REMAINDER)
	args = parser.parse_args()

	tree = dendropy.Tree.get(
		path=args.args[0],
		schema='nexus',
		rooting="default-rooted")
 
	y = [0]
	root, y1 = process_node(tree.seed_node, 0.0, y)

	with open(args.args[1], 'w') as outfile:
		json.dump(root, outfile, indent=4)


if __name__ == "__main__":
	main()
