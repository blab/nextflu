# shared variables and functions

import os, json
from Bio import SeqIO

OUTGROUP = 'A/Beijing/32/1992'

def read_json(file_name):
	try:
		handle = open(file_name, 'r')  
	except IOError:
		pass
	else:	
  		data = json.load(handle)
  		handle.close()
	return data
	
def write_json(data, file_name):
	try:
		handle = open(file_name, 'w') 
	except IOError:
		pass
	else:				
		json.dump(data, handle, indent=2)
		handle.close()
		
# fasta only has
# strain, seq

def read_fasta(file_name):
	alignment = []
	try:
		handle = open(file_name, 'r')
	except IOError:
		pass
	else:
		for record in SeqIO.parse(handle, "fasta"):
			v = {
				"strain": record.description,
				"seq": str(record.seq)				
			}
			alignment.append(v)
		handle.close()
	return alignment

def write_fasta(viruses, file_name):
	try:
		handle = open(file_name, 'w') 
	except IOError:
		pass
	else:				
  		for v in viruses:
  			handle.write(">" + v['strain'] + "\n")
			handle.write(v['seq'] + "\n")  			
		handle.close()
