# take filtered sequences and align using muscle
# be careful to keep inline with canonical ordering (no extra gaps)

import os, json, subprocess
from Bio import SeqIO

OUTGROUP = 'A/Beijing/32/1992'

def read_viruses():
	try:
		handle = open('v_filter.json', 'r')  
	except IOError:
		pass
	else:	
  		viruses = json.load(handle)
  		handle.close()
	return viruses

def write_fasta(viruses):
	try:
		handle = open('temp_in.fasta', 'w') 
	except IOError:
		pass
	else:				
  		for v in viruses:
  			handle.write(">" + v['strain'] + "\n")
			handle.write(v['nt'] + "\n")  			
		handle.close()
		
def read_alignment():
	alignment = []
	try:
		handle = open('temp_out.fasta', 'r')
	except IOError:
		print "temp_out.fasta not found"
	else:
		for record in SeqIO.parse(handle, "fasta"):
			v = {
				"strain": record.description,
				"nt": str(record.seq)				
			}
			alignment.append(v)
		handle.close()
	return alignment
	
def mask_from_outgroup(alignment):
	outgroup_seq = filter(lambda v: v['strain'] == OUTGROUP, alignment)[0]['nt']	
	mask = outgroup_seq.replace('A', '1').replace('T', '1').replace('G', '1').replace('C', '1').replace('-', '0')
	for pair in alignment:
		filtered = ""
		for (c, m) in zip(list(pair['nt']), list(mask)):
			if m == '1':
				filtered += str(c)
		pair['nt'] = filtered
	
def update_viruses(alignment, viruses):
	for v in viruses:
		seq = filter(lambda x: x['strain'] == v['strain'], alignment)[0]['nt']
		v['nt'] = seq
	
def write_viruses(viruses):
	try:
		handle = open('v_align.json', 'w') 
	except IOError:
		pass
	else:				
		json.dump(viruses, handle, indent=2)
		handle.close()	
	
def cleanup():
	try:
		os.remove('temp_in.fasta')
	except OSError:
		pass
	try:
		os.remove('temp_out.fasta')
	except OSError:
		pass		

def main():

	print "--- Virus align ---"

	viruses = read_viruses()
	write_fasta(viruses)
	subprocess.call("muscle -in temp_in.fasta -out temp_out.fasta -diags -maxiters 2", shell=True)
	alignment = read_alignment()
	mask_from_outgroup(alignment)
	update_viruses(alignment, viruses)
	write_viruses(viruses)
	cleanup()
  		
if __name__ == "__main__":
    main()