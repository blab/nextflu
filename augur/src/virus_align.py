# take filtered sequences and align using muscle

import os, time
from io_util import *

def update_viruses(alignment, viruses):
	strain_to_sequence_map = {x['strain'].lower(): x['seq'] for x in alignment}
	for v in viruses:
		v['seq'] = strain_to_sequence_map[v['strain'].lower()]

def cleanup():
	for tmp_file in ['temp_in.fasta', 'temp_out.fasta']:
		try:
			os.remove(tmp_file)
		except OSError:
			pass

def main(viruses):

	print "--- Align at " + time.strftime("%H:%M:%S") + " ---"
	write_fasta(viruses, 'temp_in.fasta')
	os.system("mafft --nofft temp_in.fasta > temp_out.fasta")
	alignment = read_fasta('temp_out.fasta')
	update_viruses(alignment, viruses)
	cleanup()
	return viruses

if __name__ == "__main__":
	main()