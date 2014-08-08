# create phylogeny from alignment
# raxml will complain about identical sequences, but still includes them all in the resulting tree
# writes out newick.tree file

import os, re, json, time, seqmagick
from Bio import Phylo
from share import *

RAXML = 'raxmlHPC-PTHREADS-AVX'

def cleanup():
	try:
		os.remove('temp.fasta')
	except OSError:
		pass
	try:
		os.remove('temp.phyx')
	except OSError:
		pass
	try:
		os.remove('temp.phyx.reduced')
	except OSError:
		pass		
	try:
		os.remove('RAxML_info.out')
	except OSError:
		pass	
	try:
		os.remove('RAxML_log.out')
	except OSError:
		pass		
	try:
		os.remove('RAxML_parsimonyTree.out')
	except OSError:
		pass
	try:
		os.remove('RAxML_result.out')
	except OSError:
		pass
										

def main():

	print "--- Tree at " + time.strftime("%H:%M:%S") + " ---"
	
	viruses = read_viruses('virus_clean.json')
	write_fasta(viruses, 'temp.fasta')
	os.system("seqmagick convert temp.fasta temp.phyx")
	os.system(RAXML + " -T 6 -s temp.phyx -n out -c 25 -f d -m GTRCAT -p 344312987")
	os.rename('RAxML_bestTree.out', 'newick.tree')
	cleanup()

	tree = Phylo.read('newick.tree', 'newick')
	tree.ladderize(reverse=True) 
#	Phylo.draw_ascii(tree)
	
if __name__ == "__main__":
    main()