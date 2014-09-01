# create phylogeny from alignment
# raxml will complain about identical sequences, but still includes them all in the resulting tree
# writes out newick.tree file

import os, re, time, glob
import subprocess
import dendropy
from io_util import *
							
RAXML_LIMIT = 1.0 # in hours			
							
def cleanup():
	for file in glob.glob("RAxML_*"):
		try:
			os.remove(file)
		except OSError:
			pass    	
		
def delimit_newick(infile_name, outfile_name):
	with open(infile_name, 'r') as file:
		newick = file.read().replace('\n', '')	
		newick = re.sub(r'(A/[^\:^,]+)', r"'\1'", newick)
	with open(outfile_name, 'w') as file:
		file.write(newick)	
									
def main():

	print "--- Tree infer at " + time.strftime("%H:%M:%S") + " ---"
		
	cleanup()
	viruses = read_json('data/virus_clean.json')
	write_fasta(viruses, 'temp.fasta')

	print "Building initial tree with FastTree"
	os.system("fasttree -gtr -nt -gamma -nosupport -mlacc 2 -slownni temp.fasta > initial_tree.newick")
	delimit_newick("initial_tree.newick", "temp.newick")
	tree = dendropy.Tree.get_from_path("temp.newick", "newick")	
	tree.resolve_polytomies()
	tree.write_to_path("initial_tree.newick", "newick")

	print "RAxML tree optimization with time limit " + str(RAXML_LIMIT) + " hours"
	os.system("seqmagick convert temp.fasta temp.phyx")
	process = subprocess.Popen("raxml -f d -T 6 -j -s temp.phyx -n topology -c 25 -m GTRCAT -p 344312987 -t initial_tree.newick", shell=True)
	time.sleep(int(RAXML_LIMIT*3600))
	process.terminate()
	
	last_tree_file = [file for file in glob.glob("RAxML_checkpoint*")][-1]	
	os.rename(last_tree_file, 'final_tree.newick')
		
	print "RAxML branch length optimization"
	os.system("raxml -f e -T 6 -j -s temp.phyx -n branches -c 25 -m GTRGAMMA -p 344312987 -t final_tree.newick")
	os.rename('RAxML_result.branches', 'data/tree.newick')
	cleanup()	

if __name__ == "__main__":
    main()