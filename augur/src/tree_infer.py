# create phylogeny from alignment
# raxml will complain about identical sequences, but still includes them all in the resulting tree
# writes out newick.tree file

import os, re, time, glob, shutil
import subprocess
import dendropy
from io_util import *
from tree_util import delimit_newick

def cleanup():
	for file in glob.glob("RAxML_*") + glob.glob("temp*") + ["raxml_tree.newick", "initial_tree.newick"]:
		try:
			os.remove(file)
		except OSError:
			pass

def main(viruses, raxml_time_limit, outgroup):

	print "--- Tree infer at " + time.strftime("%H:%M:%S") + " ---"

	cleanup()
	write_fasta(viruses, 'temp.fasta')
	print "Building initial tree with FastTree"
	os.system("fasttree -gtr -nt -gamma -nosupport -mlacc 2 -slownni temp.fasta > initial_tree.newick")
	delimit_newick("initial_tree.newick", "temp.newick")
	tree = dendropy.Tree.get_from_path("temp.newick", "newick")
	tree.resolve_polytomies()
	tree.write_to_path("initial_tree.newick", "newick")

	os.system("seqmagick convert temp.fasta temp.phyx")
	if raxml_time_limit>0:
		print "RAxML tree optimization with time limit " + str(raxml_time_limit) + " hours"
		# using exec to be able to kill process
		end_time = time.time() + int(raxml_time_limit*3600)
		process = subprocess.Popen("exec raxml -f d -T 6 -j -s temp.phyx -n topology -c 25 -m GTRCAT -p 344312987 -t initial_tree.newick", shell=True)
		while (time.time() < end_time):
			if os.path.isfile('RAxML_result.topology'):
				break
			time.sleep(10)
		process.terminate()

		checkpoint_files = [file for file in glob.glob("RAxML_checkpoint*")]
		if os.path.isfile('RAxML_result.topology'):
			checkpoint_files.append('RAxML_result.topology')
		if len(checkpoint_files) > 0:
			last_tree_file = checkpoint_files[-1]
			shutil.copy(last_tree_file, 'raxml_tree.newick')
		else:
			shutil.copy("initial_tree.newick", 'raxml_tree.newick')
	else:
		shutil.copy("initial_tree.newick", 'raxml_tree.newick')

	print "RAxML branch length optimization and rooting"
	os.system("raxml -f e -T 6 -s temp.phyx -n branches -c 25 -m GTRGAMMA -p 344312987 -t raxml_tree.newick -o " + outgroup)

	out_fname = "data/tree_infer.newick"
	os.rename('RAxML_result.branches', out_fname)
	cleanup()	
	return out_fname;

if __name__ == "__main__":
	main()
