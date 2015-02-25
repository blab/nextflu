# clean, reroot, ladderize newick tree
# output to tree.json
import numpy as np
from io_util import *
from tree_util import *
from tree_ancestral import *
from seq_util import *
import time
from collections import defaultdict
from matplotlib import pyplot as plt

min_titer = 10.0

def strain_name_fixing(name):
	return name.replace(' ','').lower()

def parse_HI_matrix(fname):
	import csv
	def titer_to_number(val):
		try:
			if '<' in val:
				return min_titer
			else:
				return float(val)
		except:
			print "Bad HI measurement:", val
			return np.nan

	with open(fname) as infile:
		csv_reader = csv.reader(infile)

		# parse sera
		row1 = csv_reader.next()
		row2 = csv_reader.next()
		print row1
		print row2
		ref_sera = [(e1+'/'+e2).replace(' ','') for e1,e2 in zip(row1, row2)[4:]]
		print ref_sera
		for row in csv_reader: # advance until the reference virus
			if row[0].startswith('REFERENCE'):
				break

		ref_strains = []
		ref_matrix = []
		for row in csv_reader: 
			if row[0].startswith('TEST'):
				break
			else: # load matrices until the test virus section starts
				ref_strains.append(strain_name_fixing(row[0]))
				ref_matrix.append(map(titer_to_number, row[4:]))

		test_strains = []
		test_matrix = []
		for row in csv_reader: # load test viruses until it is no longer an A/ flu  name
			if not row[0].startswith('A/'):
				break
			else:
				test_strains.append(strain_name_fixing(row[0]))
				test_matrix.append(map(titer_to_number, row[4:]))
		return ref_strains, np.array(ref_matrix), test_strains, np.matrix(test_matrix)


def normalize_HI_matrices(ref_matrix, test_matrix):
	self_avidity = np.diag(ref_matrix)
	norm_test_matrix = np.array(test_matrix[:,:len(self_avidity)])
	norm_ref_matrix = np.array(ref_matrix[:,:len(self_avidity)])
	for ai, av in enumerate(self_avidity):
		print av
		norm_test_matrix[:,ai] -= av
		norm_ref_matrix[:,ai]  -= av

	return -norm_ref_matrix, -norm_test_matrix

def read_tables():
	import glob
	from itertools import product
	flist = glob.glob('../../../flu_HI_data/tables/Flu*csv')
	all_names = set()
	all_measurements = defaultdict(list)
	HI_matrices = []
	for fname in flist:
		ref_names, ref_matrix, test_names, test_matrix = parse_HI_matrix(fname)
#		try:
#			ref_matrix, test_matrix = normalize_HI_matrices(ref_matrix, test_matrix)
#		except:
#			continue
		HI_matrices.append( (ref_names, ref_matrix[:,:len(ref_names)], test_names, test_matrix[:,:len(ref_names)] ))
		all_names.update(ref_names)	
		all_names.update(test_names)

		for test, ref in product(test_names, ref_names):
			all_measurements[(test, ref)].append(test_matrix[test_names.index(test), ref_names.index(ref)])
		for ref_test, ref in product(ref_names, ref_names):
			all_measurements[(ref_test, ref)].append(ref_matrix[ref_names.index(ref_test), ref_names.index(ref)])

		print fname, len(all_measurements), "measurements"

	print "total from tables:", len(all_measurements), "measurements"
	trevor_table = '../../../flu_HI_data/tables/trevor_elife_H3N2_HI_data.tsv'
	import csv
	with open(trevor_table) as infile:
		table_reader = csv.reader(infile, delimiter="\t")
		header = table_reader.next()
		for row in table_reader:
			try:
				all_names.add(strain_name_fixing(row[1]))
				all_measurements[(strain_name_fixing(row[1]), strain_name_fixing(row[4]))].append(float(row[6]))
			except:
				pass #print row

	print "grand total:", len(all_measurements), "measurements"
	return all_names, all_measurements, HI_matrices

def get_normalized_HI_titers():
	names, measurements, tables = read_tables()
	normalized_measurements = {}
	sera = set()
	HI_strains = set()
	for (test, ref), val in measurements.iteritems():
		HI_strains.add(test.lower())
		HI_strains.add(ref.lower())
		if (ref,ref) in measurements:
			sera.add(ref)
			normalized_val = np.log2(measurements[(ref, ref)]).mean() - np.log2(val).mean()
			if normalized_val<=10:
				normalized_measurements[(test, ref)] = max(0,normalized_val)
	return normalized_measurements, HI_strains, sera

def get_path_dendropy(tree, v1, v2):
	p1 = [v1]
	p2 = [v2]
	while p1[-1].parent_node != tree.seed_node:
		p1.append(p1[-1].parent_node)
	p1.append(tree.seed_node)
	p1.reverse()

	while p2[-1].parent_node != tree.seed_node:
		p2.append(p2[-1].parent_node)
	p2.append(tree.seed_node)
	p2.reverse()

#	import pdb; pdb.set_trace()
	for pi, (tmp_v1, tmp_v2) in enumerate(izip(p1,p2)):
		if tmp_v1!=tmp_v2:
			break
	path = p1[pi:] + p2[pi:]
	return path

def get_path_biopython(tree, v1, v2):
	#print "getting path for", v1.name, v2.name
	p1 = tree.get_path(v1)
	p2 = tree.get_path(v2)
	path = None
	for ii, (c1,c2) in enumerate(zip(p1, p2)):
		if c1!=c2:
			#print "found match", ii
			path = p1[ii:] + p2[ii:]
			break

	return path

def calc_tree_HI_distance(v1,v2, tree, lca):
	dist = 0
	tmp_lca = lca[(v1,v2)]
	for tmp_v in [v1, v2]:
		while tmp_v.parent_node != tmp_lca:
			dist+=tmp_v.dHI
			tmp_v1 = tmp_v.parent_node

	return dist

def fit_func(dHI, tree_graph, distances):
	return np.sum( (distances - np.dot(tree_graph, dHI))**2 )


def prune_tree(tree, taxa_with_HI):
	seq_strains = set([leaf.label for leaf in tree.leaf_iter()])
	print seq_strains.intersection(taxa_with_HI)
	tree.retain_taxa_with_labels(taxa_with_HI)


def main(tree_fname = 'data/tree_ancestral.json', HI_fname='data/HI_titers.txt'):

	print "--- Fitting HI titers at " + time.strftime("%H:%M:%S") + " ---"

	tree =  json_to_dendropy(read_json(tree_fname))
	HI_distances = load_HI_distances(HI_fname)

def get_strains_with_HI_and_sequence():
	names, measurements, HI_matrices = read_tables()
	from Bio import SeqIO
	good_strains = set()
	with open("data/strains_with_HI.fasta", 'w') as outfile, \
		 open("data/20150222_all_H3N2_HA1.fasta", 'r') as infile:
		for seq_rec in SeqIO.parse(infile, 'fasta'):
			reduced_name = strain_name_fixing(seq_rec.name)
			if reduced_name in names and (reduced_name not in good_strains):
				SeqIO.write(seq_rec, outfile,'fasta')
				good_strains.add(reduced_name)
				print seq_rec.name

def mark_HI_strains(tree, HI_strains, sera):
	for leaf in tree.leaf_iter():
		if leaf.strain.lower() in HI_strains:
			leaf.serum = leaf.strain.lower() in sera
			leaf.HI_info= True
		else:
			leaf.serum, leaf.HI_info=False, False

	for node in tree.postorder_node_iter():
		if not node.is_leaf():
			node.HI_info = any([c.HI_info for c in node.child_nodes()])

	# combine sets of branches that span identical sets of HI measurements
	branch_count = 0
	HI_split_count = 0  # HI measurment split counter
	HI_split_to_branch = defaultdict(list)
	for node in tree.preorder_node_iter():
		node.dHI, node.cHI, node.constrains, node.branch_index =0, 0, 0, branch_count
		if node.is_internal(): node.serum=False
		branch_count+=1
		if node.HI_info:
			node.HI_branch_index = HI_split_count
			HI_split_to_branch[node.HI_branch_index].append(node.branch_index)
			if node.is_leaf() or sum([c.HI_info for c in node.child_nodes()])>1:
				HI_split_count+=1

	return HI_split_to_branch, HI_split_count


def add_mutations(tree):
	'''
	add amino acid mutations to the tree
	'''
	tree.seed_node.mutations= ''
	for node in tree.postorder_node_iter():
		if node is not tree.seed_node:
			node.mutations = [a+str(pos-15)+b for pos, (a,b) in 
							enumerate(izip(node.parent_node.aa_seq, node.aa_seq)) if a!=b]


def map_HI_to_tree(tree, measurements, method = 'nnls', lam=10):
	names_to_clades = {leaf.strain.lower(): leaf for leaf in tree.leaf_iter()}
	sera = list(set([x[1].lower() for x in measurements]))
	HI_names = list(set([x[0].lower() for x in measurements] + sera))
	# assign indices to branches
	HI_split_map, HI_sc = mark_HI_strains(tree, HI_names, sera)
	add_mutations(tree)
	print "# of reference strains:",len(sera), "# of branches with HI contraint", HI_sc

	tree_graph = []
	HI_diff = []
	for (test, ref), val in measurements.iteritems():
		if not np.isnan(val):
			try:
				if test != ref and ref in names_to_clades and test in names_to_clades:
					path = get_path_dendropy(tree, names_to_clades[test], names_to_clades[ref])
					tmp = np.zeros(HI_sc + len(sera))
					branches = np.unique([c.HI_branch_index for c in path])
					tmp[branches] = 1
					tmp[HI_sc+sera.index(ref)] = 1
					tree_graph.append(tmp)
					HI_diff.append(val)
#				else:
#					print test, ref, "not found"
			except:
				print test, ref, "ERROR"

	HI_diff =  np.array(HI_diff)
	tree_graph= np.array(tree_graph)
	print "matrix dimensions: ", tree_graph.shape, "measurements x parameters"

	if method=='l1reg':  # l1 regularized fit, no constraint on sign of effect
		from l1regls import l1regls
		A = matrix(tree_graph)
		b = matrix(HI_diff)
		w = np.array([x for x in l1regls(A/np.sqrt(lam),b/np.sqrt(lam))])
		print 'l1reg', fit_func(w, tree_graph, HI_diff), np.sum((HI_diff-1)**2)
	elif method=='nnls':  # non-negative least square, not regularized
		from scipy.optimize import nnls
		w = nnls(tree_graph, HI_diff)[0]
		print 'nnls', fit_func(w, tree_graph, HI_diff), np.sum((HI_diff-1)**2)
	elif method=='nnl2reg':	# non-negative L2 norm regularized fit
		from cvxopt import matrix, solvers
		P = matrix(np.dot(tree_graph.T, tree_graph) + lam*np.eye(tree_graph.shape[1]))
		q = matrix( -np.dot( HI_diff, tree_graph))
		h = matrix(np.zeros(tree_graph.shape[1])) # Gw <=h
		G = matrix(-np.eye(tree_graph.shape[1]))
		W = solvers.qp(P,q,G,h)
		w = np.array([x for x in W['x']])
		print 'QP', fit_func(w, tree_graph, HI_diff), np.sum((HI_diff-1)**2)
	elif method=='nnl1reg':  # non-negative fit, branch terms L1 regularized, avidity terms L2 regularized
		from cvxopt import matrix, solvers
		P1 = np.zeros((tree_graph.shape[1]+HI_sc,tree_graph.shape[1]+HI_sc))
		P1[:tree_graph.shape[1], :tree_graph.shape[1]] = np.dot(tree_graph.T, tree_graph)
		for ii in xrange(HI_sc, tree_graph.shape[1]):
			P1[ii,ii]+=lam
		P = matrix(P1)

		q1 = np.zeros(tree_graph.shape[1]+HI_sc)
		q1[:tree_graph.shape[1]] = -np.dot( HI_diff, tree_graph)
		q1[tree_graph.shape[1]:] = lam
		q = matrix(q1)

		h = matrix(np.zeros(2*HI_sc)) 	# Gw <=h
		G1 = np.zeros((2*HI_sc,tree_graph.shape[1]+HI_sc))
		G1[:HI_sc, :HI_sc] = -np.eye(HI_sc)
		G1[HI_sc:, :HI_sc] = np.eye(HI_sc)
		G1[HI_sc:, tree_graph.shape[1]:] = -np.eye(HI_sc)
		G = matrix(G1)
		W = solvers.qp(P,q,G,h)
		w = np.array([x for x in W['x']])[:tree_graph.shape[1]]
		print 'QP l1', fit_func(w, tree_graph, HI_diff), np.sum((HI_diff-1)**2)

	# for each set of branches with HI constraints, pick the branch with most aa mutations
	# and assign the dHI to that one, record the number of constraints
	index_to_branch = {n.branch_index:n for n in tree.postorder_node_iter()}
	for HI_split, branches in HI_split_map.iteritems():
		likely_branch = branches[np.argmax([len(index_to_branch[b].mutations) for b in branches])]
		index_to_branch[likely_branch].dHI = w[HI_split]
		index_to_branch[likely_branch].constraints = tree_graph[:,HI_split].sum()

	# integrate the HI change dHI into a cumulative antigentic evolution score cHI
	for node in tree.preorder_node_iter():
		if node!=tree.seed_node:
			node.cHI = node.parent_node.cHI + node.dHI
	return tree, {serum:w[HI_sc+ii] for ii, serum in enumerate(sera)}

def predict_HI(virus, serum, tree,avidities):
	v = tree.find_node_with_taxon_label(virus)
	s = tree.find_node_with_taxon_label(serum)
	if v is not None and s is not None:
		path = get_path_dendropy(tree, v,s)
		return avidities[serum] + np.sum(b.dHI for b in path)
	else:
		return None

def plot_tree(tree):
	btree = to_Biopython(tree)
	color_BioTree_by_attribute(btree,"cHI", transform = lambda x:x)
	Phylo.draw(btree, label_func = lambda  x: 'X' if x.serum else 'o' if x.HI_info else '', 
		show_confidence= False) #, branch_labels = lambda x:x.mutations)

def plot_dHI_distribution(tree):
	plt.figure()
	mut_dHI = sorted([(c.dHI, c.mutations) for c in tree.postorder_node_iter() if c.HI_info], 
					reverse=True, key=lambda x:x[0])
	thres = [0,1,2,3,100]
	for lower, upper in zip(thres[:-1], thres[1:]):
		tmp = [x[0] for x in mut_dHI if len(x[1])>=lower and len(x[1])<upper]
		plt.plot(sorted(tmp), np.linspace(1,0, len(tmp)), label='#aa='+str(lower)+', total '+str(len(tmp)))
	plt.legend(loc=1)
	plt.show()

if __name__ == "__main__":
	from Bio import Phylo
	reg = 5
	#tree_fname = 'data/tree_long_HI.json'
	tree_fname = 'data/tree_refine.json'
	tree =  json_to_dendropy(read_json(tree_fname))
	normalized_measurements, HI_strains, sera = get_normalized_HI_titers()
	test_HI, train_HI = {}, {}
	for key, val in normalized_measurements.iteritems():
		if np.random.uniform()<0.2:
			test_HI[key]=val
		else:
			train_HI[key]=val

	tree,avidities = map_HI_to_tree(tree, train_HI, method='nnl1reg', lam = reg)
	validation = {}
	for key, val in test_HI.iteritems():
		pred_HI = predict_HI(key[0], key[1], tree, avidities)
		if pred_HI is not None:
			validation[key] = (val, pred_HI)
	a = np.array(validation.values())
	plt.scatter(a[:,0], a[:,1])
	plt.xlabel("measured")
	plt.ylabel("predicted")
	plt.xlabel("measured")
	plt.title('reg='+str(reg)+', avg error '+str(round(np.mean(np.abs(a[:,0]-a[:,1])),3)))

	out_tree_fname = 'data/tree_HI.json'
	write_json(dendropy_to_json(tree.seed_node), out_tree_fname, indent=None)

	#plot_tree(tree)
	#plot_dHI_distribution(tree)

