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


class HI_tree(object):

	def __init__(self, tree, HI_measurments):
		self.tree = tree
		self.HI = HI_measurments
		self.normalize_HI()
		self.add_mutations()
		self.mark_HI_splits()
		self.names_to_clades = {leaf.strain.lower(): leaf for leaf in self.tree.leaf_iter()}

	def normalize_HI(self):
		self.HI_normalized = {}
		sera = set()
		HI_strains = set()
		for (test, ref), val in self.HI.iteritems():
			HI_strains.add(test.lower())
			HI_strains.add(ref.lower())
			if (ref,ref) in self.HI:
				sera.add(ref)
				normalized_val = np.log2(self.HI[(ref, ref)]).mean() - np.log2(val).mean()
				self.HI_normalized[(test, ref)] = max(0,normalized_val)
		self.sera = list(sera)
		self.HI_strains = list(HI_strains)

	def add_mutations(self):
		'''
		add amino acid mutations to the tree
		'''
		self.tree.seed_node.mutations= ''
		for node in self.tree.postorder_node_iter():
			if node is not tree.seed_node:
				node.mutations = [a+str(pos-15)+b for pos, (a,b) in 
								enumerate(izip(node.parent_node.aa_seq, node.aa_seq)) if a!=b]

	def mark_HI_splits(self):
		for leaf in self.tree.leaf_iter():
			if leaf.strain.lower() in self.HI_strains:
				leaf.serum = leaf.strain.lower() in self.sera
				leaf.HI_info= True
			else:
				leaf.serum, leaf.HI_info=False, False

		for node in self.tree.postorder_node_iter():
			if not node.is_leaf():
				node.HI_info = any([c.HI_info for c in node.child_nodes()])

		# combine sets of branches that span identical sets of HI measurements
		self.HI_split_count = 0  # HI measurment split counter
		self.HI_split_to_branch = defaultdict(list)
		for node in self.tree.preorder_node_iter():
			node.dHI, node.cHI, node.constraints =0, 0, 0
			if node.is_internal(): node.serum=False
			if node.HI_info:
				node.HI_branch_index = self.HI_split_count
				self.HI_split_to_branch[node.HI_branch_index].append(node)
				if node.is_leaf() or sum([c.HI_info for c in node.child_nodes()])>1:
					self.HI_split_count+=1

		print "# of reference strains:",len(self.sera), "# of branches with HI contraint", self.HI_split_count

	def get_path(self, v1, v2):
		if v1 in self.names_to_clades and v2 in self.names_to_clades:
			p1 = [self.names_to_clades[v1]]
			p2 = [self.names_to_clades[v2]]
			for tmp_p in [p1,p2]:
				while tmp_p[-1].parent_node != self.tree.seed_node:
					tmp_p.append(tmp_p[-1].parent_node)
				tmp_p.append(tree.seed_node)
				tmp_p.reverse()

			for pi, (tmp_v1, tmp_v2) in enumerate(izip(p1,p2)):
				if tmp_v1!=tmp_v2:
					break
			path = p1[pi:] + p2[pi:]
		else:
			path = None
		return path

	def make_treegraph(self):
		tree_graph = []
		HI_dist = []
		for (test, ref), val in self.train_HI.iteritems():
			if not np.isnan(val):
				if True: #try:
					if test != ref  and ref in self.names_to_clades \
									and test in self.names_to_clades:
						path = self.get_path(test, ref)
						tmp = np.zeros(self.HI_split_count + len(self.sera))
						branches = np.unique([c.HI_branch_index for c in path])
						tmp[branches] = 1
						tmp[self.HI_split_count+self.sera.index(ref)] = 1
						tree_graph.append(tmp)
						HI_dist.append(val)
#				except:
#					print test, ref, "ERROR"

		self.HI_dist =  np.array(HI_dist)
		self.tree_graph= np.array(tree_graph)
		print "Found", self.tree_graph.shape, "measurements x parameters"

	def fit_func(self):
		return np.mean( (self.HI_dist - np.dot(self.tree_graph, self.params))**2 )

	def fit_l1reg(self):
		from l1regls import l1regls
		A = matrix(self.tree_graph)
		b = matrix(self.HI_dist)
		return np.array([x for x in l1regls(A/np.sqrt(self.lam),b/np.sqrt(self.lam))])

	def fit_nnls(self):
		from scipy.optimize import nnls
		return nnls(self.tree_graph, self.HI_dist)[0]

	def fit_nnl2reg(self):
		from cvxopt import matrix, solvers
		n_params = self.tree_graph.shape[1]
		P = matrix(np.dot(self.tree_graph.T, self.tree_graph) + self.lam*np.eye(n_params))
		q = matrix( -np.dot( self.HI_dist, self.tree_graph))
		h = matrix(np.zeros(n_params)) # Gw <=h
		G = matrix(-np.eye(n_params))
		W = solvers.qp(P,q,G,h)
		return np.array([x for x in W['x']])

	def fit_nnl1reg(self):
		from cvxopt import matrix, solvers
		n_params = self.tree_graph.shape[1]
		HI_sc = self.HI_split_count
		P1 = np.zeros((n_params+HI_sc,n_params+HI_sc))
		P1[:n_params, :n_params] = np.dot(self.tree_graph.T, self.tree_graph)
		for ii in xrange(HI_sc, n_params):
			P1[ii,ii]+=self.lam
		P = matrix(P1)

		q1 = np.zeros(n_params+HI_sc)
		q1[:n_params] = -np.dot( self.HI_dist, self.tree_graph)
		q1[n_params:] = self.lam
		q = matrix(q1)

		h = matrix(np.zeros(2*HI_sc)) 	# Gw <=h
		G1 = np.zeros((2*HI_sc,n_params+HI_sc))
		G1[:HI_sc, :HI_sc] = -np.eye(HI_sc)
		G1[HI_sc:, :HI_sc] = np.eye(HI_sc)
		G1[HI_sc:, n_params:] = -np.eye(HI_sc)
		G = matrix(G1)
		W = solvers.qp(P,q,G,h)
		return np.array([x for x in W['x']])[:n_params]

	def map_HI_to_tree(self, training_fraction = 1.0, method = 'nnls', lam=10):
		self.lam = lam
		if training_fraction<1.0:
			self.test_HI, self.train_HI = {}, {}
			for key, val in self.HI_normalized.iteritems():
				if np.random.uniform()>training_fraction:
					self.test_HI[key]=val
				else:
					self.train_HI[key]=val
		else:
			self.train_HI = self.HI_normalized

		self.make_treegraph()
		if method=='l1reg':  # l1 regularized fit, no constraint on sign of effect
			self.params = self.fit_l1reg()
		elif method=='nnls':  # non-negative least square, not regularized
			self.params = self.fit_nnls()
		elif method=='nnl2reg':	# non-negative L2 norm regularized fit
			self.params = self.fit_nnl2reg()
		elif method=='nnl1reg':  # non-negative fit, branch terms L1 regularized, avidity terms L2 regularized
			self.params = self.fit_nnl1reg()

		print "method",method, "regularized by", self.lam, "squared deviation=",self.fit_func()
		# for each set of branches with HI constraints, pick the branch with most aa mutations
		# and assign the dHI to that one, record the number of constraints
		for HI_split, branches in self.HI_split_to_branch.iteritems():
			likely_branch = branches[np.argmax([len(b.mutations) for b in branches])]
			likely_branch.dHI = self.params[HI_split]
			likely_branch.constraints = self.tree_graph[:,HI_split].sum()

		# integrate the HI change dHI into a cumulative antigentic evolution score cHI
		for node in self.tree.preorder_node_iter():
			if node!=self.tree.seed_node:
				node.cHI = node.parent_node.cHI + node.dHI
		self.avidities = {serum:self.params[self.HI_split_count+ii] for ii, serum in enumerate(self.sera)}

	def validate(self, plot=False):
		self.validation = {}
		for key, val in self.test_HI.iteritems():
			pred_HI = self.predict_HI(key[0], key[1])
			if pred_HI is not None:
				self.validation[key] = (val, pred_HI)
		if plot:
			import matplotlib.pyplot as plt
			a = np.array(self.validation.values())
			plt.figure()
			plt.scatter(a[:,0], a[:,1])
			plt.xlabel("measured")
			plt.ylabel("predicted")
			plt.xlabel("measured")
			plt.title('reg='+str(reg)+', avg error '+str(round(np.mean(np.abs(a[:,0]-a[:,1])),3)))


	def predict_HI(self, virus, serum):
		path = self.get_path(virus,serum)
		if path is not None:
			return self.avidities[serum] + np.sum(b.dHI for b in path)
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


def main(tree_fname = 'data/tree_ancestral.json', HI_fname='data/HI_titers.txt'):

	print "--- Fitting HI titers at " + time.strftime("%H:%M:%S") + " ---"

	tree =  json_to_dendropy(read_json(tree_fname))
	HI_distances = load_HI_distances(HI_fname)


if __name__ == "__main__":
	from Bio import Phylo
	reg = 10
	#tree_fname = 'data/tree_long_HI.json'
	tree_fname = 'data/tree_refine.json'
	tree =  json_to_dendropy(read_json(tree_fname))
	names, measurements, tables = read_tables()
	HI_map = HI_tree(tree, measurements)
	HI_map.map_HI_to_tree(training_fraction=0.9, method = 'nnl1reg', lam=reg)
	HI_map.validate(plot=True)

	out_tree_fname = 'data/tree_HI.json'
	write_json(dendropy_to_json(tree.seed_node), out_tree_fname, indent=None)

	#plot_tree(tree)
	#plot_dHI_distribution(tree)

