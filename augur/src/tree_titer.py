# clean, reroot, ladderize newick tree
# output to tree.json
import numpy as np
import time, os, gzip
from collections import defaultdict
from matplotlib import pyplot as plt
from itertools import izip

def myopen(fname, mode='r'):
	if fname[-2:]=='gz':
		return gzip.open(fname, mode)
	else:
		return open(fname, mode)

class HI_tree(object):

	def __init__(self, HI_fname = 'source-data/HI_titers.txt', min_aamuts = 0,**kwargs):
		self.HI, tmp = self.read_HI_titers(HI_fname)
		self.tree_graph = None
		self.min_aamuts = min_aamuts

	def read_HI_titers(self, fname):
		strains = set()
		measurements = {}
		with myopen(fname, 'r') as infile:
			for line in infile:
				entries = line.strip().split()
				test, ref_virus, serum, val = entries[0], entries[1],entries[2], map(float, entries[3:])
				ref = (ref_virus, serum)
				if len(val):
					measurements[(test, (ref_virus, serum))] = val
					strains.update([test, ref])
				else:
					print line.strip()
		return measurements, strains

	def normalize_HI(self):
		'''
		convert the HI measurements into the log2 difference between the average 
		HI titer measured between test virus and reference serum and the average 
		homologous titer. all measurements relative to sera without homologous titer
		are excluded
		'''
		consensus_func = np.mean
		self.HI_normalized = {}
		sera = set()
		ref_strains = set()
		HI_strains = set()
		for (test, ref), val in self.HI.iteritems():
			if test.upper() in self.node_lookup and ref[0].upper() in self.node_lookup:
				HI_strains.add(test.upper())
				HI_strains.add(ref[0].upper())
				if (ref[0],ref) in self.HI:
					sera.add(ref)
					ref_strains.add(ref[0])
					normalized_val = consensus_func(np.log2(self.HI[(ref[0], ref)])) - consensus_func(np.log2(val))
					self.HI_normalized[(test, ref)] = normalized_val
				else:
					print "no homologous titer found:", ref
		self.sera = list(sera)
		self.ref_strains = list(ref_strains)
		self.HI_strains = list(HI_strains)

	def add_mutations(self):
		'''
		add amino acid mutations to the tree
		'''
		self.tree.seed_node.mutations= ''
		for node in self.tree.postorder_node_iter():
			if node is not self.tree.seed_node:
				node.mutations = [a+str(pos+1)+b for pos, (a,b) in 
								enumerate(izip(node.parent_node.aa_seq, node.aa_seq)) if a!=b]

	def mark_HI_splits(self):
		# flag all branches on the tree with HI_strain = True if they lead to strain with titer data
		for leaf in self.tree.leaf_iter():
			if leaf.strain.upper() in self.HI_strains:
				leaf.serum = leaf.strain.upper() in self.ref_strains
				leaf.HI_info= 1
			else:
				leaf.serum, leaf.HI_info=False, 0

		for node in self.tree.postorder_internal_node_iter():
			node.HI_info = sum([c.HI_info for c in node.child_nodes()])
			node.serum= False

		# combine sets of branches that span identical sets of HI measurements
		self.HI_split_count = 0  # HI measurment split counter
		self.HI_split_to_branch = defaultdict(list)
		for node in self.tree.preorder_node_iter():
			node.dHI, node.cHI, node.constraints =0, 0, 0
			if node.HI_info>1:
				node.HI_branch_index = self.HI_split_count
				self.HI_split_to_branch[node.HI_branch_index].append(node)
				if sum([c.HI_info>0 for c in node.child_nodes()])>1:
					self.HI_split_count+=1
				elif node.is_leaf():
					self.HI_split_count+=1

		print "# of reference strains:",len(self.sera), "# of branches with HI constraint", self.HI_split_count

	def get_path_no_terminals(self, v1, v2):
		'''
		returns the path between two tips in the tree excluding the terminal branches. 
		'''
		if v1 in self.node_lookup and v2 in self.node_lookup:
			p1 = [self.node_lookup[v1]]
			p2 = [self.node_lookup[v2]]
			for tmp_p in [p1,p2]:
				while tmp_p[-1].parent_node != self.tree.seed_node:
					tmp_p.append(tmp_p[-1].parent_node)
				tmp_p.append(self.tree.seed_node)
				tmp_p.reverse()

			for pi, (tmp_v1, tmp_v2) in enumerate(izip(p1,p2)):
				if tmp_v1!=tmp_v2:
					break
			path = [n for n in p1[pi:] if n.HI_info>1] + [n for n in p2[pi:] if n.HI_info>1]
		else:
			path = None
		return path

	def make_treegraph(self):
		'''
		code the path between serum and test virus of each HI measurement into a matrix
		the matrix has dimensions #measurements x #tree branches with HI info
		if the path between test and serum goes through a branch, the corresponding matrix element is 1, 0 otherwise
		'''
		tree_graph = []
		HI_dist = []
		n_params = self.HI_split_count + len(self.sera) + len(self.HI_strains)
		for (test, ref), val in self.train_HI.iteritems():
			if not np.isnan(val):
				try:
					if ref[0] in self.node_lookup and test in self.node_lookup\
						and self.node_lookup[ref[0]].parent_node is not None\
						and self.node_lookup[test].parent_node is not None:
						path = self.get_path_no_terminals(test, ref[0])
						tmp = np.zeros(n_params)
						# determine branch indices on path
						if type(self.min_aamuts)==int:
							branches = np.unique([c.HI_branch_index for c in path 
							                     if (hasattr(c, 'HI_branch_index') and len(c.aa_muts)>=self.min_aamuts)])
						elif self.min_aamuts=='epi':
							branches = np.unique([c.HI_branch_index for c in path if (hasattr(c, 'HI_branch_index') and c.parent_node.ep<c.ep)])
						elif self.min_aamuts=='rbs':
							branches = np.unique([c.HI_branch_index for c in path if (hasattr(c, 'HI_branch_index') and c.parent_node.rb<c.rb)])
						else:
							branches = np.unique([c.HI_branch_index for c in path if hasattr(c, 'HI_branch_index') ])
						if len(branches): tmp[branches] = 1
						# add serum effect
						tmp[self.HI_split_count+self.sera.index(ref)] = 1
						# add virus effect
						tmp[self.HI_split_count+len(self.sera)+self.HI_strains.index(test)] = 1
						# append model and fit value to lists tree_graph and HI_dist
						tree_graph.append(tmp)
						HI_dist.append(val)
				except:
					import pdb; pdb.set_trace()
					print test, ref, "ERROR"

		# convert to numpy arrays and save product of tree graph with its transpose for future use
		self.HI_dist =  np.array(HI_dist)
		self.tree_graph = np.array(tree_graph)
		self.TgT = np.dot(self.tree_graph.T, self.tree_graph)
		print "Found", self.tree_graph.shape, "measurements x parameters"

	def fit_func(self):
		return np.mean( (self.HI_dist - np.dot(self.tree_graph, self.params))**2 )

	def fit_l1reg(self):
		from cvxopt import matrix
		from l1regls import l1regls
		A = matrix(self.tree_graph)
		b = matrix(self.HI_dist)
		return np.array([x for x in l1regls(A/np.sqrt(self.lam_HI),b/np.sqrt(self.lam_HI))])

	def fit_nnls(self):
		from scipy.optimize import nnls
		return nnls(self.tree_graph, self.HI_dist)[0]

	def fit_nnl2reg(self):
		from cvxopt import matrix, solvers
		n_params = self.tree_graph.shape[1]
		P = matrix(np.dot(self.tree_graph.T, self.tree_graph) + self.lam_HI*np.eye(n_params))
		q = matrix( -np.dot( self.HI_dist, self.tree_graph))
		h = matrix(np.zeros(n_params)) # Gw <=h
		G = matrix(-np.eye(n_params))
		W = solvers.qp(P,q,G,h)
		return np.array([x for x in W['x']])

	def fit_nnl1reg(self):
		from cvxopt import matrix, solvers
		n_params = self.tree_graph.shape[1]
		HI_sc = self.HI_split_count
		n_sera = len(self.sera)
		n_v = len(self.HI_strains)

		# set up the quadratic matrix containing the deviation term (linear xterm below)
		# and the l2-regulatization of the avidities and potencies
		P1 = np.zeros((n_params+HI_sc,n_params+HI_sc))
		P1[:n_params, :n_params] = self.TgT
		for ii in xrange(HI_sc, HI_sc+n_sera):
			P1[ii,ii]+=self.lam_pot*0.5
		for ii in xrange(HI_sc+n_sera, n_params):
			P1[ii,ii]+=self.lam_avi*0.5
		P = matrix(P1)

		# set up cost for auxillary parameter and the linear cross-term
		q1 = np.zeros(n_params+HI_sc)
		q1[:n_params] = -np.dot( self.HI_dist, self.tree_graph)
		q1[n_params:] = self.lam_HI
		q = matrix(q1)

		# set up linear constraint matrix to enforce positivity of the
		# dHIs and bounding of dHI by the auxillary parameter
		h = matrix(np.zeros(2*HI_sc)) 	# Gw <=h
		G1 = np.zeros((2*HI_sc,n_params+HI_sc))
		G1[:HI_sc, :HI_sc] = -np.eye(HI_sc)
		G1[HI_sc:, :HI_sc] = np.eye(HI_sc)
		G1[HI_sc:, n_params:] = -np.eye(HI_sc)
		G = matrix(G1)
		W = solvers.qp(P,q,G,h)
		sol = np.array([x for x in W['x']])[:n_params]
		self.params=sol
		print "squared deviation prior to relax=",self.fit_func()
		# redo the linear cost relaxing terms that seem to be relevant to avoid 
		# compression of the fit. 0.2 seems to be a good cut-off, linear tune to zero
		q1[n_params:] = self.lam_HI*(1-5.0*np.minimum(0.2,sol[:HI_sc]))
		q = matrix(q1)
		W = solvers.qp(P,q,G,h)
		sol = np.array([x for x in W['x']])[:n_params]
		self.params=sol
		print "squared deviation after relax=",self.fit_func()
		return sol

	def prepare_HI_map(self):
		'''
		normalize the HI measurements, split the data into training and test sets
		and determine which branches on the tree are transversed by HI measurements
		'''
		from random import sample
		self.normalize_HI()
		self.add_mutations()
		self.mark_HI_splits()
		if self.training_fraction<1.0:
			self.test_HI, self.train_HI = {}, {}
			if self.subset_strains:
				training_strains = sample(self.HI_strains, int(self.training_fraction*len(self.HI_strains)))
				for tmpstrain in self.ref_strains:
					if tmpstrain not in training_strains:
						training_strains.append(tmpstrain)
				for key, val in self.HI_normalized.iteritems():
					if key[0] in training_strains:
						self.train_HI[key]=val
					else:
						self.test_HI[key]=val
			else:
				for key, val in self.HI_normalized.iteritems():
					if np.random.uniform()>self.training_fraction:
						self.test_HI[key]=val
					else:
						self.train_HI[key]=val
		else:
			self.train_HI = self.HI_normalized

		if self.cutoff_date is not None:
			self.train_HI = {key:val for key,val in self.train_HI.iteritems()
							if self.node_lookup[key[0]].num_date<cutoff_date and 
							   self.node_lookup[key[1][0]].num_date<cutoff_date}

		self.make_treegraph()		

	def map_HI_to_tree(self, training_fraction = 1.0, method = 'nnls', lam_HI=5.0, 
						lam_pot = 5.0, lam_avi = 5.0, cutoff_date = None, subset_strains = False):
		self.training_fraction = training_fraction
		self.subset_strains=subset_strains
		self.lam_pot = lam_pot
		self.lam_avi = lam_avi
		self.lam_HI = lam_HI
		self.cutoff_date = cutoff_date
		if self.tree_graph is None:
			self.prepare_HI_map()

		if method=='l1reg':  # l1 regularized fit, no constraint on sign of effect
			self.params = self.fit_l1reg()
		elif method=='nnls':  # non-negative least square, not regularized
			self.params = self.fit_nnls()
		elif method=='nnl2reg':	# non-negative L2 norm regularized fit
			self.params = self.fit_nnl2reg()
		elif method=='nnl1reg':  # non-negative fit, branch terms L1 regularized, avidity terms L2 regularized
			self.params = self.fit_nnl1reg()

		print "method",method, "regularized by", self.lam_HI, "squared deviation=",self.fit_func()
		# for each set of branches with HI constraints, pick the branch with most aa mutations
		# and assign the dHI to that one, record the number of constraints
		for node in self.tree.postorder_node_iter():
			node.dHI=0
		for HI_split, branches in self.HI_split_to_branch.iteritems():
			likely_branch = branches[np.argmax([len(b.mutations) for b in branches])]
			likely_branch.dHI = self.params[HI_split]
			likely_branch.constraints = self.tree_graph[:,HI_split].sum()

		# integrate the HI change dHI into a cumulative antigentic evolution score cHI
		for node in self.tree.preorder_node_iter():
			if node!=self.tree.seed_node:
				node.cHI = node.parent_node.cHI + node.dHI
		self.serum_potency = {serum:self.params[self.HI_split_count+ii] 
							  for ii, serum in enumerate(self.sera)}
		self.virus_effect = {strain:self.params[self.HI_split_count+len(self.sera)+ii]
							  for ii, strain in enumerate(self.HI_strains)}

	def validate(self, plot=False):
		self.validation = {}
		for key, val in self.test_HI.iteritems():
			pred_HI = self.predict_HI(key[0], key[1])
			if pred_HI is not None:
				self.validation[key] = (val, pred_HI)
		from scipy.stats import linregress, pearsonr
		a = np.array(self.validation.values())
		self.abs_error = np.mean(np.abs(a[:,0]-a[:,1]))
		self.rms_error = np.sqrt(np.mean((a[:,0]-a[:,1])**2))
		self.slope, self.intercept, tmpa, tmpb, tmpc = linregress(a[:,0], a[:,1])
		print "error (abs/rms): ",self.abs_error, self.rms_error
		print "slope, intercept:", self.slope, self.intercept
		print "pearson correlation:", pearsonr(a[:,0], a[:,1])
		if plot:
			import matplotlib.pyplot as plt
			plt.figure()
			plt.plot([-1,6], [-1,6], 'k')
			plt.scatter(a[:,0], a[:,1])
			plt.xlabel("measured")
			plt.ylabel("predicted")
			plt.xlabel("measured")
			plt.title('reg HI/pot/avi ='+str(self.lam_HI)+'/'+str(self.lam_pot)+'/'+str(self.lam_avi)+', avg abs/rms '\
						+str(round(self.abs_error, 3))\
					+'/'+str(round(self.rms_error,3)))

	def add_titers(self):
		for ref in self.ref_strains:
			self.node_lookup[ref].HI_titers_perserum= defaultdict(dict)
			self.node_lookup[ref].potency={}
		for ref in self.sera:
			self.node_lookup[ref[0]].potency[ref[1]] = self.serum_potency[ref]
		for (test, ref), val in self.HI_normalized.iteritems():
			self.node_lookup[ref[0]].HI_titers_perserum[self.node_lookup[test].clade][ref[1]] = val
		for test in self.HI_strains:
			self.node_lookup[test].avidity = self.virus_effect[test]
		for ref in self.ref_strains:
			self.node_lookup[ref].HI_titers = {key:np.mean(titers.values()) for key, titers in 
			 									self.node_lookup[ref].HI_titers_perserum.iteritems()}
			self.node_lookup[ref].mean_potency = np.mean(self.node_lookup[ref].potency.values())




	def predict_HI(self, virus, serum):
		path = self.get_path_no_terminals(virus,serum[0])
		if path is not None:
			return self.serum_potency[serum] + self.virus_effect[virus] + np.sum(b.dHI for b in path)
		else:
			return None

######################################################################################
####  utility functions for reading and writing HI Tables, plotting trees, etc
######################################################################################

def plot_tree(tree):
	from Bio import Phylo
	from tree_util import to_Biopython, color_BioTree_by_attribute
	btree = to_Biopython(tree)
	color_BioTree_by_attribute(btree,"cHI", transform = lambda x:x)
	Phylo.draw(btree, label_func = lambda  x: 'X' if x.serum else '' if x.HI_info else '', 
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

min_titer = 10.0

def strain_name_fixing(name):
	return name.replace(' ','').upper().strip().strip('*').lstrip('0123456789')
def titer_to_number(val):
	try:
		if '<' in val:
			return np.nan
		if len(val.split())>1:
			return float(val.split()[0])
		else:
			return float(val)
	except:
		#print "Bad HI measurement:", val
		return np.nan

def parse_HI_matrix(fname):
	import pandas as pd
	from string import strip
	import csv
	name_abbrev = {'HK':"HONGKONG", 'SWITZ':"SWITZERLAND", 'VIC':"VICTORIA", 'STOCK':"STOCKHOLM",
					'STHAFR':"SOUTHAFRICA", 'SAFRICA':"SOUTHAFRICA", "ENG":"ENGLAND", "NIB-85":"A/ALMATY/2958/2013", 'NOR':'NORWAY',
					'NTHCAROL':"NORTHCAROLINA",'ALA':"ALABAMA", 'NY':"NEWYORK", "GLAS":"GLASGOW", "AL":"ALABAMA",
					"NETH":"NETHERLANDS", "FIN":"FINLAND", "BRIS":"BRISBANE", "MARY":"MARYLAND",	
					"ST.P'BURG":"ST.PETERSBURG", 'CAL':'CALIFORNIA', 'AUCK':'AUCKLAND', "C'CHURCH":'CHRISTCHURCH',
					'CHCH':'CHRISTCHURCH', 'ASTR':'ASTRAKHAN', 'ASTRAK':'ASTRAKHAN', 'ST.P':"ST.PETERSBURG",
					'JHB':'JOHANNESBURG', 'FOR':'FORMOSA','MAL':'MALAYSIA', 'STHAUS':'SOUTHAUSTRALIA',
					'FL':'FLORIDA', 'MASS':'MASSACHUSETTS','NOVO':'NOVOSIBIRSK','WIS':'WISCONSIN','BANG':'BANGLADESH','EG':'EGYPT' 	}
	src_id = fname.split('/')[-1]
	print fname
	with myopen(fname) as infile:
		csv_reader = csv.reader(infile)

		# parse sera
		row1 = csv_reader.next()
		row2 = csv_reader.next()
		row3 = csv_reader.next()
		ref_sera = [[strain_name_fixing(e1+'/'+e2), e3.replace(' ','')] for e1,e2,e3 in zip(row1, row2, row3)[4:]]
		for ri in xrange(len(ref_sera)):
			abbr = ref_sera[ri][0].split('/')[1].rstrip('01234566789')
			if abbr in name_abbrev:
				ref_sera[ri][0] = strain_name_fixing(ref_sera[ri][0].replace(abbr, name_abbrev[abbr]))
			else:
				ref_sera[ri][0] = strain_name_fixing(ref_sera[ri][0])
			# strip numbers
			tmp = ref_sera[ri][0].split('/')
			ref_sera[ri][0] = '/'.join([tmp[0], tmp[1].rstrip('0123456789')]+tmp[2:])
			try:
				y = int(ref_sera[ri][0].split('/')[-1])
				if y<100:
					if y<20:
						ref_sera[ri][0] = '/'.join(ref_sera[ri][0].split('/')[:-1])+'/'+str(2000+y)
					else:
						ref_sera[ri][0] = '/'.join(ref_sera[ri][0].split('/')[:-1])+'/'+str(1900+y)
			except:
				print ref_sera[ri]

		fields = ['source','ref/test', 'genetic group', 'collection date', 'passage history']+map(tuple, ref_sera)
		#print fields
		for row in csv_reader: # advance until the reference virus
			if row[0].startswith('REFERENCE'):
				break

		ref_strains = []
		ref_matrix = []
		for row in csv_reader: 
			if row[0].startswith('TEST'):
				break
			else: # load matrices until the test virus section starts
				ref_strains.append(strain_name_fixing(row[0].strip()))
				ref_matrix.append([src_id,'ref']+map(strip, row[1:4])+map(titer_to_number, row[4:]))

		test_strains = []
		test_matrix = []
		for row in csv_reader: # load test viruses until it is no longer an A/ flu  name
			if not row[0].startswith('A/'):
				break
			else:
				test_strains.append(strain_name_fixing(row[0].strip()))
				test_matrix.append([src_id,'test']+map(strip,row[1:4])+map(titer_to_number, row[4:]))

		print len(ref_sera), ref_sera
		print len(ref_strains), len(test_strains)
		HI_table  = pd.DataFrame(ref_matrix+test_matrix, index = ref_strains+test_strains, columns= fields)

		return HI_table


def read_tables(flutype = 'H3N2'):
	import glob
	import pandas as pd
	from itertools import product
	flist = glob.glob('../../flu_HI_data/'+flutype+'_tables/NIMR*csv')
	all_names = set()
	all_measurements = defaultdict(list)
	HI_matrices = pd.DataFrame()
	for fname in flist:
		tmp = parse_HI_matrix(fname)
		HI_matrices = HI_matrices.append(tmp)
	return HI_matrices

def read_trevor_table(flutype):
	trevor_table = 'source-data/'+flutype+'_HI.tsv'
	import csv
	measurements = defaultdict(list)
	sera = set()
	strains = set()
	if os.path.isfile(trevor_table):
		with myopen(trevor_table) as infile:
			table_reader = csv.reader(infile, delimiter="\t")
			header = table_reader.next()
			for row in table_reader:
	#			try:
					val = titer_to_number(row[6])
					if not np.isnan(val):
						strains.add(strain_name_fixing(row[1]))
						serum = (strain_name_fixing(row[4]), row[3])
						sera.add(serum)
						measurements[(strain_name_fixing(row[1]), serum)].append(val)
	#			except:
	#				print row
	else:
		print trevor_table, "not found"
	print "trevor total:", len(measurements), "measurements"
	return strains, sera, measurements


def table_to_flat(HI_table):
	flat_measurements = defaultdict(list)
	for ref_serum in HI_table.columns[5:]:
		sub_set = HI_table[ref_serum][~np.isnan(HI_table[ref_serum])]
		for virus, val in izip(sub_set.index, sub_set):
			flat_measurements[(virus, ref_serum)].append(val)
	print "NIMR total:", len(flat_measurements), "measurements"
	return flat_measurements

def get_all_titers_flat(flutype='H3N2'):
	HI_titers = read_tables(flutype)
	HI_titers_flat = table_to_flat(HI_titers)
	HI_trevor = read_trevor_table(flutype)[2]
	HI_titers_flat.update(HI_trevor)
	return HI_titers_flat


def write_strains_with_HI_and_sequence(flutype='H3N2'):
	HI_titers = read_tables(flutype)
	HI_trevor = read_trevor_table(flutype)
	HI_strains = set(HI_titers.index)
	HI_strains.update([v[0] for v in HI_trevor[2]])
	from Bio import SeqIO
	good_strains = set()
	with myopen("data/"+flutype+"_strains_with_HI.fasta", 'w') as outfile, \
		 myopen("source-data/"+flutype+"_HI_strains.txt", 'w') as HI_strain_outfile, \
		 myopen("data/"+flutype+"_gisaid_epiflu_sequence.fasta", 'r') as infile:
		for seq_rec in SeqIO.parse(infile, 'fasta'):
			tmp_name = seq_rec.description.split('|')[0].strip()
			reduced_name = strain_name_fixing(tmp_name)
			if reduced_name in HI_strains and (reduced_name not in good_strains):
				SeqIO.write(seq_rec, outfile,'fasta')
				good_strains.add(reduced_name)
				HI_strain_outfile.write(tmp_name+'\n')
				#print seq_rec.name


def write_flat_HI_titers(flutype = 'H3N2', fname = None):
	measurements = get_all_titers_flat(flutype)
	with myopen('source-data/'+flutype+'_HI_strains.txt') as infile:
		strains = [strain_name_fixing(line.strip()) for line in infile]	
	if fname is None:
		fname = 'source-data/'+flutype+'_HI_titers.txt'
	with myopen(fname, 'w') as outfile:
		for (test, ref), val in measurements.iteritems():
			if strain_name_fixing(test) in strains and strain_name_fixing(ref[0]) in strains:
				outfile.write(test+'\t'+ref[0]+'\t'+ref[1]+'\t'+'\t'.join(map(str,val))+'\n')


def main(tree, HI_fname='source-data/HI_titers.txt', training_fraction = 1.0, reg=5):
	print "--- Fitting HI titers at " + time.strftime("%H:%M:%S") + " ---"
	measurements, strains = read_HI_titers(HI_fname)
	HI_map = HI_tree(tree, measurements)
	HI_map.map_HI_to_tree(training_fraction=training_fraction, method = 'nnl1reg', lam=reg)
	return HI_map


if __name__ == "__main__":
	from Bio import Phylo
	from io_util import json_to_dendropy
	reg = 10
	tree_fname = 'data/tree_refine.json'
	tree =  json_to_dendropy(read_json(tree_fname))
	HI_map = main(tree, training_fraction=0.8, reg=reg)
	HI_map.validate(plot=True)
	plot_tree(tree)
	plot_dHI_distribution(tree)

	out_tree_fname = 'data/tree_HI.json'
	write_json(dendropy_to_json(tree.seed_node), out_tree_fname, indent=None)
