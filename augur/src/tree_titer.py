# clean, reroot, ladderize newick tree
# output to tree.json
import numpy as np
import time, os, gzip
from collections import defaultdict
from matplotlib import pyplot as plt
from itertools import izip
from virus_filter import fix_name
import pandas as pd
from diagnostic_figures import fs, fmts, figheight

def myopen(fname, mode='r'):
	if fname[-2:]=='gz':
		return gzip.open(fname, mode)
	else:
		return open(fname, mode)

def HI_fix_name(name):
	tmp_name = fix_name(name)
	return tmp_name.upper().lstrip('*')


class HI_tree(object):

	def __init__(self, HI_fname = 'source-data/HI_titers.txt', min_aamuts = 0,**kwargs):
		self.HI_fname = HI_fname
		if "excluded_tables" in kwargs:
			self.excluded_tables = kwargs["excluded_tables"]
		else:
			self.excluded_tables = []
		self.HI, tmp, sources = self.read_HI_titers(HI_fname)
		self.sources = list(sources)
		self.tree_graph = None
		self.min_aamuts = min_aamuts
		self.serum_potency = {}
		self.virus_effect = {}

	def read_HI_titers(self, fname):
		strains = set()
		measurements = defaultdict(list)
		sources = set()
		with myopen(fname, 'r') as infile:
			for line in infile:
				entries = line.strip().split()
				test, ref_virus, serum, src_id, val = (entries[0], entries[1],entries[2],
														entries[3], float(entries[4]))
				ref = (ref_virus, serum)
				if src_id not in self.excluded_tables:
					try:
						measurements[(test, (ref_virus, serum))].append(val)
						strains.update([test, ref])
						sources.add(src_id)
					except:
						print line.strip()
		return measurements, strains, sources

	def normalize(self, ref, val):
		consensus_func = np.mean
		return consensus_func(np.log2(self.autologous_titers[ref]['val'])) - consensus_func(np.log2(val))

	def normalize_HI(self):
		'''
		convert the HI measurements into the log2 difference between the average
		HI titer measured between test virus and reference serum and the average
		homologous titer. all measurements relative to sera without homologous titer
		are excluded
		'''
		self.HI_normalized = {}
		self.HI_raw = {}
		sera = set()
		ref_strains = set()
		HI_strains = set()
		all_per_serum = defaultdict(list)
		autologous = defaultdict(list)
		for (test, ref), val in self.HI.iteritems():
			if test.upper() in self.node_lookup and ref[0].upper() in self.node_lookup:
				HI_strains.add(test.upper())
				HI_strains.add(ref[0].upper())
				all_per_serum[ref].append(val)
				if ref[0]==test:
					autologous[ref].append(val)

		self.autologous_titers = {}
		for serum in all_per_serum:
			if serum in autologous:
				self.autologous_titers[serum] = {'val':autologous[serum], 'autologous':True}
				print("autologous titer found for",serum)
			else:
				if len(all_per_serum[serum])>10:
					self.autologous_titers[serum] = {'val':np.max(all_per_serum[serum]), 'autologous':False}
					print(serum,": using max titer instead of autologous,", np.max(all_per_serum[serum]))
				else:
					print("discarding",serum,"since there are only ",len(all_per_serum[serum]),'measurements')

		for (test, ref), val in self.HI.iteritems():
			if test.upper() in self.node_lookup and ref[0].upper() in self.node_lookup:
				if ref in self.autologous_titers:
					sera.add(ref)
					ref_strains.add(ref[0])
					self.HI_normalized[(test, ref)] = self.normalize(ref, val)
					self.HI_raw[(test, ref)] = np.median(val)
				else:
					pass
					#print "no homologous titer found:", ref

		self.sera = list(sera)
		self.ref_strains = list(ref_strains)
		self.HI_strains = list(HI_strains)

	def add_mutations(self):
		'''
		add amino acid mutations to the tree
		'''
		self.tree.seed_node.mutations= []
		for node in self.tree.postorder_node_iter():
			if node is not self.tree.seed_node:
				node.mutations = []
				for prot in node.aa_muts:
					if node.aa_muts[prot]!='':
						node.mutations.extend(node.aa_muts[prot].split(','))

	def mark_HI_splits(self):
		# flag all branches on the tree with HI_strain = True if they lead to strain with titer data
		for leaf in self.tree.leaf_iter():
			if leaf.strain.upper() in self.HI_strains:
				leaf.serum = leaf.strain.upper() in self.ref_strains
				leaf.HI_info = 1
			else:
				leaf.serum, leaf.HI_info=False, 0

		for node in self.tree.postorder_internal_node_iter():
			node.HI_info = sum([c.HI_info for c in node.child_nodes()])
			node.serum= False

		# combine sets of branches that span identical sets of HI measurements
		self.HI_split_count = 0  # HI measurment split counter
		self.HI_split_to_branch = defaultdict(list)
		for node in self.tree.preorder_node_iter():
			if self.map_to_tree:
				node.dHI, node.cHI, node.constraints = 0, 0, 0
			if node.HI_info>1:
				node.HI_branch_index = self.HI_split_count
				self.HI_split_to_branch[node.HI_branch_index].append(node)
				# at a bi- or multifurcation, increase the split count and HI index
				# either individual child branches have enough HI info be counted,
				# or the pre-order node iteraction will move towards the root
				if sum([c.HI_info>0 for c in node.child_nodes()])>1:
					self.HI_split_count+=1
				elif node.is_leaf():
					self.HI_split_count+=1

		self.genetic_params = self.HI_split_count
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

	def get_mutations(self, strain1, strain2):
		if strain1 in self.node_lookup and strain2 in self.node_lookup:
			node1 = self.node_lookup[strain1].parent_node
			node2 = self.node_lookup[strain2].parent_node
			if node1 is None or node2 is None:
				return None
			else:
				return self.get_mutations_nodes(node1, node2)
		else:
			return None

	def get_mutations_nodes(self, node1, node2):
		muts = []
		for prot in node1.aa_seq:
			seq1 = node1.aa_seq[prot]
			seq2 = node2.aa_seq[prot]
			muts.extend([(prot, aa1+str(pos+1)+aa2) for pos, (aa1, aa2) in enumerate(izip(seq1, seq2)) if aa1!=aa2])
		return muts

	def make_seqgraph(self, colin_thres = 5):
		'''
		code amino acid differences between sequences into a matrix
		the matrix has dimensions #measurements x #observed mutations
		'''
		seq_graph = []
		HI_dist = []
		# list all mutations
		mutation_counter = defaultdict(int)
		for (test, ref), val in self.train_HI.iteritems():
			muts = self.get_mutations(ref[0], test)
			if muts is None:
				continue
			for mut in muts:
				mutation_counter[mut]+=1

		relevant_muts = []
		min_count= 10
		min_freq = 1.0*min_count/len(self.viruses)
		for mut, count in mutation_counter.iteritems():
			gene = mut[0]
			pos = int(mut[1][1:-1])-1
			aa1, aa2 = mut[1][0],mut[1][-1]
			if gene=='HA1' and count>min_count and \
				self.aa_frequencies[gene][self.aa_alphabet.index(aa1),pos]>min_freq and\
				self.aa_frequencies[gene][self.aa_alphabet.index(aa2),pos]>min_freq:
				relevant_muts.append(mut)

		relevant_muts.sort(key = lambda x:int(x[1][1:-1]))

		self.relevant_muts = relevant_muts
		self.genetic_params = len(relevant_muts)
		n_params = self.genetic_params + len(self.sera) + len(self.HI_strains)

		for (test, ref), val in self.train_HI.iteritems():
			if not np.isnan(val):
				try:
					muts = self.get_mutations(ref[0], test)
					if muts is None:
						continue
					tmp = np.zeros(n_params)
					# determine branch indices on path
					branches = np.unique([relevant_muts.index(mut) for mut in muts
					                     if mut in relevant_muts])
					if len(branches): tmp[branches] = 1
					# add serum effect
					tmp[len(relevant_muts)+self.sera.index(ref)] = 1
					# add virus effect
					tmp[len(relevant_muts)+len(self.sera)+self.HI_strains.index(test)] = 1
					# append model and fit value to lists seq_graph and HI_dist
					seq_graph.append(tmp)
					HI_dist.append(val)
				except:
					import pdb; pdb.set_trace()
					print test, ref, "ERROR"

		# convert to numpy arrays and save product of tree graph with its transpose for future use
		self.HI_dist =  np.array(HI_dist)
		self.tree_graph = np.array(seq_graph)
		if colin_thres is not None:
			self.collapse_colinear_mutations(colin_thres)
		self.TgT = np.dot(self.tree_graph.T, self.tree_graph)
		print "Found", self.tree_graph.shape, "measurements x parameters"

	def collapse_colinear_mutations(self, colin_thres):
		'''
		find colinear columns of the design matrix, collapse them into clusters
		'''
		TT = self.tree_graph[:,:self.genetic_params].T
		mutation_clusters = []
		n_measurements = self.tree_graph.shape[0]
		for col, mut in izip(TT, self.relevant_muts):
			col_found = False
			for cluster in mutation_clusters:
				if np.sum(col==cluster[0])>=n_measurements-colin_thres:
					cluster[1].append(mut)
					col_found=True
					print("adding",mut,"to cluster ",cluster[1])
					break
			if not col_found:
				mutation_clusters.append([col, [mut]])
		print("dimensions of old design matrix",self.tree_graph.shape)
		self.tree_graph = np.hstack((np.array([c[0] for c in mutation_clusters]).T,
									 self.tree_graph[:,self.genetic_params:]))
		self.genetic_params = len(mutation_clusters)
		# use the first mutation of a cluster to index the effect
		# make a dictionary that maps this effect to the cluster
		self.mutation_clusters = {c[1][0]:c[1] for c in mutation_clusters}
		self.relevant_muts = [c[1][0] for c in mutation_clusters]
		print("dimensions of new design matrix",self.tree_graph.shape)

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
							                     if (hasattr(c, 'HI_branch_index') and
							                         len(c.mutations)>=self.min_aamuts)])
						elif self.min_aamuts=='epi':
							branches = np.unique([c.HI_branch_index for c in path
							                     if (hasattr(c, 'HI_branch_index') and c.parent_node.ep<c.ep)])
						elif self.min_aamuts=='rbs':
							branches = np.unique([c.HI_branch_index for c in path
							                     if (hasattr(c, 'HI_branch_index') and c.parent_node.rb<c.rb)])
						else:
							branches = np.unique([c.HI_branch_index for c in path
							                     if hasattr(c, 'HI_branch_index') ])
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
		from cvxopt import matrix, solvers
		n_params = self.tree_graph.shape[1]
		HI_sc = self.genetic_params if self.map_to_tree else len(self.relevant_muts)
		n_sera = len(self.sera)
		n_v = len(self.HI_strains)

		# set up the quadratic matrix containing the deviation term (linear xterm below)
		# and the l2-regulatization of the avidities and potencies
		P1 = np.zeros((n_params+HI_sc,n_params+HI_sc))
		P1[:n_params, :n_params] = self.TgT
		for ii in xrange(HI_sc, HI_sc+n_sera):
			P1[ii,ii]+=self.lam_pot
		for ii in xrange(HI_sc+n_sera, n_params):
			P1[ii,ii]+=self.lam_avi
		P = matrix(P1)

		# set up cost for auxillary parameter and the linear cross-term
		q1 = np.zeros(n_params+HI_sc)
		q1[:n_params] = -np.dot( self.HI_dist, self.tree_graph)
		q1[n_params:] = self.lam_HI
		q = matrix(q1)

		# set up linear constraint matrix to regularize the HI parametesr
		h = matrix(np.zeros(2*HI_sc)) 	# Gw <=h
		G1 = np.zeros((2*HI_sc,n_params+HI_sc))
		G1[:HI_sc, :HI_sc] = -np.eye(HI_sc)
		G1[:HI_sc:, n_params:] = -np.eye(HI_sc)
		G1[HI_sc:, :HI_sc] = np.eye(HI_sc)
		G1[HI_sc:, n_params:] = -np.eye(HI_sc)
		G = matrix(G1)
		W = solvers.qp(P,q,G,h)
		self.params = np.array([x for x in W['x']])[:n_params]
		print "rms deviation prior to relax=",np.sqrt(self.fit_func())
		return self.params

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
		HI_sc = self.genetic_params if self.map_to_tree else len(self.relevant_muts)
		n_sera = len(self.sera)
		n_v = len(self.HI_strains)

		# set up the quadratic matrix containing the deviation term (linear xterm below)
		# and the l2-regulatization of the avidities and potencies
		P1 = np.zeros((n_params,n_params))
		P1[:n_params, :n_params] = self.TgT
		for ii in xrange(HI_sc, HI_sc+n_sera):
			P1[ii,ii]+=self.lam_pot
		for ii in xrange(HI_sc+n_sera, n_params):
			P1[ii,ii]+=self.lam_avi
		P = matrix(P1)

		# set up cost for auxillary parameter and the linear cross-term
		q1 = np.zeros(n_params)
		q1[:n_params] = -np.dot(self.HI_dist, self.tree_graph)
		q1[:HI_sc] += self.lam_HI
		q = matrix(q1)

		# set up linear constraint matrix to enforce positivity of the
		# dHIs and bounding of dHI by the auxillary parameter
		h = matrix(np.zeros(HI_sc)) 	# Gw <=h
		G1 = np.zeros((HI_sc,n_params))
		G1[:HI_sc, :HI_sc] = -np.eye(HI_sc)
		G = matrix(G1)
		W = solvers.qp(P,q,G,h)
		self.params = np.array([x for x in W['x']])[:n_params]
		print "rms deviation prior to relax=",np.sqrt(self.fit_func())
		# redo the linear cost relaxing terms that seem to be relevant to avoid
		# compression of the fit. 0.2 seems to be a good cut-off, linear tune to zero
		#q1[n_params:] = self.lam_HI*(1-5.0*np.minimum(0.2,sol[:HI_sc]))
		#q = matrix(q1)
		#W = solvers.qp(P,q,G,h)
		#sol = np.array([x for x in W['x']])[:n_params]
		#self.params=sol
		#print "rms deviation after relax=",np.sqrt(self.fit_func())
		return self.params

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
				tmp = set(self.HI_strains)
				tmp.difference_update(self.ref_strains)
				training_strains = sample(tmp, int(self.training_fraction*len(tmp)))
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

		if self.map_to_tree:
			self.make_treegraph()
		else:
			self.make_seqgraph()

	def map_HI(self, training_fraction = 1.0, method = 'nnls', lam_HI=1.0, map_to_tree = True,
			lam_pot = 0.5, lam_avi = 3.0, cutoff_date = None, subset_strains = False, force_redo = False):
		self.map_to_tree = map_to_tree
		self.training_fraction = training_fraction
		self.subset_strains=subset_strains
		self.lam_pot = lam_pot
		self.lam_avi = lam_avi
		self.lam_HI = lam_HI
		self.cutoff_date = cutoff_date
		if self.tree_graph is None or force_redo:
			self.prepare_HI_map()

		if method=='l1reg':  # l1 regularized fit, no constraint on sign of effect
			self.params = self.fit_l1reg()
		elif method=='nnls':  # non-negative least square, not regularized
			self.params = self.fit_nnls()
		elif method=='nnl2reg':	# non-negative L2 norm regularized fit
			self.params = self.fit_nnl2reg()
		elif method=='nnl1reg':  # non-negative fit, branch terms L1 regularized, avidity terms L2 regularized
			self.params = self.fit_nnl1reg()

		print "method",method, "regularized by", self.lam_HI, "rms deviation=",np.sqrt(self.fit_func())
		# for each set of branches with HI constraints, pick the branch with most aa mutations
		# and assign the dHI to that one, record the number of constraints
		if self.map_to_tree:
			for node in self.tree.postorder_node_iter():
				node.dHI=0
			for HI_split, branches in self.HI_split_to_branch.iteritems():
				likely_branch = branches[np.argmax([len(b.mutations) for b in branches])]
				likely_branch.dHI = self.params[HI_split]
				likely_branch.constraints = self.tree_graph[:,HI_split].sum()

			# integrate the HI change dHI into a cumulative antigentic evolution score cHI
			for node in self.tree.preorder_node_iter():
				if node.parent_node is not None:
					node.cHI = node.parent_node.cHI + node.dHI
				else:
					node.cHI=0
		else:
			self.mutation_effects={}
			for mi, mut in enumerate(self.relevant_muts):
				self.mutation_effects[mut] = self.params[mi]

		self.serum_potency['tree' if self.map_to_tree else 'mutation'] =\
					{serum:self.params[self.genetic_params+ii]
					  for ii, serum in enumerate(self.sera)}
		self.virus_effect['tree' if self.map_to_tree else 'mutation'] = \
					{strain:self.params[self.genetic_params+len(self.sera)+ii]
				  for ii, strain in enumerate(self.HI_strains)}

	def generate_validation_figures(self, method = 'nnl1reg'):
		import matplotlib.pyplot as plt

		lam_pot = self.lam_pot
		lam_avi = self.lam_avi
		lam_HI =  self.lam_HI
		# summary figures using previously determined models
		for map_to_tree, model in [(True, 'tree'), (False,'mutation')]:
			try:
				plt.figure(figsize=(1.3*figheight,figheight))
				ax = plt.subplot(121)
				plt.hist(self.virus_effect[model].values(), bins=np.linspace(-2,2,21), normed=True)
				plt.xlabel('avidity', fontsize=fs)
				plt.text(0.05, 0.93,  ('tree model' if model=='tree' else 'mutation model'),
				         weight='bold', fontsize=fs, transform=plt.gca().transAxes)
				ax.set_xticks([-2,-1,0,1,2])
				ax.tick_params(axis='both', labelsize=fs)
				ax = plt.subplot(122)
				plt.hist(self.serum_potency[model].values(), bins=10, normed=True)
				plt.xlabel('potency', fontsize=fs)
				plt.tight_layout()
				ax.set_xticks([-4,-2,0,2,4])
				ax.tick_params(axis='both', labelsize=fs)
				for fmt in fmts: plt.savefig(self.htmlpath()+'HI_effects_'+model+fmt)
			except:
				print "can't plot effect distributions"


		for map_to_tree, model in [(True, 'tree'), (False,'mutation')]:
			self.map_HI(training_fraction=0.9, method=method,lam_HI=lam_HI, lam_avi=lam_avi,
						lam_pot = lam_pot, force_redo=True, map_to_tree=map_to_tree, subset_strains=True)

			self.validate(plot=True)
			for fmt in fmts: plt.savefig(self.htmlpath()+'HI_prediction_virus_'+model+fmt)

			self.map_HI(training_fraction=0.9, method=method,lam_HI=lam_HI, lam_avi=lam_avi,
						lam_pot = lam_pot, force_redo=True, map_to_tree=map_to_tree)
			self.validate(plot=True)
			for fmt in fmts: plt.savefig(self.htmlpath()+'HI_prediction_'+model+fmt)

		self.save_trunk_cHI()

	def validate(self, plot=False, cutoff=0.0, validation_set = None, incl_ref_strains='yes'):
		if validation_set is None:
			validation_set=self.test_HI
		from scipy.stats import linregress, pearsonr
		self.validation = {}
		for key, val in validation_set.iteritems():
			if self.map_to_tree:
				pred_HI = self.predict_HI_tree(key[0], key[1], cutoff=cutoff)
			else:
				pred_HI = self.predict_HI_mutations(key[0], key[1], cutoff=cutoff)

			if pred_HI is not None:
				if any([incl_ref_strains=='yes',
						incl_ref_strains=='no' and (key[0].upper() not in self.ref_strains),
						incl_ref_strains=='only' and (key[0].upper() in self.ref_strains)]):
					self.validation[key] = (val, pred_HI)

		a = np.array(self.validation.values())
		print "number of prediction-measurement pairs",a.shape
		self.abs_error = np.mean(np.abs(a[:,0]-a[:,1]))
		self.rms_error = np.sqrt(np.mean((a[:,0]-a[:,1])**2))
		self.slope, self.intercept, tmpa, tmpb, tmpc = linregress(a[:,0], a[:,1])
		print "error (abs/rms): ",self.abs_error, self.rms_error
		print "slope, intercept:", self.slope, self.intercept
		print "pearson correlation:", pearsonr(a[:,0], a[:,1])

		if plot:
			import matplotlib.pyplot as plt
			import seaborn as sns
			sns.set_style('darkgrid')
			plt.figure(figsize=(figheight,figheight))
			ax = plt.subplot(111)
			plt.plot([-1,6], [-1,6], 'k')
			plt.scatter(a[:,0], a[:,1])
			plt.ylabel("predicted log2 distance", fontsize = fs)
			plt.xlabel("measured log2 distance" , fontsize = fs)
			ax.tick_params(axis='both', labelsize=fs)
			plt.ylim([-3,8])
			plt.xlim([-3,7])
			plt.text(-2.5,7.3, ('tree model' if self.map_to_tree else 'mutation model'), weight='bold', fontsize=fs)
			plt.text(-2.5,6,'regularization:\nprediction error:', fontsize = fs-2)
			plt.text(1.2,6, str(self.lam_HI)+'/'+str(self.lam_pot)+'/'+str(self.lam_avi)+' (HI/pot/avi)'
			         +'\n'+str(round(self.abs_error, 2))\
					 +'/'+str(round(self.rms_error, 2))+' (abs/rms)', fontsize = fs-2)
			plt.tight_layout()
		return a.shape[0]


	def add_titers(self):
		for ref in self.ref_strains:
			self.node_lookup[ref].HI_titers= defaultdict(dict)
			self.node_lookup[ref].HI_titers_raw= defaultdict(dict)
			self.node_lookup[ref].potency_mut={}
			self.node_lookup[ref].potency_tree={}
			self.node_lookup[ref].autologous_titers = {}
		for ref in self.sera:
			self.node_lookup[ref[0]].autologous_titers[ref[1]] = np.median(self.autologous_titers[ref]['val'])
			if 'mutation' in self.virus_effect:
				self.node_lookup[ref[0]].potency_mut[ref[1]] = self.serum_potency['mutation'][ref]
			if 'tree' in self.virus_effect:
				self.node_lookup[ref[0]].potency_tree[ref[1]] = self.serum_potency['tree'][ref]
		for (test, ref), val in self.HI_normalized.iteritems():
			try:
				self.node_lookup[ref[0]].HI_titers[self.node_lookup[test].clade][ref[1]] = val
			except:
				print("Can't assign",test, ref)
		for (test, ref), val in self.HI_raw.iteritems():
			try:
				self.node_lookup[ref[0]].HI_titers_raw[self.node_lookup[test].clade][ref[1]] = val
			except:
				print("Can't assign",test, ref)
		for test in self.HI_strains:
			if 'mutation' in self.virus_effect:
				self.node_lookup[test].avidity_tree = self.virus_effect['mutation'][test]
			if 'tree' in self.virus_effect:
				self.node_lookup[test].avidity_mut = self.virus_effect['tree'][test]
		for ref in self.ref_strains:
			self.node_lookup[ref].mean_HI_titers = {key:np.mean(titers.values()) for key, titers in
			 									self.node_lookup[ref].HI_titers.iteritems()}
			self.node_lookup[ref].mean_potency_tree = np.mean(self.node_lookup[ref].potency_tree.values())
			self.node_lookup[ref].mean_potency_mut = np.mean(self.node_lookup[ref].potency_mut.values())

	def check_sources(self):
		self.source_HIs = defaultdict(dict)
		with myopen(self.HI_fname, 'r') as infile:
			for line in infile:
				test, ref_virus, serum, src_id, val_str = line.strip().split()
				try:
					val = float(val_str)
					if not np.isnan(val):
						self.source_HIs[src_id][test, (ref_virus, serum)] = self.normalize((ref_virus, serum), float(val))
				except:
					print test, ref_virus, serum, src_id, float(val)

		self.source_validation = {}
		for src_id in self.source_HIs:
			print '\n############### \n',src_id,'\n############### \n'
			print "number of measurements:",len(self.source_HIs[src_id])
			try:
				n_checks = self.validate(validation_set=self.source_HIs[src_id], incl_ref_strains='no')
				self.source_validation[src_id] = [self.abs_error, self.rms_error, self.slope, self.intercept, n_checks]
			except:
				print "skipped due to too few measurements"

	def save_trunk_cHI(self):
		cHI_trunk = []
		trunk_muts = []
		co=0.5
		tmp_pivots = self.tree.seed_node.pivots
		for node in self.tree.postorder_internal_node_iter():
			node.num_date = np.min([c.num_date for c in node.child_nodes()])
			if node.trunk:
				tmp_freq = node.freq['global']
				if tmp_freq[0]>0.5:
					continue
				ii = np.argmax(tmp_freq>co)
				slope = (tmp_freq[ii]-tmp_freq[ii-1])/(tmp_pivots[ii]-tmp_pivots[ii-1])
				tp = tmp_pivots[ii-1] + (co-tmp_freq[ii-1])/slope

				tmp_muts = [x for x in node.aa_muts.items() if x[1]!='']
				for muts in tmp_muts:
					gene = muts[0]
					for pos in muts[1].split(','):
						tmp_mut = (gene, pos)
						if tmp_mut in self.mutation_effects:
							trunk_muts.append([tp,self.mutation_effects[tmp_mut]])
				cHI_trunk.append([tp, node.cHI] )

		cHI_trunk = np.array(cHI_trunk)[::-1]
		trunk_muts = np.array(trunk_muts)[::-1]
		np.savetxt('data/'+self.prefix+self.resolution+'_cHI.txt', cHI_trunk)
		np.savetxt('data/'+self.prefix+self.resolution+'_trunk_muts.txt', trunk_muts)

		from random import sample
		n_leafs = 1
		leaf_sample = []
		for y in range(int(self.time_interval[0]), int(self.time_interval[1])):
			tmp_leafs = [leaf for leaf in self.tree.leaf_iter() if leaf.num_date>y and leaf.num_date<y+1]
			leaf_sample.extend(sample(tmp_leafs, min(n_leafs, len(tmp_leafs))))

		cHI_trunk = []
		for node in leaf_sample:
			tmp = []
			while node.parent_node is not None:
				tmp.append([node.num_date, node.cHI])
				node = node.parent_node
			cHI_trunk.append(np.array(tmp)[::-1])
		import cPickle as pickle
		with open('data/'+self.prefix+self.resolution+'_cHI_path.pkl', 'w') as ofile:
			pickle.dump(cHI_trunk, ofile)

	def predict_HI_tree(self, virus, serum, cutoff=0.0):
		path = self.get_path_no_terminals(virus,serum[0])
		if path is not None:
			return self.serum_potency['tree'][serum] \
					+ self.virus_effect['tree'][virus] \
					+ np.sum([b.dHI for b in path if b.dHI>cutoff])
		else:
			return None

	def predict_HI_mutations(self, virus, serum, cutoff=0.0):
		muts= self.get_mutations(serum[0], virus)
		if muts is not None:
			return self.serum_potency['mutation'][serum] \
					+ self.virus_effect['mutation'][virus] \
					+ np.sum([self.mutation_effects[mut] for mut in muts
					if (mut in self.mutation_effects and self.mutation_effects[mut]>cutoff)])
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
		ref_sera = [[HI_fix_name(e1+'/'+e2), e3.replace(' ','')] for e1,e2,e3 in zip(row1, row2, row3)[4:]]
		for ri in xrange(len(ref_sera)):
			abbr = ref_sera[ri][0].split('/')[1].rstrip('01234566789')
			if abbr in name_abbrev:
				ref_sera[ri][0] = HI_fix_name(ref_sera[ri][0].replace(abbr, name_abbrev[abbr]))
			else:
				ref_sera[ri][0] = HI_fix_name(ref_sera[ri][0])
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
				ref_strains.append(HI_fix_name(row[0].strip()))
				ref_matrix.append([src_id,'ref']+map(strip, row[1:4])+map(titer_to_number, row[4:]))

		test_strains = []
		test_matrix = []
		for row in csv_reader: # load test viruses until it is no longer an A/ flu  name
			if not (row[0].startswith('A/') or row[0].startswith('B/')):
				break
			else:
				test_strains.append(HI_fix_name(row[0].strip()))
				test_matrix.append([src_id,'test']+map(strip,row[1:4])+map(titer_to_number, row[4:]))

		print len(ref_sera), ref_sera
		print len(ref_strains), len(test_strains)
		HI_table  = pd.DataFrame(ref_matrix+test_matrix, index = ref_strains+test_strains, columns= fields)

		return HI_table


def read_tables(flutype = 'H3N2'):
	import glob
	from itertools import product
	flist = glob.glob('../../HI_titers/'+flutype+'_tables/NIMR*csv')
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
	measurements = []
	sera = set()
	strains = set()
	if os.path.isfile(trevor_table):
		with myopen(trevor_table) as infile:
			table_reader = csv.reader(infile, delimiter="\t")
			header = table_reader.next()
			for row in table_reader:
				val = titer_to_number(row[6])
				if not np.isnan(val):
					strains.add(HI_fix_name(row[1]))
					serum = (HI_fix_name(row[4]), row[3])
					src_id = row[-1]
					sera.add(serum)
					measurements.append([HI_fix_name(row[1]), serum, src_id, val])
	else:
		print trevor_table, "not found"
	print "trevor total:", len(measurements), "measurements"
	return strains, sera, pd.DataFrame(measurements)


def table_to_flat(HI_table):
	flat_measurements = list()
	for ref_serum in HI_table.columns[5:]:
		sub_set_vals = HI_table[ref_serum][~np.isnan(HI_table[ref_serum])]
		sub_set_source = HI_table['source'][~np.isnan(HI_table[ref_serum])]
		for virus, val, src_id in izip(sub_set_vals.index, sub_set_vals, sub_set_source):
			flat_measurements.append([virus, ref_serum, src_id, val])
	print "NIMR total:", len(flat_measurements), "measurements"
	return pd.DataFrame(flat_measurements)

def get_all_titers_flat(flutype='H3N2'):
	HI_titers = read_tables(flutype)
	HI_titers_flat = table_to_flat(HI_titers)
	HI_trevor = read_trevor_table(flutype)[2]
	HI_titers_flat = pd.concat((HI_titers_flat, HI_trevor))
	return HI_titers_flat


def write_strains_with_HI_and_sequence(flutype='H3N2'):
	HI_titers = read_tables(flutype)
	HI_trevor = read_trevor_table(flutype)
	HI_strains = set(HI_titers.index)
	HI_strains.update(HI_trevor[0])
	from Bio import SeqIO
	good_strains = set()
	with myopen("data/"+flutype+"_strains_with_HI.fasta", 'w') as outfile, \
		 myopen("source-data/"+flutype+"_HI_strains.txt", 'w') as HI_strain_outfile, \
		 myopen("data/"+flutype+"_gisaid_epiflu_sequence.fasta.gz", 'r') as infile:
		for seq_rec in SeqIO.parse(infile, 'fasta'):
			tmp_name = seq_rec.description.split('|')[0].strip()
			reduced_name = HI_fix_name(tmp_name)
			if reduced_name in HI_strains and (reduced_name not in good_strains):
				SeqIO.write(seq_rec, outfile,'fasta')
				good_strains.add(reduced_name)
				HI_strain_outfile.write(tmp_name+'\n')
				if fix_name(tmp_name)!=tmp_name:
					HI_strain_outfile.write(fix_name(tmp_name)+'\n')
				#print seq_rec.name


def write_flat_HI_titers(flutype = 'H3N2', fname = None):
	measurements = get_all_titers_flat(flutype)
	with myopen('source-data/'+flutype+'_HI_strains.txt') as infile:
		strains = [HI_fix_name(line.strip()).upper() for line in infile]
	if fname is None:
		fname = 'source-data/'+flutype+'_HI_titers.txt'
	written = 0
	skipped = 0
	with myopen(fname, 'w') as outfile:
		for ii, rec in measurements.iterrows():
			test, ref, src_id, val = rec
			if HI_fix_name(test).upper() in strains and HI_fix_name(rec[1][0]).upper() in strains:
				outfile.write('\t'.join(map(str, [test, ref[0], ref[1], src_id, val]))+'\n')
				written+=1
			else:
				skipped+=1
	print "written",written,"records"
	print "skipped",skipped,"records"

def main(tree, HI_fname='source-data/HI_titers.txt', training_fraction = 1.0, reg=5):
	print "--- Fitting HI titers at " + time.strftime("%H:%M:%S") + " ---"
	measurements, strains, sources = read_HI_titers(HI_fname)
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
