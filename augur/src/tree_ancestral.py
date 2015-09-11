'''Reconstruct ancestral sequences in a tree

The module consists of a single class, ancestral_sequences.

Usage:

a = ancestral_sequences(tree, aln, ...)
a.calc_ancestral_sequences()

or, slightly less Pythonic,

(ancestral_sequences(tree, aln, ...)
	.calc_ancestral_sequences())

'''


import numpy as np
from Bio import Seq
import copy, time

class ancestral_sequences:
	'''
	class that generates a biopython tree dressed with ancestral sequences
	and the marginal probabilities of different states in the tree
	NOTE: THIS TREATS ANYTHING OTHER THAN ACGT AS MISSING DATA. NO GAPHANDLING
	should be extended to handle ambiguous symbols correctly.
	'''

	def __init__(self, tree, aln, alphabet = 'ACGT', sub_matrix = None,
				 eps_branch_length = 1e-7, copy_tree = False,
				 attrname='seq', seqtype='Seq'):
		'''
		arguments:
		tree  -- a biopython tree
		aln   -- a biopython alignment with the same names as the terminal nodes

		keyword arguments:
		alphabet   -- allows character states
		sub_matrix -- substitution matrix. defaults to flat matrix
		eps_branch_length -- minimal branch length to prevent division by zero exception
		copy_tree -- if true, a new tree object is constructed
		'''
		if copy_tree:
			self.T = copy.deepcopy(tree)
		else:
			self.T = tree
		self.alphabet = np.array(list(alphabet))
		self.nstates = self.alphabet.shape[0]
		self.seq_len = aln.get_alignment_length()
		self.pseudo_branch_length = eps_branch_length
		self.attrname = attrname
		self.seqtype = seqtype

		# construct substitution matrix if not provided
		if sub_matrix:
			self.sub_matrix = np.array(sub_matrix)
		else:
			# every mutation equally likely. subtract diagonal, normalize to rate 1.
			# this matrix is symmetric
			self.sub_matrix = (np.ones((self.nstates, self.nstates))-
							   self.nstates*np.eye(self.nstates))*1.0/(self.nstates-1)
			self.calc_eigendecomp()

		names_to_seqs = {seq.id:seq  for seq in aln}
		for leaf in self.T.leaf_iter():
			if leaf.taxon.label in names_to_seqs:
				leaf.prob = self.get_state_array()

				seqtmp = names_to_seqs[leaf.taxon.label].seq

				if self.seqtype != 'Seq':
					seqtmp = ''.join(seqtmp).upper()
				setattr(leaf, self.attrname, seqtmp)

				# convert to a numpy array for convenient slicing
				tmp_seq_array = np.fromstring(''.join(getattr(leaf, self.attrname)), 'S1')

				# code the sequences as a 0/1 probability matrix  (redundant but convenient)
				for ni in xrange(self.nstates):
					leaf.prob[:,ni] = tmp_seq_array == self.alphabet[ni]
				missing_prob = (1.0-leaf.prob.sum(axis=1))/leaf.prob.shape[1]
				for ni in xrange(self.nstates):
					leaf.prob[:,ni]+=missing_prob
			else:
				print('ancestral sequences: Leaf '+leaf.taxon.label+' has no sequence')

		if self.seqtype == 'Seq':
			self.biopython_alphabet = leaf.seq.alphabet

		# dress each internal node with a probability vector for the ancestral state
		for node in self.T.postorder_internal_node_iter():
			node.prob = self.get_state_array()
			if node.edge_length is None:
				node.edge_length = 0.0
		# there is no parental information to the root (the root node is artificial)
		# hence init the message with ones
		self.T.seed_node.up_message = self.get_state_array()


	def get_state_array(self):
		'''
		provide a unified function that returns an all one array
		for messages and probabilities of the right size
		'''
		return np.ones((self.seq_len,self.nstates))

	def calc_eigendecomp(self):
		'''
		calculates eigenvalues and eigenvectors.
		note that this assumes that the substitution matrix is symmetric
		'''
		self.evals, self.evecs = np.linalg.eigh(self.sub_matrix)

	def calc_state_probabilites(self,P, t):
		'''
		input: initial state, time interval
		return the solution of the character evolution equaiton
		'''
		transition_matrix = np.dot(self.evecs, np.dot(np.diag(
					np.exp((t+self.pseudo_branch_length)*self.evals)), self.evecs.T))
		return np.dot(P, transition_matrix.T)

	def normalize(self, clade):
		'''
		normalize the distribution of states at each position such that the sum equals one
		'''
		clade.prob/=np.repeat(np.array([np.sum(clade.prob, axis=1)]).T, self.nstates, axis=1)

		if np.isnan(np.sum(clade.prob[:])):
			print "encountered nan in ancestral inference in clade ", clade.taxon.label
			print np.isnan(clade.prob).nonzero()

	def log_normalize(self, clade):
		'''
		convert the unnormalized and logarithmic probabilites array to linear and then normaliz
		normalize the distribution of states at each position such that the sum equals one
		'''
		# substract the maximum value in each column
		clade.prob -= np.repeat(np.array([np.max(clade.prob, axis=1)]).T, self.nstates, axis=1)
		# exponentiate
		clade.prob = np.exp(clade.prob)
		#normalize
		self.normalize(clade)

	def calc_down_messages(self,clade):
		'''
		recursively calculates the messages passed on the parents of each node
		input: clade whose down_message is to be calculated
		'''
		if clade.is_leaf():
			# if clade is terminal, the sequence is fix and we can emit the state probabilities
			clade.down_message = self.calc_state_probabilites(clade.prob, clade.edge_length)
			#print "down clade", clade.taxon.label, 'min:', np.min(clade.down_message)
			clade.down_message[clade.down_message<1e-30] = 1e-30
		else:
			#otherwise, multiply all down messages from children, normalize and pass down
			clade.prob[:]=0
			for child in clade.child_nodes():
				self.calc_down_messages(child)
				clade.prob+=np.log(child.down_message)

			self.log_normalize(clade)
			clade.down_message = self.calc_state_probabilites(clade.prob, clade.edge_length)
			#print "down clade", clade.taxon.label, 'min:', np.min(clade.down_message)
			clade.down_message[clade.down_message<1e-30] = 1e-30

	def calc_up_messages(self,clade):
		'''
		calculate the messages that are passed on to the children
		input calde for which these are to calculated
		'''
		if clade.is_leaf():
			#nothing to be done for terminal nodes
			return
		else:
			#else, loop over children and calculate the message for each of the children
			for child in clade.child_nodes():
				# initialize with the message comming from the parent
				clade.prob[:]=np.log(clade.up_message)
				for child2  in clade.child_nodes():
					if child2 != child:
						#multiply with the down message from each of the children, but skip child 1
						clade.prob+=np.log(child2.down_message)

				# normalize, adjust for modifications along the branch, and save.
				self.log_normalize(clade)
				child.up_message = self.calc_state_probabilites(clade.prob, child.edge_length)
				#print "up clade", clade.taxon.label, 'min:', np.min(child.up_message)
				child.up_message[child.up_message<1e-30] = 1e-30
			# do recursively for all children
			for child in clade.child_nodes():
				self.calc_up_messages(child)

	def calc_marginal_probabilities(self,clade):
		'''
		calculate the marginal probabilities by multiplying all incoming messages
		'''
		if clade.is_internal():
			clade.prob[:]=np.log(clade.up_message)
			for child in clade.child_nodes():
				clade.prob+=np.log(child.down_message)
			# normalize and continue for all children
			self.log_normalize(clade)
			#print clade.taxon.label, np.max(1.0-np.max(clade.prob, axis=1))
			for child in clade.child_nodes():
				self.calc_marginal_probabilities(child)
		else:
			tmp_seq_array = np.fromstring(''.join(getattr(clade, self.attrname)), 'S1')
			# code the sequences as a 0/1 probability matrix  (redundant but convenient)
			for ni in xrange(self.nstates):
				clade.prob[:,ni] = tmp_seq_array == self.alphabet[ni]
			missing_seq = clade.prob.sum(axis=1)==0
			clade.prob[missing_seq] = clade.up_message[missing_seq]


	def calc_most_likely_sequences(self, clade):
		'''
		recursively calculate the most likely sequences for each node
		'''
		seq = "".join(self.alphabet[np.argmax(clade.prob, axis=1)])
		if self.seqtype == 'Seq':
			seq = Seq.Seq(seq, alphabet = self.biopython_alphabet)

		setattr(clade, self.attrname, seq)

		# repeat for all children
		for child in clade.child_nodes():
			if child.is_internal():
				self.calc_most_likely_sequences(child)


	def calc_ancestral_sequences(self):
		'''
		given the initialized instance, calculates the most likely ancestral sequences
		and the marginal probabilities for each position at each internal node.
		'''
		print "--- Forward pass at " + time.strftime("%H:%M:%S") + " ---"
		self.calc_down_messages(self.T.seed_node)
		print "--- Backward pass at " + time.strftime("%H:%M:%S") + " ---"
		self.calc_up_messages(self.T.seed_node)
		print "--- Calculating marginals at " + time.strftime("%H:%M:%S") + " ---"
		self.calc_marginal_probabilities(self.T.seed_node)
		print "--- Most likely nucleotides at " + time.strftime("%H:%M:%S") + " ---"
		self.calc_most_likely_sequences(self.T.seed_node)
		self.cleanup_tree()


	def cleanup_tree(self, attrnames=['prob', 'down_message', 'up_message']):
		'''Clean up pollution attributes of leaves'''
		nodes = self.T.postorder_node_iter()
		for leaf in nodes:
			for attrname in attrnames:
				if hasattr(leaf, attrname):
					delattr(leaf, attrname)

def main(tree, viruses):
	from seq_util import json_to_Bio_alignment
	from tree_util import json_to_dendropy
	print "--- Ancestral inference at " + time.strftime("%H:%M:%S") + " ---"
	aln = json_to_Bio_alignment(viruses)
	print "--- Set-up ancestral inference at " + time.strftime("%H:%M:%S") + " ---"
	anc_seq = ancestral_sequences(tree, aln, seqtype='str')
	anc_seq.calc_ancestral_sequences()
	anc_seq.cleanup_tree()
	out_fname = "data/tree_ancestral.json"
	return json_to_dendropy(anc_seq.T.seed_node)

if __name__=="__main__":
	tree = main()
