"""
Tests for the `fitness_model` module.
"""
import Bio.SeqIO
import dendropy
import pytest

try:
	from StringIO import StringIO
except ImportError:
	from io import StringIO

#
# Fixtures
#

@pytest.fixture
def tree():
	"""Returns a tree with three sequences: a root and two direct descendents with
	one modification each.
	"""
	# Build simple tree.
	tree = dendropy.Tree(stream=StringIO("(A,B);"), schema="newick")

	# Build sequences for tree nodes. One leaf has a Koel and epitope site
	# mutation. The other leaf has a signal peptide mutation.
	root = sequence()
	leaf_a = modify_sequence_at_site(root, 145 + 16 - 1)
	leaf_b = modify_sequence_at_site(root, 14)

	# Assign sequences to nodes.
	sequences = [root, leaf_a, leaf_b]
	index = 0
	for node in tree.preorder_node_iter():
		node.aa = sequences[index]
		index += 1

	return tree

@pytest.fixture
def fitness_model(tree):
	from src.fitness_model import fitness_model
	return fitness_model(
		pivots_per_year=12,
		time_interval=(2012.0, 2015.0),
		tree=tree
	)

@pytest.fixture
def sequence():
	"""Returns an amino acid sequence for an ancestral H3N2 virus (Hong Kong 1968).
	"""
	with open("tests/data/AAK51718.fasta", "r") as handle:
		record = list(Bio.SeqIO.parse(handle, "fasta"))[0]

	aa = str(record.seq)
	return aa

#
# Utility functions
#

def modify_sequence_at_site(sequence, site):
	"""Returns the given sequence with a modified base at the given site.
	"""
	other_sequence_list = list(sequence)
	other_sequence_list[site] = "Z"
	return "".join(other_sequence_list)

#
# Tests
#

class TestFitnessModel(object):
	def test_prep_nodes(self, fitness_model):
		assert not hasattr(fitness_model, "nodes")
		fitness_model.prep_nodes()
		assert hasattr(fitness_model, "nodes")
		assert hasattr(fitness_model, "rootnode")
		assert hasattr(fitness_model.rootnode, "pivots")

	def test_calc_node_frequencies(self, fitness_model):
		fitness_model.prep_nodes()
		assert not hasattr(fitness_model, "freq_arrays")
		fitness_model.calc_node_frequencies()
		assert hasattr(fitness_model, "freq_arrays")
		assert len(fitness_model.freq_arrays) > 0
