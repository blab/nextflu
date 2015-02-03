from itertools import izip
import numpy as np
epitope_mask = np.fromstring("00000000000000000000000000000000000000000000000000000000000011111011011001010011000100000001001011110011100110101000001100000100000001000110101011111101011010111110001010011111000101011011111111010010001111101110111001010001110011111111000000111110000000101010101110000000000011100100000001011011100000000000001001011000110111111000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000", dtype='S1')

def partition_string(string, length):
	return list(string[0+i:length+i] for i in range(0, len(string), length))

def translate(nuc):
	"""Translate nucleotide sequence to amino acid"""
	from Bio import Seq
	return Seq.translate(nuc) #returns string when argument is a string, Bio.Seq otherwise

def epitope_sites(aa):
	aaa = np.fromstring(aa, 'S1')
	return ''.join(aaa[epitope_mask[:len(aa)]=='1'])

def nonepitope_sites(aa):
	aaa = np.fromstring(aa, 'S1')
	return ''.join(aaa[epitope_mask[:len(aa)]=='0'])

def receptor_binding_sites(aa):
	"""Receptor binding site mutations from Koel et al. 2014"""
	"""These are (145, 155, 156, 158, 159, 189, 193) in canonical HA numbering"""
	"""When counting from ATG/M, need to offset by 16, giving (161, 171, 172, 174, 175, 205, 209)"""
	"""When indexing from 0, these are (160, 170, 171, 173, 174, 204, 208)"""
	sites = [160, 170, 171, 173, 174, 204, 208]
	return ''.join([aa[pos] for pos in sites])


def epitope_distance(aaA, aaB):
	"""Return distance of sequences aaA and aaB by comparing epitope sites"""
	epA = epitope_sites(aaA)
	epB = epitope_sites(aaB)
	distance = sum(a != b for a, b in izip(epA, epB))
	return distance

def nonepitope_distance(aaA, aaB):
	"""Return distance of sequences aaA and aaB by comparing non-epitope sites"""
	neA = nonepitope_sites(aaA)
	neB = nonepitope_sites(aaB)
	distance = sum(a != b for a, b in izip(neA, neB))
	return distance
def receptor_binding_distance(aaA, aaB):
	"""Return distance of sequences aaA and aaB by comparing receptor binding sites"""
	neA = receptor_binding_sites(aaA)
	neB = receptor_binding_sites(aaB)
	distance = sum(a != b for a, b in izip(neA, neB))
	return distance

def json_to_Bio_alignment(seq_json):
	from Bio.Align import MultipleSeqAlignment
	from Bio.SeqRecord import SeqRecord
	from Bio.Seq import Seq
	aln = MultipleSeqAlignment([])
	for seq in seq_json:
		aln.append(SeqRecord(name=seq['strain'], id=seq['strain'], seq=Seq(seq['seq'])))
	return aln

def main():
	"""Testing with Hong Kong/68"""
	nuc = "ATGAAGACCATCATTGCTTTGAGCTACATTTTCTGTCTGGCTCTCGGCCAAGACCTTCCAGGAAATGACAACAGCACAGCAACGCTGTGCCTGGGACATCATGCGGTGCCAAACGGAACACTAGTGAAAACAATCACAGATGATCAGATTGAAGTGACTAATGCTACTGAGCTAGTTCAGAGCTCCTCAACGGGGAAAATATGCAACAATCCTCATCGAATCCTTGATGGAATAGACTGCACACTGATAGATGCTCTATTGGGGGACCCTCATTGTGATGTTTTTCAAAATGAGACATGGGACCTTTTCGTTGAACGCAGCAAAGCTTTCAGCAACTGTTACCCTTATGATGTGCCAGATTATGCCTCCCTTAGGTCACTAGTTGCCTCGTCAGGCACTCTGGAGTTTATCACTGAGGGTTTCACTTGGACTGGGGTCACTCAGAATGGGGGAAGCAATGCTTGCAAAAGGGGACCTGGTAGCGGTTTTTTCAGTAGACTGAACTGGTTGACCAAATCAGGAAGCACATATCCAGTGCTGAACGTGACTATGCCAAACAATGACAATTTTGACAAACTATACATTTGGGGGGTTCACCACCCGAGCACGAACCAAGAACAAACCAGCCTGTATGTTCAAGCATCAGGGAGAGTCACAGTCTCTACCAGAAGAAGCCAGCAAACTATAATCCCGAATATCTGGTCCAGACCCTGGGTAAGGGGTCTGTCTAGTAGAATAAGCATCTATTGGACAATAGTTAAGCCGGGAGACGTACTGGTAATTAATAGTAATGGGAACCTAATCGCTCCTCGGGGTTATTTCAAAATGCGCACTGGGAAAAGCTCAATAATGAGGTCAGATGCACCTATTGATACCTGTATTTCTGAATGCATCACTCCAAATGGAAGCATTCCCAATGACAAGCCCTTTCAAAACGTAAACAAGATCACATATGGAGCATGCCCCAAGTATGTTAAGCAAAACACC"
	aa = translate(nuc)
	ep = epitope_sites(aa)
	ne = nonepitope_sites(aa)
	rb = receptor_binding_sites(aa)
	print "nuc: " + nuc
	print "aa: " + aa
	print "ep: " + ep
	print "ne: " + ne
	print "rb: " + rb

if __name__ == "__main__":
	main()