from itertools import izip
import numpy as np

def hamming_distance(seq1, seq2):
	aseq1, aseq2 = np.array(seq1), np.array(seq2)
	non_gap = (aseq1!='-')*(aseq2!='-')
	return np.mean(aseq1[non_gap]!=aseq2[non_gap])

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
	'''
	Receptor binding site mutations from Koel et al. 2014
	These are (145, 155, 156, 158, 159, 189, 193) in canonical HA numbering
	need to subtract one since python arrays start at 0
	'''
	sites = [144, 154, 155, 157, 158, 188, 192]
	return ''.join([aa[pos] for pos in sites])

def get_HA1(aa):
	'''
	return the part of the peptide corresponding to HA1, starts is 329 aa long
	'''
	return aa[:329]

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
	aa = translate(nuc[48:])
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