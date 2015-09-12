from itertools import izip
import numpy as np

def hamming_distance(seq1, seq2):
	aseq1, aseq2 = np.array(seq1), np.array(seq2)
	non_gap = (aseq1!='-')*(aseq2!='-')
	return np.mean(aseq1[non_gap]!=aseq2[non_gap])

def translate(nuc):
	"""Translate nucleotide sequence to amino acid"""
	from Bio import Seq
	try:
		tmp_aa = Seq.translate(nuc.replace('-','N')) #returns string when argument is a string, Bio.Seq otherwise
	except:
		print("translation failed",nuc)
		tmp_aa = 'X'*len(nuc)//3
	aa_seq = ""
	for i,aa in enumerate(tmp_aa):
		if nuc[i*3:(i+1)*3]=='---':
			aa_seq+='-'
		else:
			aa_seq+=aa
	return aa_seq

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