from itertools import izip

def partition_string(string, length):
	return list(string[0+i:length+i] for i in range(0, len(string), length))

def translate(nuc):
	"""Translate nucleotide sequence to amino acid"""

	to_aa = {}
	to_aa["TTT"] = "F"
	to_aa["TTC"] = "F"
	to_aa["TTA"] = "L"
	to_aa["TTG"] = "L"

	to_aa["CTT"] = "L"
	to_aa["CTC"] = "L"
	to_aa["CTA"] = "L"
	to_aa["CTG"] = "L"

	to_aa["ATT"] = "I"
	to_aa["ATC"] = "I"
	to_aa["ATA"] = "I"
	to_aa["ATG"] = "M"

	to_aa["GTT"] = "V"
	to_aa["GTC"] = "V"
	to_aa["GTA"] = "V"
	to_aa["GTG"] = "V"

	to_aa["TCT"] = "S"
	to_aa["TCC"] = "S"
	to_aa["TCA"] = "S"
	to_aa["TCG"] = "S"

	to_aa["CCT"] = "P"
	to_aa["CCC"] = "P"
	to_aa["CCA"] = "P"
	to_aa["CCG"] = "P"

	to_aa["ACT"] = "T"
	to_aa["ACC"] = "T"
	to_aa["ACA"] = "T"
	to_aa["ACG"] = "T"

	to_aa["GCT"] = "A"
	to_aa["GCC"] = "A"
	to_aa["GCA"] = "A"
	to_aa["GCG"] = "A"

	to_aa["TAT"] = "Y"
	to_aa["TAC"] = "Y"
	to_aa["TAA"] = "X"
	to_aa["TAG"] = "X"

	to_aa["CAT"] = "H"
	to_aa["CAC"] = "H"
	to_aa["CAA"] = "Q"
	to_aa["CAG"] = "Q"

	to_aa["AAT"] = "N"
	to_aa["AAC"] = "N"
	to_aa["AAA"] = "K"
	to_aa["AAG"] = "K"

	to_aa["GAT"] = "D"
	to_aa["GAC"] = "D"
	to_aa["GAA"] = "E"
	to_aa["GAG"] = "E"

	to_aa["TGT"] = "C"
	to_aa["TGC"] = "C"
	to_aa["TGA"] = "X"
	to_aa["TGG"] = "W"

	to_aa["CGT"] = "R"
	to_aa["CGC"] = "R"
	to_aa["CGA"] = "R"
	to_aa["CGG"] = "R"

	to_aa["AGT"] = "S"
	to_aa["AGC"] = "S"
	to_aa["AGA"] = "R"
	to_aa["AGG"] = "R"

	to_aa["GGT"] = "G"
	to_aa["GGC"] = "G"
	to_aa["GGA"] = "G"
	to_aa["GGG"] = "G"

	partitions = partition_string(nuc, 3)
	aa = ''.join([to_aa[codon] for codon in partitions])
	return aa

def epitope_sites(nuc):
	epitope_mask = "00000000000000000000000000000000000000000000000000000000000011111011011001010011000100000001001011110011100110101000001100000100000001000110101011111101011010111110001010011111000101011011111111010010001111101110111001010001110011111111000000111110000000101010101110000000000011100100000001011011100000000000001001011000110111111000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000"
	aa = translate(nuc)
	zipped = zip(aa, epitope_mask)
	epi = ""
	for n in zipped:
		if n[1] == '1':
			epi += n[0]
	return epi

def nonepitope_sites(nuc):
	nonepitope_mask = "11111111111111111111111111111111111111111111111111111111111100000100100110101100111011111110110100001100011001010111110011111011111110111001010100000010100101000001110101100000111010100100000000101101110000010001000110101110001100000000111111000001111111010101010001111111111100011011111110100100011111111111110110100111001000000111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111"
	aa = translate(nuc)
	zipped = zip(aa, nonepitope_mask)
	epi = ""
	for n in zipped:
		if n[1] == '1':
			epi += n[0]
	return epi

def receptor_binding_sites(nuc):
	"""Receptor binding site mutations from Koel et al. 2014"""
	"""These are (145, 155, 156, 158, 159, 189, 193) in canonical HA numbering"""
	"""When counting from ATG/M, need to offset by 16, giving (161, 171, 172, 174, 175, 205, 209)"""
	"""When indexing from 0, these are (160, 170, 171, 173, 174, 204, 208)"""
	sites = (160, 170, 171, 173, 174, 204, 208);
	aa = translate(nuc)
	rb = ""
	for site in sites:
		rb += aa[site]
	return rb

def epitope_distance(nucA, nucB):
	"""Return distance of sequences nucA and nucB by comparing epitope sites"""
	epA = epitope_sites(nucA)
	epB = epitope_sites(nucB)
	distance = sum(a != b for a, b in izip(epA, epB))
	return distance

def nonepitope_distance(nucA, nucB):
	"""Return distance of sequences nucA and nucB by comparing non-epitope sites"""
	neA = nonepitope_sites(nucA)
	neB = nonepitope_sites(nucB)
	distance = sum(a != b for a, b in izip(neA, neB))
	return distance

def receptor_binding_distance(nucA, nucB):
	"""Return distance of sequences nucA and nucB by comparing receptor binding sites"""
	neA = receptor_binding_sites(nucA)
	neB = receptor_binding_sites(nucB)
	distance = sum(a != b for a, b in izip(neA, neB))
	return distance

def main():
	"""Testing with Hong Kong/68"""
	nuc = "ATGAAGACCATCATTGCTTTGAGCTACATTTTCTGTCTGGCTCTCGGCCAAGACCTTCCAGGAAATGACAACAGCACAGCAACGCTGTGCCTGGGACATCATGCGGTGCCAAACGGAACACTAGTGAAAACAATCACAGATGATCAGATTGAAGTGACTAATGCTACTGAGCTAGTTCAGAGCTCCTCAACGGGGAAAATATGCAACAATCCTCATCGAATCCTTGATGGAATAGACTGCACACTGATAGATGCTCTATTGGGGGACCCTCATTGTGATGTTTTTCAAAATGAGACATGGGACCTTTTCGTTGAACGCAGCAAAGCTTTCAGCAACTGTTACCCTTATGATGTGCCAGATTATGCCTCCCTTAGGTCACTAGTTGCCTCGTCAGGCACTCTGGAGTTTATCACTGAGGGTTTCACTTGGACTGGGGTCACTCAGAATGGGGGAAGCAATGCTTGCAAAAGGGGACCTGGTAGCGGTTTTTTCAGTAGACTGAACTGGTTGACCAAATCAGGAAGCACATATCCAGTGCTGAACGTGACTATGCCAAACAATGACAATTTTGACAAACTATACATTTGGGGGGTTCACCACCCGAGCACGAACCAAGAACAAACCAGCCTGTATGTTCAAGCATCAGGGAGAGTCACAGTCTCTACCAGAAGAAGCCAGCAAACTATAATCCCGAATATCTGGTCCAGACCCTGGGTAAGGGGTCTGTCTAGTAGAATAAGCATCTATTGGACAATAGTTAAGCCGGGAGACGTACTGGTAATTAATAGTAATGGGAACCTAATCGCTCCTCGGGGTTATTTCAAAATGCGCACTGGGAAAAGCTCAATAATGAGGTCAGATGCACCTATTGATACCTGTATTTCTGAATGCATCACTCCAAATGGAAGCATTCCCAATGACAAGCCCTTTCAAAACGTAAACAAGATCACATATGGAGCATGCCCCAAGTATGTTAAGCAAAACACC"
	aa = translate(nuc)
	ep = epitope_sites(nuc)
	ne = nonepitope_sites(nuc)
	rb = receptor_binding_sites(nuc)
	print "nuc: " + nuc
	print "aa: " + aa
	print "ep: " + ep
	print "ne: " + ne
	print "rb: " + rb

if __name__ == "__main__":
	main()