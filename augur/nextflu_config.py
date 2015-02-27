config = {
	'virus':'H3N2',
	'alignment_file':'data/20150222_all_H3N2_HA1.fasta',
	'fasta_fields':{0:'strain', 1:"date", 4:"passage", -1:'accession'},
	'outgroup':'A/Beijing/32/1992',
	'max_global':True,   # sample as evenly as possible from different geographic regions 
	'aggregate_regions': [  ("global", None), ("NA", ["NorthAmerica"]), ("EU", ["Europe"]), 
							("AS", ["China", "SoutheastAsia", "JapanKorea"]), ("OC", ["Oceania"]) ],
	'frequency_stiffness':1.0,
	'time_interval':(2012.0, 2015.1),
	'pivots_per_year':12.0,
	'cds':[48,-1], # define the HA1 start i n 0 numbering
	# define relevant clades in canonical HA1 numbering (+1)
	'clade_designations': { "3c3.a":[(128,'A'), (142,'G'), (159,'S')],
						   "3c3":  [(128,'A'), (142,'G'), (159,'F')],
						   "3c2.a":[(144,'S'), (159,'Y'), (225,'D'), (311,'H'),(489,'N')],
						   "3c2":  [(144,'N'), (159,'F'),(225,'N'), (489,'N')]
							}
}
