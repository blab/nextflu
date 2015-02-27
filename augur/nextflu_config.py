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
	'pivots_per_year':6.0
}
