import time, argparse,re,os, socket
import matplotlib as mpl
if socket.gethostname() not in ['olt', 'rneher-iMac']:
	mpl.use('pdf')
from virus_filter import flu_filter, fix_name
from virus_clean import virus_clean
from tree_refine import tree_refine
from tree_titer import HI_tree
from fitness_model import fitness_model
from process import process, virus_config
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Align import MultipleSeqAlignment
import numpy as np
from itertools import izip

# HA2 AA sites are shifted by +329 relative to HA1
# So HA2:77V is 406V in HA1 numbering

virus_config.update({
	# data source and sequence parsing/cleaning/processing
	'virus':'H3N2',
	'alignment_file':'data/h3n2.fasta',
	'outgroup':'A/Beijing/32/1992',
	'force_include':'data/h3n2_hi_strains.tsv',
	'force_include_all':False,
	'date_spec':'year',
	'max_global':True,   # sample as evenly as possible from different geographic regions
	'cds':[0,None], # define the HA1 start i n 0 numbering
	'n_iqd':5,
	'min_mutation_frequency':0.01,
	# define relevant clades in canonical HA1 numbering (+1)
	# numbering starting at HA1 start, adding sp to obtain numbering from methionine
	'clade_designations': { "3c3.a":[('HA1',128,'A'), ('HA1',142,'G'), ('HA1',159,'S')],
						   "3c3":   [('HA1',128,'A'), ('HA1',142,'G'), ('HA1',159,'F')],
						   "3c2.a": [('HA1',144,'S'), ('HA1',159,'Y'), ('HA1',225,'D'), ('HA1',311,'H'), ('HA2',160,'N')],
						   "171K": [('HA1',144,'S'), ('HA1',159,'Y'), ('HA1',171,'K'), ('HA1',225,'D'), ('HA1',311,'H'), ('HA2',77,'V'), ('HA2',155,'E'), ('HA2',160,'N')],
						   "3c2":   [('HA1',144,'N'), ('HA1',159,'F'), ('HA1',225,'N'), ('HA2',160,'N'), ('HA1',142,'R')],
						   "3c3.b": [('HA1',83,'R'), ('HA1',261,'Q'), ('HA1',62,'K'),  ('HA1',122,'D')]
							},
	'epitope_masks_fname':'source-data/H3N2_epitope_masks.tsv',
	'epitope_mask_version':'wolf',
	'HI_fname':'data/h3n2_hi_titers.tsv',
	'html_vars': {'coloring': 'ep, ne, rb, lbi, dfreq, region, date, cHI, HI_dist',
				   'gtplaceholder': 'HA1 positions...',
					'freqdefault': '3c2.a, 3c3.a, 3c3.b'},
	'js_vars': {'LBItau': 0.0005, 'LBItime_window': 0.5, 'dfreq_dn':2},
	'excluded_tables': ['NIMR_Sep2012_08.csv'], #, 'nimr-sep-2010-table8', 'nimr-sep-2010-table8','NIMR_Sep2012_11.csv'],
	'layout':'auspice',
	'min_aamuts': 1,
	'predictors': ['ep']														# estimate
#	'predictors': { 'dfreq': [2.50, 2.84], 'cHI': [1.68, 0.45] }				# fix predictor: [value, std deviation]
	})


class H3N2_filter(flu_filter):
	def __init__(self,min_length = 900, **kwargs):
		'''
		parameters
		min_length  -- minimal length for a sequence to be acceptable
		'''
		flu_filter.__init__(self, **kwargs)
		self.min_length = min_length
		self.vaccine_strains =[
				{
					"strain": "A/Wisconsin/67/2005",
					"db": "IRD",
					"accession": "CY163984",
					"date": "2005-08-31",
					"region": "north_america",
					"country": "usa",
					"seq": "ATGAAGACTATCATTGCTTTGAGCTACATTCTATGTCTGGTTTTCGCTCAAAAACTTCCCGGAAATGACAACAGCACGGCAACGCTGTGCCTTGGGCACCATGCAGTACCAAACGGAACGATAGTGAAAACAATCACGAATGACCAAATTGAAGTTACTAATGCTACTGAGCTGGTTCAGAGTTCCTCAACAGGTGGAATATGCGACAGTCCTCATCAGATCCTTGATGGAGAAAACTGCACACTAATAGATGCTCTATTGGGAGACCCTCAGTGTGATGGCTTCCAAAATAAGAAATGGGACCTTTTTGTTGAACGCAGCAAAGCCTACAGCAACTGTTACCCTTATGATGTGCCGGATTATGCCTCCCTTAGGTCACTAGTTGCCTCATCCGGCACACTGGAGTTTAACGATGAAAGCTTCAATTGGACTGGAGTCACTCAAAATGGAACAAGCTCTTCTTGCAAAAGGAGATCTAATAACAGTTTCTTTAGTAGATTGAATTGGTTGACCCACTTAAAATTCAAATACCCAGCATTGAACGTGACTATGCCAAACAATGAAAAATTTGACAAATTGTACATTTGGGGGGTTCACCACCCGGTTACGGACAATGACCAAATCTTCCTGTATGCTCAAGCATCAGGAAGAATCACAGTCTCTACCAAAAGAAGCCAACAAACTGTAATCCCGAATATCGGATCTAGACCCAGAATAAGGAATATCCCCAGCAGAATAAGCATCTATTGGACAATAGTAAAACCGGGAGACATACTTTTGATTAACAGCACAGGGAATCTAATTGCTCCTAGGGGTTACTTCAAAATACGAAGTGGGAAAAGCTCAATAATGAGATCAGATGCACCCATTGGCAAATGCAATTCTGAATGCATCACTCCAAATGGAAGCATTCCCAATGACAAACCATTTCAAAATGTAAACAGGATCACATATGGGGCCTGTCCCAGATATGTTAAGCAAAACACTCTGAAATTGGCAACAGGGATGCGAAATGTACCAGAGAAACAAACTAGAGGCATATTTGGCGCAATCGCGGGTTTCATAGAAAATGGTTGGGAGGGAATGGTGGATGGTTGGTACGGTTTCAGGCATCAAAATTCTGAGGGAATAGGACAAGCAGCAGATCTCAAAAGCACTCAAGCAGCAATCAATCAAATCAATGGGAAGCTGAATAGGTTGATCGGGAAAACCAACGAGAAATTCCATCAGATTGAAAAAGAATTCTCAGAAGTAGAAGGGAGAATTCAGGACCTCGAGAAATATGTTGAGGACACTAAAATAGATCTCTGGTCATACAACGCGGAGCTTCTTGTTGCCCTGGAGAACCAACATACAATTGATCTAACTGACTCAGAAATGAACAAACTGTTTGAAAGAACAAAGAAGCAACTGAGGGAAAATGCTGAGGATATGGGCAATGGTTGTTTCAAAATATACCACAAATGTGACAATGCCTGCATAGGATCAATCAGAAATGGAACTTATGACCATGATGTATACAGAGATGAAGCATTAAACAACCGGTTCCAGATCAAAGGCGTTGAGCTGAAGTCAGGATACAAAGATTGGATCCTATGGATTTCCTTTGCCATATCATGTTTTTTGCTTTGTGTTGCTTTGTTGGGGTTCATCATGTGGGCCTGCCAAAAAGGCAACATTAGGTGCAACATTTGCATTTGA"
				},	{
					"strain": "A/Brisbane/10/2007",
					"db": "IRD",
					"accession": "CY113005",
					"date": "2007-02-06",
					"region": "oceania",
					"country": "australia",
					"seq": "ATGAAGACTATCATTGCTTTGAGCTACATTCTATGTCTGGTTTTCACTCAAAAACTTCCCGGAAATGACAACAGCACGGCAACGCTGTGCCTTGGGCACCATGCAGTACCAAACGGAACGATAGTGAAAACAATCACGAATGACCAAATTGAAGTTACTAATGCTACTGAGCTGGTTCAGAGTTCCTCAACAGGTGAAATATGCGACAGTCCTCATCAGATCCTTGATGGAGAAAACTGCACACTAATAGATGCTCTATTGGGAGACCCTCAGTGTGATGGCTTCCAAAATAAGAAATGGGACCTTTTTGTTGAACGCAGCAAAGCCTACAGCAACTGTTACCCTTATGATGTGCCGGATTATGCCTCCCTTAGGTCACTAGTTGCCTCATCCGGCACACTGGAGTTTAACAATGAAAGCTTCAATTGGACTGGAGTCACTCAAAACGGAACAAGCTCTGCTTGCATAAGGAGATCTAATAACAGTTTCTTTAGTAGATTGAATTGGTTGACCCACTTAAAATTCAAATACCCAGCATTGAACGTGACTATGCCAAACAATGAAAAATTTGACAAATTGTACATTTGGGGGGTTCACCACCCGGGTACGGACAATGACCAAATCTTCCCGTATGCTCAAGCATCAGGAAGAATCACAGTCTCTACCAAAAGAAGCCAACAAACTGTAATCCCGAATATCGGATCTAGACCCAGAGTAAGGAATATCCCCAGCAGAATAAGCATCTATTGGACAATAGTAAAACCGGGAGACATACTTTTGATTAACAGCACAGGGAATCTAATTGCTCCTAGGGGTTACTTCAAAATACGAAGTGGGAAAAGCTCAATAATGAGATCAGATGCACCCATTGGCAAATGCAATTCTGAATGCATCACTCCAAACGGAAGCATTCCCAATGACAAACCATTCCAAAATGTAAACAGGATCACATACGGGGCCTGTCCCAGATATGTTAAGCAAAACACTCTGAAATTGGCAACAGGGATGCGAAATGTACCAGAGAAACAAACTAGAGGCATATTTGGCGCAATCGCGGGTTTCATAGAAAATGGTTGGGAGGGAATGGTGGATGGTTGGTACGGTTTCAGGCATCAAAATTCTGAGGGAATAGGACAAGCAGCAGATCTCAAAAGCACTCAAGCAGCAATCGATCAAATCAATGGGAAGCTGAATAGGTTGATCGGGAAAACCAACGAGAAATTCCATCAGATTGAAAAAGAATTCTCAGAAGTCGAAGGGAGAATTCAGGACCTTGAGAAATATGTTGAGGACACCAAAATAGATCTCTGGTCATACAACGCGGAGCTTCTTGTTGCCCTGGAGAACCAACATACAATTGATCTAACTGACTCAGAAATGAACAAACTGTTTGAAAAAACAAAGAAGCAACTGAGGGAAAATGCTGAGGATATGGGCAATGGTTGTTTCAAAATATACCACAAATGTGACAATGCCTGCATAGGATCAATCAGAAATGGAACTTATGACCACAATGTATACAGAGATGAAGCATTAAACAACCGGTTCCAGATCAAGGGCGTTGAGCTGAAGTCAGGATACAAAGATTGGATCCTATGGATTTCCTTTGCCATATCATGTTTTTTGCTTTGTGTTGCTTTGTTGGGGTTCATCATGTGGGCCTGCCAAAAAGGCAACATTAGGTGCAACATTTGCATTTGA"
				},	{
					"strain": "A/Perth/16/2009",
					"db": "IRD",
					"accession": "GQ293081",
					"date": "2009-04-07",
					"region": "oceania",
					"country": "australia",
					"seq": "ATGAAGACTATCATTGCTTTGAGCTACATTCTATGTCTGGTTTTCGCTCAAAAACTTCCTGGAAATGACAACAGCACGGCAACGCTGTGCCTTGGGCACCATGCAGTACCAAACGGAACGATAGTGAAAACAATCACGAATGACCAAATTGAAGTTACTAATGCTACTGAGCTGGTTCAGAGTTCCTCAACAGGTGAAATATGCGACAGTCCTCATCAGATCCTTGATGGAAAAAACTGCACACTAATAGATGCTCTATTGGGAGACCCTCAGTGTGATGGCTTCCAAAATAAGAAATGGGACCTTTTTGTTGAACGCAGCAAAGCCTACAGCAACTGTTACCCTTATGATGTGCCGGATTATGCCTCCCTTAGGTCACTAGTTGCCTCATCCGGCACACTGGAGTTTAACAATGAAAGCTTCAATTGGACTGGAGTCACTCAAAACGGAACAAGCTCTGCTTGCATAAGGAGATCTAAAAACAGTTTCTTTAGTAGATTGAATTGGTTGACCCACTTAAACTTCAAATACCCAGCATTGAACGTGACTATGCCAAACAATGAACAATTTGACAAATTGTACATTTGGGGGGTTCACCACCCGGGTACGGACAAAGACCAAATCTTCCTGTATGCTCAAGCATCAGGAAGAATCACAGTCTCTACCAAAAGAAGCCAACAAACCGTAAGCCCGAATATCGGATCTAGACCCAGAGTAAGGAATATCCCTAGCAGAATAAGCATCTATTGGACAATAGTAAAACCGGGAGACATACTTTTGATTAACAGCACAGGGAATCTAATTGCTCCTAGGGGTTACTTCAAAATACGAAGTGGGAAAAGCTCAATAATGAGATCAGATGCACCCATTGGCAAATGCAATTCTGAATGCATCACTCCAAATGGAAGCATTCCCAATGACAAACCATTCCAAAATGTAAACAGGATCACATACGGGGCCTGTCCCAGATATGTTAAGCAAAACACTCTGAAATTGGCAACAGGGATGCGAAATGTACCAGAGAAACAAACTAGAGGCATATTTGGCGCAATCGCGGGTTTCATAGAAAATGGTTGGGAGGGAATGGTGGATGGTTGGTACGGTTTCAGGCATCAAAATTCTGAGGGAAGAGGACAAGCAGCAGATCTCAAAAGCACTCAAGCAGCAATCGATCAAATCAATGGGAAGCTGAATAGATTGATCGGGAAAACCAACGAGAAATTCCATCAGATTGAAAAAGAATTCTCAGAAGTCGAAGGGAGAATTCAGGACCTTGAGAAATATGTTGAGGACACTAAAATAGATCTCTGGTCATACAACGCGGAGCTTCTTGTTGCCCTGGAGAACCAACATACAATTGATCTAACTGACTCAGAAATGAACAAACTGTTTGAAAAAACAAAGAAGCAACTGAGGGAAAATGCTGAGGATATGGGCAATGGTTGTTTCAAAATATACCACAAATGTGACAATGCCTGCATAGGATCAATCAGAAATGGAACTTATGACCACGATGTATACAGAGATGAAGCATTAAACAACCGGTTTCAGATCAAGGGAGTTGAGCTGAAGTCAGGGTACAAAGATTGGATCCTATGGATTTCCTTTGCCATATCATGTTTTTTGCTTTGTGTTGCTTTGTTGGGGTTCATCATGTGGGCCTGCCAAAAAGGCAACATTAGGTGCAACATTTGCATTTGA"
				},	{
					"strain": "A/Victoria/361/2011",
					"db": "IRD",
					"accession": "GQ293081",
					"date": "2011-10-24",
					"region": "oceania",
					"country": "australia",
					"seq": "ATGAAGACTATCATTGCTTTGAGCCACATTCTATGTCTGGTTTTCGCTCAAAAACTTCCTGGAAATGACAACAGCACGGCAACGCTGTGCCTTGGGCACCATGCAGTACCAAACGGAACGATAGTGAAAACAATCACGAATGACCAAATTGAAGTTACTAATGCTACTGAGCTGGTTCAGAATTCCTCAATAGGTGAAATATGCGACAGTCCTCATCAGATCCTTGATGGAGAAAACTGCACACTAATAGATGCTCTATTGGGAGACCCTCAGTGTGATGGCTTCCAAAATAAGAAATGGGACCTTTTTGTTGAACGAAGCAAAGCCTACAGCAACTGTTACCCTTATGATGTGCCGGATTATGCCTCCCTTAGGTCACTAGTTGCCTCATCCGGCACACTGGAGTTTAACAATGAAAGCTTCAATTGGACTGGAGTCACTCAAAACGGAACAAGTTCTGCTTGCATAAGGAGATCTAATAATAGTTTCTTTAGTAGATTAAATTGGTTGACCCGCTTAAACTTCAAATACCCAGCATTGAACGTGACTATGCCAAACAATGAACAATTTGACAAATTGTACATTTGGGGGGTTCACCACCCGGTTACGGACAAGGAACAAATCTTCCTGTATGCTCAATCATCAGGAAGAATCACAGTATCTACCAAAAGAAGCCAACAAGCTGTAATCCCGAATATCGGATATAGACCCAGAATAAGGAATATCCCTAGCAGAATAAGCATCTATTGGACAATAGTAAAACCGGGAGACATACTTTTGATTAACAGCACAGGGAATCTAATTGCTCCTAGGGGTTACTTCAAAATACGAAGTGGGAAAAGCTCAATAATGAGATCAGATGCACCCATTGGCAAATGCAATTCTGAATGCATCACTCCAAATGGAAGCATTCCCAATGACAAACCATTCCAAAATGTAAACAGGATCACATACGGGGCCTGTCCCAGATATGTTAAGCAAAGCACTCTGAAATTGGCAACAGGAATGCGAAATGTACCAGAGAAACAAACTAGAGGCATATTTGGCGCAATAGCGGGTTTCATAGAAAATGGTTGGGAGGGAATGGTGGATGGTTGGTACGGTTTCAGGCATCAAAATTCTGAGGGAAGAGGACAAGCAGCAGATCTCAAAAGCACTCAAGCAGCAATCGATCAAATCAATGGGAAGCTGAATCGATTGATCGGGAAAACCAACGAGAAATTCCATCAGATTGAAAAAGAATTCTCAGAAGTCGAAGGGAGAATTCAGGACCTTGAGAAATATGTTGAGGACACTAAAATAGATCTCTGGTCATACAACGCGGAGCTTCTTGTTGCCCTGGAGAACCAACATACAATTGATCTAACTGACTCAGAAATGAACAAACTGTTTGAAAAAACAAAGAAGCAACTAAGGGAAAATGCTGAGGATATGGGCAATGGTTGTTTCAAAATATACCACAAATGTGACAATGCCTGCATAGGATCAATCAGAAATGGAACTTATGACCACGATGTATACAGAGATGAAGCATTAAACAACCGGTTCCAGATCAAGGGAGTTGAGCTGAAGTCAGGGTACAAAGATTGGATCCTATGGATTTCCTTTGCCATATCATGTTTTTTGCTTTGTGTTGCTTTGTTGGGGTTCATCATGTGGGCCTGCCAAAAGGGCAACATTAGGTGCAACATTTGCATTTGA"
				},	{
					"strain": "A/Texas/50/2012",
					"db": "GISAID",
					"isolate_id": "EPI_ISL_129858",
					"date": "2012-04-15",
					"region": "north_america",
					"country": "usa",
					"seq": "ATGAAGACTATCATTGCTTTGAGCTACATTCTATGTCTGGTTTTCGCTCAAAAACTTCCTGGAAATGACAATAGCACGGCAACGCTGTGCCTTGGGCACCATGCAGTACCAAACGGAACGATAGTGAAAACAATCACGAATGACCGAATTGAAGTTACTAATGCTACTGAACTGGTTCAGAATTCCTCAATAGGTGAAATATGCGACAGTCCTCATCAGATCCTTGATGGAGAAAACTGCACACTAATAGATGCTCTATTGGGAGACCCTCAGTGTGATGGCTTCCAAAATAAGAAATGGGACCTTTTTGTTGAACGAAGCAAAGCCTACAGCAACTGTTACCCTTATGATGTGCCGGATTATGCCTCCCTTAGGTCACTAGTTGCCTCATCCGGCACACTGGAGTTTAACAATGAAAGCTTCAATTGGAATGGAGTCACTCAAAACGGAACAAGTTCTGCTTGCATAAGGAGATCTAATAATAGTTTCTTTAGTAGATTAAATTGGTTGACCCACTTAAACTTCAAATACCCAGCATTGAACGTGACTATGCCAAACAATGAACAATTTGACAAATTGTACATTTGGGGGGTTCACCACCCGGGTACGGACAAGGACCAAATCTTCCTGTATGCTCAACCATCAGGAAGAATCACAGTATCTACCAAAAGAAGCCAACAAGCTGTAATCCCGAATATCGGATCTAGACCCAGAATAAGGAATATCCCTAGCAGAATAAGCATCTATTGGACAATAGTAAAACCGGGAGACATACTTTTGATTAACAGCACAGGGAATCTAATTGCTCCTAGGGGTTACTTCAAAATACGAAGTGGGAAAAGCTCAATAATGAGATCAGATGCACCCATTGGCAAATGCAAGTCTGAATGCATCACTCCAAATGGAAGCATTCCCAATGACAAACCATTCCAAAATGTAAACAGGATCACATACGGGGCCTGTCCCAGATATGTTAAGCAAAGCACTCTGAAATTGGCAACAGGAATGCGGAATGTACCAGAGAAACAAACTAGAGGCATATTTGGCGCAATAGCGGGTTTCATAGAAAATGGTTGGGAGGGAATGGTGGATGGTTGGTACGGTTTCAGGCATCAAAATTCTGAGGGAAGAGGACAAGCAGCAGATCTCAAAAGCACTCAAGCAGCAATCGATCAAATCAATGGGAAGCTGAATCGATTGATCGGGAAAACCAACGAGAAATTCCATCAGATTGAAAAAGAATTCTCAGAAGTAGAAGGGAGAATTCAGGACCTTGAGAAATATGTTGAGGACACTAAAATAGATCTCTGGTCATACAACGCGGAGCTTCTTGTTGCCCTGGAGAACCAACATACAATTGATCTAACTGACTCAGAAATGAACAAACTGTTTGAAAAAACAAAGAAGCAACTGAGGGAAAATGCTGAGGATATGGGCAATGGTTGTTTCAAAATATACCACAAATGTGACAATGCCTGCATAGGATCAATCAGAAATGGAACTTATGACCACGATGTATACAGAGATGAAGCATTAAACAACCGGTTCCAGATCAAGGGAGTTGAGCTGAAGTCAGGGTACAAAGATTGGATCCTATGGATTTCCTTTGCCATATCATGTTTTTTGCTTTGTGTTGCTTTGTTGGGGTTCATCATGTGGGCCTGCCAAAAGGGCAACATTAGGTGCAACATTTGCATTTGA",
				},	{
					"strain": "A/Switzerland/9715293/2013",
					"db": "GISAID",
					"isolate_id": "EPI_ISL_162149",
					"date": "2013-12-06",
					"region": "europe",
					"country": "switzerland",
					"seq": "ATGAAGACTATCATTGCTTTGAGCTACATTCTATGTCTGGTTTTCGCTCAAAAACTTCCTGGAAATGACAATAGCACGGCAACGCTGTGCCTTGGGCACCATGCAGTACCAAACGGAACGATAGTGAAAACAATCACGAATGACCGAATTGAAGTTACTAATGCTACTGAGCTGGTTCAGAATTCCTCAATAGGTGAAATATGCGACAGTCCTCATCAGATCCTTGATGGAGAAAACTGCACACTAATAGATGCTCTATTGGGAGACCCTCAGTGTGATGGCTTTCAAAATAAGAAATGGGACCTTTTTGTTGAACGAAGCAAAGCCTACAGCAACTGTTACCCTTATGATGTGCCGGATTATGCCTCCCTTAGGTCACTAGTTGCCTCATCCGGCACACTGGAGTTTAACAATGAAAGCTTCAATTGGGCTGGAGTCACTCAAAACGGAACAAGTTCTTCTTGCATAAGGGGATCTAATAGTAGTTTCTTTAGTAGATTAAATTGGTTGACCCACTTAAACTCCAAATACCCAGCATTAAACGTGACTATGCCAAACAATGAACAATTTGACAAATTGTACATTTGGGGGGTTCACCACCCGGGTACGGACAAGGACCAAATCTTCCTGTATGCACAATCATCAGGAAGAATCACAGTATCTACCAAAAGAAGCCAACAAGCTGTAATCCCGAATATCGGATCTAGACCCAGAATAAGGGATATCCCTAGCAGAATAAGCATCTATTGGACAATAGTAAAACCGGGAGACATACTTTTGATTAACAGCACAGGGAATCTAATTGCTCCTAGGGGTTACTTCAAAATACGAAGTGGGAAAAGCTCAATAATGAGATCAGATGCACCCATTGGCAAATGCAAGTCTGAATGCATCACTCCAAATGGAAGCATTCCCAATGACAAACCATTCCAAAATGTAAACAGGATCACATACGGGGCCTGTCCCAGATATGTTAAGCAAAGCACTCTGAAATTGGCAACAGGAATGCGAAATGTACCAGAGAGACAAACTAGAGGCATATTTGGCGCAATAGCGGGTTTCATAGAAAATGGTTGGGAGGGAATGGTGGATGGTTGGTACGGCTTCAGGCATCAAAATTCTGAGGGAAGAGGACAAGCAGCAGATCTCAAAAGCACTCAAGCAGCAATCGATCAAATCAATGGGAAGCTGAATCGATTGATCGGGAAAACCAACGAGAAATTCCATCAGATTGAAAAAGAATTCTCAGAAGTAGAAGGGAGAATTCAGGACCTTGAGAAATATGTTGAGGACACAAAAATAGATCTCTGGTCATACAACGCGGAGCTTCTTGTTGCCCTGGAGAACCAACATACAATTGATCTAACTGACTCAGAAATGAACAAACTGTTTGAAAAAACAAAGAAGCAACTGAGGGAAAATGCTGAGGATATGGGCAATGGTTGTTTCAAAATATACCACAAATGTGACAATGCCTGCATAGGATCAATCAGAAATGGAACTTATGACCACGATGTATACAGGGATGAAGCATTAAACAACCGGTTCCAGATCAAGGGAGTTGAGCTGAAGTCAGGGTACAAAGATTGGATCCTATGGATTTCCTTTGCCATATCATGTTTTTTGCTTTGTGTTGCTTTGTTGGGGTTCATCATGTGGGCCTGCCAAAAGGGCAACATTAGGTGCAACATTTGCATTTGA",
				},  {
					"strain": "A/HongKong/4801/2014",
					"db": "GISAID",
					"isolate_id": "EPI_ISL_165554",
					"date": "2014-02-26",
					"region": "china",
					"country": "hong_kong",
					"seq": "ATGAAGACTATCATTGCTTTGAGCTACATTCTATGTCTGGTTTTCGCTCAAAAAATTCCTGGAAATGACAATAGCACGGCAACGCTGTGCCTTGGGCACCATGCAGTACCAAACGGAACGATAGTGAAAACAATCACGAATGACCGAATTGAAGTTACTAATGCTACTGAGCTGGTTCAGAATTCCTCAATAGGTGAAATATGCGACAGTCCTCATCAGATCCTTGATGGAGAAAACTGCACACTAATAGATGCTCTATTGGGAGACCCTCAGTGTGATGGCTTTCAAAATAAGAAATGGGACCTTTTTGTTGAACGAAGCAAAGCCTACAGCAACTGTTACCCTTATGATGTGCCGGATTATGCCTCCCTTAGGTCACTAGTTGCCTCATCCGGCACACTGGAGTTTAACAATGAAAGCTTCAATTGGACTGGAGTCACTCAAAACGGAACAAGTTCTGCTTGCATAAGGAGATCTAGTAGTAGTTTCTTTAGTAGATTAAATTGGTTGACCCACTTAAACTACACATACCCAGCATTGAACGTGACTATGCCAAACAATGAACAATTTGACAAATTGTACATTTGGGGGGTTCACCACCCGGGTACGGACAAGGACCAAATCTTCCTGTATGCTCAATCATCAGGAAGAATCACAGTATCTACCAAAAGAAGCCAACAAGCTGTAATCCCAAATATCGGATCTAGACCCAGAATAAGGGATATCCCTAGCAGAATAAGCATCTATTGGACAATAGTAAAACCGGGAGACATACTTTTGATTAACAGCACAGGGAATCTAATTGCTCCTAGGGGTTACTTCAAAATACGAAGTGGGAAAAGCTCAATAATGAGATCAGATGCACCCATTGGCAAATGCAAGTCTGAATGCATCACTCCAAATGGAAGCATTCCCAATGACAAACCATTCCAAAATGTAAACAGGATCACATACGGGGCCTGTCCCAGATATGTTAAGCATAGCACTCTGAAATTGGCAACAGGAATGCGAAATGTACCAGAGAAACAAACTAGAGGCATATTTGGCGCAATAGCGGGTTTCATAGAAAATGGTTGGGAGGGAATGGTGGATGGTTGGTACGGTTTCAGGCATCAAAATTCTGAGGGAAGAGGACAAGCAGCAGATCTCAAAAGCACTCAAGCAGCAATCGATCAAATCAATGGGAAGCTGAATCGATTGATCGGGAAAACCAACGAGAAATTCCATCAGATTGAAAAAGAATTCTCAGAAGTAGAAGGAAGAATTCAGGACCTTGAGAAATATGTTGAGGACACTAAAATAGATCTCTGGTCATACAACGCGGAGCTTCTTGTTGCCCTGGAGAACCAACATACAATTGATCTAACTGACTCAGAAATGAACAAACTGTTTGAAAAAACAAAGAAGCAACTGAGGGAAAATGCTGAGGATATGGGCAATGGTTGTTTCAAAATATACCACAAATGTGACAATGCCTGCATAGGATCAATAAGAAATGGAACTTATGACCACAATGTGTACAGGGATGAAGCATTAAACAACCGGTTCCAGATCAAGGGAGTTGAGCTGAAGTCAGGGTACAAAGATTGGATCCTATGGATTTCCTTTGCCATATCATGTTTTTTGCTTTGTGTTGCTTTGTTGGGGTTCATCATGTGGGCCTGCCAAAAGGGCAACATTAGGTGCAACATTTGCATTTGA",
				},  {
					"strain": "A/Alaska/232/2015",
					"db": "GISAID",
					"isolate_id": "EPI787411",
					"date": "2015-09-09",
					"region": "north_america",
					"country": "usa",
					"seq": "GGATAATTCTATTAACCATGAAGACTATCATTGCTTTGAGCTACATTCTATGTCTGGTTTTCGCTCAAAAAATTCCTGGAAATGACAATAGCACGGCAACGCTGTGCCTTGGGCACCATGCAGTACCAAACGGAACGATAGTGAAAACAATCACAAATGACCGAATTGAAGTTACTAATGCTACTGAGTTGGTTCAGAATTCCTCAATAGGTGAAATATGCGACAGTCCTCATCAGATCCTTGATGGAGAGAACTGCACACTAATAGATGCTCTATTGGGAGACCCTCAGTGTGATGGCTTTCAAAATAAGAAATGGGACCTTTTTGTTGAACGAAGCAAAGCCTACAGCAACTGTTACCCTTATGATGTGCCGGATTATGCCTCCCTTAGGTCACTAGTTGCCTCATCCGGCACACTGGAGTTTAACAATGAAAGCTTCAATTGGACTGGAGTCACTCAAAACGGAACAAGTTCTTCTTGCATAAGGAGATCTAGTAGTAGTTTCTTTAGTAGATTAAATTGGTTGACCCACTTAAACTACAAATATCCAGCATTGAACGTGACTATGCCAAACAAGGAACAATTTGACAAATTGTACATTTGGGGGGTTCACCACCCGGGTACGGACAAGGACCAAATCTACCCGTATGCTCAATCATCAGGAAGAATCACAGTATCTACCAAAAGAAGCCAACAAGCTGTAATCCCAAATATCGGATCTAGACCCAGAATAAGGGATATCCCTAGCAGAATAAGCATCTATTGGACAATAGTAAAACCGGGAGACATACTTTTGATTAACAGCACAGGGAATCTAATTGCTCCTAGGGGTTACTTCAAAATACGAAGTGGGAAAAGCTCAATAATGAGATCAGATGCACCCATTGGCAAATGCAAGTCTGAATGCATCACTCCAAATGGAAGCATTCCCAATGACAAACCATTCCAAAATGTAAACAGGATCACATACGGGGCCTGTCCCAGATATGTTAAGCATAGCACTCTGAAATTGGCAACAGGAATGCGAAATGTACCAGAGAAACAAACTAGAGGCATATTTGGCGCAATAGCGGGTTTCATAGAAAATGGTTGGGAGGGAATGGTGGATGGTTGGTACGGTTTCAGGCATCAAAATTCTGAGGGAAGAGGACAAGCAGCAGATCTCAAAAGCACTCAAGCAGCAATCGATCAAATCAATGGGAAGCTGAATCGGTTGATCGGGAAAACCAACGAGAAATTCCATCAGATTGAAAAAGAATTCTCAGAAGTAGAAGGAAGAGTTCAAGACCTTGAGAAATATGTTGAGGACACTAAAATAGATCTCTGGTCATACAACGCGGAGCTTCTTGTTGCCCTGGAGAACCAACATACAATTGATCTAACTGACTCAGAAATGAACAAACTGTTTGAAAAAACAAAGAAGCAACTGAGGGAAAATGCTGAGGATATGGGAAATGGTTGTTTCAAAATATACCACAAATGTGACAATGCCTGCATAGGATCAATAAGAAATGAAACTTATGACCACAATGTGTACAGGGATGAAGCATTAAACAACCGGTTCCAGATCAAGGGAGTTGAGCTGAAGTCAGGGTACAAAGATTGGATCCTATGGATTTCCTTTGCCATATCATGTTTTTTGCTTTGTGTTGCTTTGTTGGGGTTCATCATGTGGGCCTGCCAAAAGGGCAACATTAGATGCAACATTTGCATTTGAGTGCATTAATTAAAAACAC"
				}
			]
		tmp_outgroup = SeqIO.read('source-data/H3N2_outgroup.gb', 'genbank')
		genome_annotation = tmp_outgroup.features
		self.cds = {x.qualifiers['gene'][0]:x for x in genome_annotation
				if 'gene' in x.qualifiers and x.type=='CDS' and
				x.qualifiers['gene'][0] in ['SigPep', 'HA1', 'HA2']}
		self.outgroup = {
			'strain': 'A/Beijing/32/1992',
			'db': 'IRD',
			'accession': 'U26830',
			'date': '1992-01-01',
			'country': 'china',
			'region': 'china',
			'seq': str(tmp_outgroup.seq).upper()
		}

class H3N2_clean(virus_clean):
	def __init__(self,**kwargs):
		virus_clean.__init__(self, **kwargs)

	def clean_outbreaks(self):
		"""Remove duplicate strains, where the geographic location, date of sampling and sequence are identical"""
		virus_hashes = set()
		new_viruses = []
		for v in self.viruses:
			try:
				geo = re.search(r'A/([^/]+)/', v.strain).group(1)
			except:
				print "clean outbreaks:, couldn't parse geo of ",v.strain
				continue
			if geo:
				vhash = (geo, v.date, str(v.seq))
				if vhash not in virus_hashes:
					new_viruses.append(v)
					virus_hashes.add(vhash)

		self.viruses = MultipleSeqAlignment(new_viruses)

	def clean_outliers(self):
		"""Remove single outlying viruses"""
		new_viruses = []
		outlier_strains = [
			"A/Chile/8266/2003", "A/Louisiana/4/2003", "A/Lousiana/4/2003", "A/OSAKA/31/2005",
			"A/Sari/388/2006", "A/HongKong/HK1/2008", "A/HongKong/HK1MA21-1/2008",
			"A/HongKong/HK1MA21-2/2008", "A/HongKong/HK2/2008", "A/HongKong/HK2MA21-1/2008",
			"A/HongKong/HK2MA21-2/2008", "A/HongKong/HK2MA21-3/2008", "A/HongKong/HK4/2008",
			"A/HongKong/HK5/2008", "A/HongKong/HK5MA21-1/2008", "A/HongKong/HK5MA21-3/2008",
			"A/HongKong/HK6/2008", "A/HongKong/HK6MA21-2/2008", "A/HongKong/HK6MA21-3/2008",
			"A/HongKong/HK7/2008", "A/HongKong/HK8/2008", "A/HongKong/HK8MA21-2/2008",
			"A/HongKong/HK8MA21-3/2008", "A/HongKong/HK8MA21-4/2008", "A/HongKong/HK9/2008",
			"A/HongKong/HK9MA21-1/2008", "A/HongKong/HK9MA21-2/2008", "A/HongKong/HK9MA21-3/2008",
			"A/HongKong/HK10/2008", "A/HongKong/HK10MA21-1/2008", "A/HongKong/HK10MA21-2/2008",
			"A/HongKong/HK10MA21-4/2008", "A/HongKong/HK11MA21-1/2008", "A/HongKong/HK11MA21-4/2008",
			"A/HongKong/HK12/2008", "A/HongKong/HK12MA21-2/2008", "A/HongKong/HKMA12/2008",
			"A/HongKong/HKMA12A/2008", "A/HongKong/HKMA12B/2008", "A/HongKong/HKMA12D/2008",
			"A/HongKong/HKMA12E/2008", "A/HongKong/HKMA20B/2008", "A/HongKong/HKMA20E/2008",
			"A/Busan/15453/2009", "A/Pennsylvania/14/2010", "A/Pennsylvania/40/2010",
			"A/Guam/AF2771/2011", "A/Indiana/8/2011", "A/Kenya/155/2011", "A/Kenya/168/2011",
			"A/Kenya/170/2011", "A/Nepal/142/2011", "A/Pennsylvania/09/2011", "A/Pennsylvania/9/2011",
			"A/Quebec/167936/2011", "A/Quebec/170658/2011", "A/India/6352/2012", "A/Indiana/08/2012",
			"A/Indiana/13/2012", "A/Ohio/34/2012", "A/Helsinki/942/2013", "A/Indiana/06/2013",
			"A/Indiana/11/2013", "A/Indiana/21/2013", "A/Jiangsu-Tianning/1707/2013",
			"A/HuNan/01/2014", "A/Jiangsu-Chongchuan/1830/2014", "A/Jiangsu-Chongchuan/12179/2014",
			"A/Ohio/2/2014", "A/Ohio/4319/2014", "A/SaoPaulo/3-34891/2014", "A/NewJersey/53/2015",
			"A/SaoPaulo/36178/2015", "A/SaoPaulo/61282/2015", "A/SaoPaulo/65194/2015",
			"A/Sydney/53/2015", "A/Michigan/82/2016", "A/Michigan/83/2016", "A/Michigan/84/2016",
			"A/Michigan/87/2016", "A/Michigan/89/2016", "A/Michigan/90/2016", "A/Michigan/91/2016",
			"A/Michigan/93/2016", "A/Michigan/94/2016", "A/Michigan/95/2016", "A/Michigan/96/2016",
			"A/Ohio/27/2016", "A/Ohio/28/2016", "A/Ohio/32/2016", "A/Ohio/33/2016", "A/Ohio/35/2016",
			"A/Zhejiang-Wuxin/1300/2016"]
		for v in self.viruses:
			if v.strain in outlier_strains:
				if self.verbose > 1:
					print "\tremoving", v.strain
			else:
				new_viruses.append(v)
		self.viruses = MultipleSeqAlignment(new_viruses)

	def clean(self):
		self.clean_generic()
		self.clean_outbreaks()
		print "Number of viruses after outbreak filtering:",len(self.viruses)
		self.clean_outliers()
		print "Number of viruses after outlier filtering:",len(self.viruses)

class H3N2_refine(tree_refine):
	def __init__(self, **kwargs):
		tree_refine.__init__(self, **kwargs)
		self.epitope_mask = ""
		if "epitope_masks_fname" in self.kwargs and "epitope_mask_version" in self.kwargs:
			epitope_map = {}
			with open(self.kwargs["epitope_masks_fname"]) as f:
				for line in f:
					(key, value) = line.split()
					epitope_map[key] = value
			if self.kwargs["epitope_mask_version"] in epitope_map:
				self.epitope_mask = epitope_map[self.kwargs["epitope_mask_version"]]
		self.epitope_mask = np.fromstring(self.epitope_mask, dtype='S1')				# epitope_mask is numpy array

	def refine(self):
		self.refine_generic()  # -> all nodes now have aa_seq, xvalue, yvalue, trunk, and basic virus properties
		self.add_H3N2_attributes()

	def epitope_sites(self, aa):
		aaa = np.fromstring(aa, 'S1')
		return ''.join(aaa[self.epitope_mask[:len(aa)]=='1'])

	def nonepitope_sites(self, aa):
		aaa = np.fromstring(aa, 'S1')
		return ''.join(aaa[self.epitope_mask[:len(aa)]=='0'])

	def receptor_binding_sites(self, aa):
		'''
		Receptor binding site mutations from Koel et al. 2014
		These are (145, 155, 156, 158, 159, 189, 193) in canonical HA numbering
		need to subtract one since python arrays start at 0
		'''
		sp = 16
		rbs = map(lambda x:x+sp-1, [145, 155, 156, 158, 159, 189, 193])
		return ''.join([aa[pos] for pos in rbs])

	def get_total_peptide(self, node):
		'''
		the concatenation of signal peptide, HA1, HA1
		'''
		return node.aa_seq['SigPep']+node.aa_seq['HA1'] + node.aa_seq['HA2']

	def epitope_distance(self, aaA, aaB):
		"""Return distance of sequences aaA and aaB by comparing epitope sites"""
		epA = self.epitope_sites(aaA)
		epB = self.epitope_sites(aaB)
		distance = sum(a != b for a, b in izip(epA, epB))
		return distance

	def nonepitope_distance(self, aaA, aaB):
		"""Return distance of sequences aaA and aaB by comparing non-epitope sites"""
		neA = self.nonepitope_sites(aaA)
		neB = self.nonepitope_sites(aaB)
		distance = sum(a != b for a, b in izip(neA, neB))
		return distance

	def receptor_binding_distance(self, aaA, aaB):
		"""Return distance of sequences aaA and aaB by comparing receptor binding sites"""
		neA = self.receptor_binding_sites(aaA)
		neB = self.receptor_binding_sites(aaB)
		distance = sum(a != b for a, b in izip(neA, neB))
		return distance

	def add_H3N2_attributes(self):
		root = self.tree.seed_node
		root_total_aa_seq = self.get_total_peptide(root)
		for node in self.tree.postorder_node_iter():
			total_aa_seq = self.get_total_peptide(node)
			node.ep = self.epitope_distance(total_aa_seq, root_total_aa_seq)
			node.ne = self.nonepitope_distance(total_aa_seq, root_total_aa_seq)
			node.rb = self.receptor_binding_distance(total_aa_seq, root_total_aa_seq)


class H3N2_HI(HI_tree):
	def __init__(self, **kwargs):
		HI_tree.__init__(self, **kwargs)

class H3N2_fitness(fitness_model):
	def __init__(self, **kwargs):
		if 'predictors' in self.kwargs:
			predictor_input = self.kwargs['predictors']
			fitness_model.__init__(self, predictor_input = predictor_input, **kwargs)
		else:
			fitness_model.__init__(self, **kwargs)

	def annotate_fitness(self, estimate_frequencies = True):
		self.predict(estimate_frequencies=estimate_frequencies)


class H3N2_process(process, H3N2_filter, H3N2_clean, H3N2_refine, H3N2_HI, H3N2_fitness):
	"""docstring for H3N2_process, H3N2_filter"""
	def __init__(self,verbose = 0, force_include = None,
				force_include_all = False, max_global= True, **kwargs):
		self.force_include = force_include
		self.force_include_all = force_include_all
		self.max_global = max_global
		process.__init__(self, **kwargs)
		H3N2_filter.__init__(self,**kwargs)
		H3N2_clean.__init__(self,**kwargs)
		H3N2_refine.__init__(self,**kwargs)
		H3N2_HI.__init__(self,**kwargs)
		H3N2_fitness.__init__(self,**kwargs)
		self.verbose = verbose

	def run(self, steps, viruses_per_month=50, raxml_time_limit = 1.0, lam_HI=2.0, lam_pot=0.3, lam_avi=2.0):
		if 'filter' in steps:
			print "--- Virus filtering at " + time.strftime("%H:%M:%S") + " ---"
			self.filter()
			if self.force_include is not None:
				with open(self.force_include) as infile:
					forced_strains = [fix_name(line.strip().split('\t')[0]).upper() for line in infile]
			else:
				forced_strains = []
			self.subsample(viruses_per_month,
				prioritize=forced_strains, all_priority=self.force_include_all,
				region_specific = self.max_global)
			self.add_older_vaccine_viruses(dt = 3)
			self.dump()
		else:
			self.load()
		if 'align' in steps:
			self.align()   	# -> self.viruses is an alignment object
		if 'clean' in steps:
			print "--- Clean at " + time.strftime("%H:%M:%S") + " ---"
			self.clean()   # -> every node as a numerical date
			self.dump()
		if 'tree' in steps:
			print "--- Tree	 infer at " + time.strftime("%H:%M:%S") + " ---"
			self.infer_tree(raxml_time_limit)  # -> self has a tree
			self.dump()
		if 'ancestral' in steps:
			print "--- Infer ancestral sequences " + time.strftime("%H:%M:%S") + " ---"
			self.infer_ancestral()  # -> every node has a sequence
			self.dump()
		if 'refine' in steps:
			print "--- Tree refine at " + time.strftime("%H:%M:%S") + " ---"
			self.refine()
			self.dump()
		if 'frequencies' in steps:
			print "--- Estimating frequencies at " + time.strftime("%H:%M:%S") + " ---"
			self.determine_variable_positions()
			self.estimate_frequencies(tasks = ["mutations", "tree"])
			if 'genotype_frequencies' in steps:
					self.estimate_frequencies(tasks = ["genotypes"])
			self.dump()

		method = 'nnl1reg'
		if 'HI' in steps:
			print "--- Adding HI titers to the tree " + time.strftime("%H:%M:%S") + " ---"
			try:
				self.determine_variable_positions()
				self.map_HI(training_fraction=1.0, method = 'nnl1reg',
					lam_HI=lam_HI, lam_avi=lam_avi, lam_pot=lam_pot, map_to_tree=True)
				self.map_HI(training_fraction=1.0, method = 'nnl1reg', force_redo=True,
					lam_HI=lam_HI, lam_avi=lam_avi, lam_pot=lam_pot, map_to_tree=False)
			except:
				print("HI modeling failed!")
			#freqs = self.determine_HI_mutation_frequencies(threshold = 0.1)
			#self.frequencies["mutations"]["global"].update(freqs)
			self.dump()

		if 'fitness' in steps:
			print "--- Estimating fitnesses at " + time.strftime("%H:%M:%S") + " ---"
			self.annotate_fitness()
			self.dump()

		if 'export' in steps:
			self.add_titers()
			self.temporal_regional_statistics()
			# exporting to json, including the H3N2 specific fields
			self.export_to_auspice(tree_fields = [
				'ep', 'ne', 'rb', 'aa_muts','accession','isolate_id', 'lab', 'db', 'country', 'dfreq', 'fitness', 'pred_distance',
				'dHI', 'cHI', 'mHI', 'mean_HI_titers', 'HI_titers', 'HI_titers_raw', 'serum', 'HI_info',
				'avidity_tree', 'avidity_mut', 'potency_mut', 'potency_tree', 'mean_potency_mut', 'mean_potency_tree', 'autologous_titers'],
				   annotations = ['3c2.a', '3c3.a', '3c3.b', '171K'])
			if params.html:
				self.generate_indexHTML()
			self.export_HI_mutation_effects()
			#self.export_clade_frequencies()
			#self.export_viruses()

		if 'HIvalidate' in steps:
			from diagnostic_figures import tree_additivity_symmetry, fmts

			print "--- generating validation figures " + time.strftime("%H:%M:%S") + " ---"
			print "-- number of non-zero branch parameters: ",np.sum([n.dHI>1e-3 for n in self.tree.postorder_node_iter()])
			print "-- number of non-zero mutation parameters: ",np.sum([val>1e-3 for val in self.mutation_effects.values()])
			for model in ['tree', 'mutation']:
				try:
					tree_additivity_symmetry(self, model)
					for fmt in fmts: plt.savefig(self.htmlpath()+'HI_symmetry_'+model+fmt)
				except:
					print("Can't generate symmetry/additivity figures")
			try:
				self.slopes_muts = slope_vs_mutation(self)
			except:
				print("Couldn't derive slopes, probably to small time interval")
			self.generate_validation_figures(method)


if __name__=="__main__":
	all_steps = ['filter', 'align', 'clean', 'tree', 'ancestral', 'refine',
				 'frequencies', 'HI', 'fitness', 'export', 'HIvalidate']

	from process import parser
	import matplotlib.pyplot as plt
	plt.ion()
	params = parser.parse_args()

	lt = time.localtime()
	num_date = round(lt.tm_year+(lt.tm_yday-1.0)/365.0,2)
	params.time_interval = (num_date-params.years_back, num_date)
	if params.interval is not None and len(params.interval)==2 and params.interval[0]<params.interval[1]:
		params.time_interval = (params.interval[0], params.interval[1])
	dt= params.time_interval[1]-params.time_interval[0]
	params.pivots_per_year = 12.0 if dt<5 else 6.0
	steps = all_steps[all_steps.index(params.start):(all_steps.index(params.stop)+1)]
	if params.skip is not None:					# params.skip will be a string ("genotype_frequencies HIvalidate") if called from make_all, and a list (["genotype_frequencies", "HIvalidate"]) if called directly from process
		if type(params.skip) is str:
			params.skip = params.skip.split()	# params.skip is definitely a list
		for tmp_step in params.skip:
			if tmp_step in steps:
				print "skipping", tmp_step
				steps.remove(tmp_step)

	# add all arguments to virus_config (possibly overriding)
	virus_config.update(params.__dict__)
	virus_config['serum_Kc'] = 0.003
	# pass all these arguments to the processor: will be passed down as kwargs through all classes
	myH3N2 = H3N2_process(**virus_config)
	if params.test:
		myH3N2.load()
	else:
		myH3N2.run(steps, viruses_per_month = virus_config['viruses_per_month'],
				   raxml_time_limit = virus_config['raxml_time_limit'],
				   lam_HI = virus_config['lam_HI'],
				   lam_avi = virus_config['lam_avi'],
				   lam_pot = virus_config['lam_pot'],
				   )
