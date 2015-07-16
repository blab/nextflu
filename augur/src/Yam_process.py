import time, re, os
from virus_filter import flu_filter, fix_name
from virus_clean import virus_clean
from tree_refine import tree_refine
from tree_titer import HI_tree
from H3N2_process import H3N2_refine as BYam_refine
from process import process, virus_config
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Align import MultipleSeqAlignment
import numpy as np
from itertools import izip

sp = 15
epitope_mask = np.array(['1' if pos in [141,142,145,146,172,176,178,179,180,181,183,184,185, #Sa
										170,173,174,177,206,207,210,211,212,214,216,		 #Sb
										183,187,191,196,221,225,254,258,288,				 #Ca1
										154,157,158,159,161,163,238,239,242,243,			 #Ca2
										87, 88, 90, 91, 92, 95, 96, 98, 99, 100, 132, 139	 #Cb
									   ]
						else '0' for pos in xrange(1,1725)])

receptor_binding_sites = [159,169,170,172,173,203,207]


virus_config.update({
	# data source and sequence parsing/cleaning/processing
	'virus':'Yam',
	'alignment_file':'data/Yam_gisaid_epiflu_sequence.fasta.gz',
	'outgroup':'B/Singapore/11/94',
	'force_include':'source-data/Yam_HI_strains.txt',
	'force_include_all':True,
	'date_spec':'year',
	'max_global':True,   # sample as evenly as possible from different geographic regions 
	'cds':[11,None], # define the translation start i n 0 numbering
	# define relevant clades in canonical HA1 numbering (+1)
	# numbering starting at methionine including the signal peptide
	'clade_designations': {
		'2':  [('HA1', 48,'K'), ('HA1', 108, 'A'), ('HA1', 150, 'S')],
		'3':  [('HA1', 48,'R'), ('HA1', 108, 'P'), ('HA1', 150, 'I')],
		'3a': [('HA1', 37,'A'), ('HA1', 298, 'E'), ('HA1', 48,'R'), ('HA1', 105, 'P'), ('HA1', 150, 'I')],
	},
	'auspice_prefix':'Yam_',
	'HI_fname':'source-data/Yam_HI_titers.txt',
	'html_vars': {'coloring': 'lbi, dfreq, region, date, cHI, HI_dist',
				  'gtplaceholder': 'HA1 positions...',
				  'freqdefault': '2, 3, 3a'},
	'js_vars': {'LBItau': 0.0005, 'LBItime_window': 0.5, 'dfreq_dn':2},	
	})


class BYam_filter(flu_filter):
	def __init__(self,min_length = 987, **kwargs):
		'''
		parameters
		min_length  -- minimal length for a sequence to be acceptable
		'''
		flu_filter.__init__(self, **kwargs)
		self.min_length = min_length
		self.vaccine_strains =[
			{
				'strain':    	'B/Beijing/184/93',
				'isolate_id':	'EPI_ISL_969',
				'date':    		'1993-07-01', #(Month and day unknown)
				'region':   	'China', 
				'seq':     		'GATCGAATCTGTACTGGGATAACATCTTCAAACTCACCTCATGTGGTCAAAACAGCTACTCAAGGGGAGGTCAATGTGACTGGTGTGATACCACTGACAACAACACCAACAAAATCTCATTTTGGAAATCTCAAAGGAACAAAGACCAGAGGGAAACTATGCCCAAACTGTCTCAACTGCACAGATCTGGATGTGGCCTTGGGCAGACCAATGTGTGTGGGGACCACACCTTCGGCAAAAGCTTCAATACTCCACGAAGTCAGACCTGTTACATCCGGGTGCTTTCCTATAATGCACGACAGAACAAAAATCAGACAGCTACCCAATCTTCTCAGAGGATATGAAAATATCAGATTATCAACCCAAAACGTTATCAACGCAGAAAAGGCACCAGGAGGACCCTACAGGCTTGGAACCTCAGGATCTTGCCCTAACGCTACCAGTAGAAGCGGATTTTTCGCAACAATGGCTTGGGCTGTCCCAAGGGACAACAACAAAACAGCAACAAATCCACTAACAGTAGAAGTACCATACATTTGTACAAAAGGAGAAGACCAAATTACTGTTTGGGGGTTCCATTCTGATAACAAAATCCAAATGAAAAACCTCTATGGAGACTCAAATCCTCAAAAGTTCACCTCATCTGCTAATGGAGTAACCACACATTATGTTTCTCAGATTGGCGGCTTCCCAGATCAAACAGAAGACGGAGGGCTACCACAAAGCGGCAGAATTGTTGTTGATTACATGGTGCAAAAACCTGGGAAAACAGGAACAATTGTCTATCAAAGAGGTGTTTTGTTGCCTCAAAAGGTGTGGTGCGCGAGTGGCAGGAGCAAAGTAATAAAAGGGTCCTTGCCTTTAATTGGTGAAGCAGATTGCCTTCACGAAAAATACGGTGGATTAAACAAAAGCAAGCCTTACTACACAGGAGAACATGCAAAAGCCATAGGAAATTGCCCAATATGGGTGAAAACACCTTTGAAGCTTGCCAATGGAACCAAATATAGACCTCCTGCAAAACTATTAAAGGAAAGG',
			},
			{	
				'strain':    	'B/Sichuan/379/99',
				'isolate_id': 	'EPI_ISL_21113',
				'date':    		'1999-07-01', # (Month and day unknown)	
				'region':   	'China',
				'seq':     		'GAGGCAATAATTGTACTACTCATGGTAGTAACATCCAATGCAGATCGAATCTGCACTGGGATAACATCGTCAAACTCACCTCATGTGGTCAAAACAGCTACTCAAGGGGAGGTCAATGTAACTGGTGCGATACCACTGACAACAACACCAACAAAATCTCATTTTGCAAATCTCAAAGGAACAAAGACCAGAGGGAAACTATGCCCAACCTGTCTCAACTGCACAGATCTGGATGTGGCCTTGGGCAGACCAATGTGTGTGGGGATCACACCTTCGGCAAAAGCTTCAATACTCCACGAAATCAAACCTGTTACATCCGGATGCTTTCCTATAATGCACGACAGAACAAAAATCAGACAGCTACCCAATCTTCTCAGAGGATATGAAAAAATCAGATTATCAACCCAAAACGTTATCAACGCAGAAAAGGCACCAGGAGGACCTTACAGACTTGGAACTTCAGGATCTTGCCCTAACGCTACCAGTAAAAGCGGATTTTTCGCAACAATGGCTTGGGCTGTCCCAAGGGACAACAACAAAACAGCAACGAATCCACTAACAGTAGAAGTACCACACATCTGTACAAAAGAAGAAGACCAAATTACTGTTTGGGGGTTCCATTCTGATGACAAAACCCAAATGAAAAACCTCTATGGAGACTCAAATCCTCAAAAGTTCACCTCATCTGCTAATGGAATAACCACACATTATGTTTCTCAGATTGGCGGCTTCCCGGACCAAACAGAGGACGGAGGGCTACCACAAAGCGGCAGAATTGTTGTTGATTACATGGTGCAAAAACCTGGGAAAACAGGAACAATTGTCTATCAAAGAGGGATTTTGTTGCCTCAAAAGGTGTGGTGCGCGAGTGGTAGGAGCAAAGTAATAAAAGGGTCCTTGCCTTTAATTGGTGAAGCAGATTGCCTTCATGAAAAATACGGTGGATTAAACAAAAGCAAGCCTTACTACACAGGAGAACATGCAAAAGCCATAGGAAATTGCCCAATATGGGTGAAAACACCTTTGAAGCTTGCCAATGGAACCAAATATAGACCTCCTGCAAAACTATTAAAGGAAAGGGGTTTCT',
			},
			{
				'strain':    	'B/Shanghai/361/2002',
				'isolate_id': 	'EPI_ISL_2842',
				'date':    		'2002-06-12',
				'region':   	'China',
				'seq': 			'AATGCAGATCGAATCTGCACTGGGATAACATCTTCAAACTCACCTCATGTGGTCAAAACAGCTACTCAAGGGGAGGTCAATGTGACTGGTGTGATACCACTGACAACAACTCCAATAAAATCTCATTTTGCAAATCTCAAAGGAACAAGGACTAGAGGGAAACTATGCCCAGATTGTCTCAACTGCACAGATCTGGATGTGGCCTTGGGCAGACCAATGTGTGTGGGGACCACACCTTCGGCAAAAGCTTCAATACTCCACGAAGTCAGACCTGTTACATCCGGGTGCTTTCCTATAATGCACGACAGAACAAAAATCAGACAACTACCCAATCTTCTCAGAGGATATGAAAATATCAGGTTATCAACCCAAAACGTTATCGATGCAGAAAAGGCCCTAGGAGGACCCTACAGACTTGGAACCTCAGGATCTTGCCCTAACGCCACCAGTAAAAGCGGATTTTTCGCAACAATGGCTTGGGCTGTCCCAAAGGACAACAACAAAAATGCAACGAACCCACTAACAGTAGAAGTACCATACATCTGTACAGAAGGGGAAGACCAAATTACTGTTTGGGGGTTCCATTCAGATGACAAAACCCAAATGAAAAACCTCTATGGAGACTCAAATCCTCAAAAGTTCACCTCATCTGCTAATGGAGTAACCACACATTATGTTTCTCAGATTGGCGGCTTCCCAGATCAAACAGAAGACGGAGGACTACCACAAAGCGGCAGAATTGTTGTTGATTACATGGTGCAAAAACCTGGGAAAACAGGAACAATTGTCTATCAAAGAGGTGTTTTGTTGCCTCAAAAGGTGTGGTGCGCGAGTGGCAGGAGCAAAGTAATAAAAGGGTCCTTGCCTTTAATTGGTGAAGCAGATTGCCTTCATGAAAAATACGGTGGGTTAAACAAAAGCAAGCCTTACTACACAGGAGAACATGCAAAAGCCATAGGAAATTGCCCAATATGGGTGAAAACACCTTTGAAGCTTGCCAATGGAACCAAATATAGACCTCCTGCAAAACTATTAAAGGAAAGGGGTTTC',
			},
			{
				'strain':		'B/Florida/4/2006',
				'isolate_id':	'EPI_ISL_21307',
				'date':			'2006-11-01',
				'region':		'NorthAmerica',
				'seq':			'ATGAAGGCAATAATTGTACTACTCATGGTAGTAACATCCAATGCAGATCGAATCTGCACTGGAATAACATCTTCAAACTCACCTCATGTGGTCAAAACAGCCACTCAAGGGGAGGTCAATGTGACTGGTGTGATACCACTAACAACAACACCAACAAAATCTTATTTTGCAAATCTCAAAGGAACAAGGACCAGAGGGAAACTATGCCCAGACTGTCTCAACTGCACAGATCTGGATGTGGCTTTGGGCAGACCAATGTGTGTGGGGACCACACCTTCGGCGAAAGCTTCAATACTCCACGAAGTCAAACCTGTTACATCCGGGTGCTTTCCTATAATGCACGACAGAACAAAAATCAGGCAACTACCCAATCTTCTCAGAGGATATGAAAATATCAGGCTATCAACCCAAAACGTCATCGATGCGGAAAAGGCACCAGGAGGACCCTACAGACTTGGAACCTCAGGATCTTGCCCTAACGCTACCAGTAAGAGCGGATTTTTCGCAACAATGGCTTGGGCTGTCCCAAAGGACAACAACAAAAATGCAACGAACCCACTAACAGTAGAAGTACCATACATTTGTACAGAAGGGGAAGACCAAATCACTGTTTGGGGGTTCCATTCAGATGACAAAACCCAAATGAAGAACCTCTATGGAGACTCAAATCCTCAAAAGTTCACCTCATCTGCTAATGGAGTAACCACACACTATGTTTCTCAGATTGGCAGCTTCCCAGATCAAACAGAAGACGGAGGACTACCACAAAGCGGCAGGATTGTTGTTGATTACATGATGCAAAAACCTGGGAAAACAGGAACAATTGTCTACCAAAGAGGTGTTTTGTTGCCTCAAAAGGTGTGGTGCGCGAGTGGCAGGAGCAAAGTAATAAAAGGGTCCTTGCCTTTAATTGGTGAAGCAGATTGCCTTCATGAAAAATACGGTGGATTAAACAAAAGCAAGCCTTACTACACAGGAGAACATGCAAAAGCCATAGGAAATTGCCCAATATGGGTGAAAACACCTTTGAAGCTTGCCAATGGAACCAAATATAGACCTCCTGCAAAACTATTAAAGGAAAGGGGTTTCTTCGGAGCTATTGCTGGTTTCCTAGAAGGAGGATGGGAAGGAATGATTGCAGGCTGGCACGGATACACATCTCACGGAGCACATGGAGTGGCAGTGGCGGCGGACCTTAAGAGTACGCAAGAAGCTATAAACAAGATAACAAAAAATCTCAATTCTTTGAGTGAGCTAGAAGTAAAGAATCTTCAAAGACTAAGTGGTGCCATGGATGAACTCCACAACGAAATACTCGAGCTGGATGAGAAAGTGGATGATCTCAGAGCTGACACTATAAGCTCGCAAATAGAACTTGCAGTCTTGCTTTCCAACGAAGGAATAATAAACAGTGAAGATGAGCATCTATTGGCACTTGAGAGAAAACTAAAGAAAATGCTGGGTCCCTCTGCTGTAGAGATAGGAAATGGATGCTTCGAAACCAAACACAAGTGCAACCAGACCTGCTTAGACAGGATAGCTGCTGGCACCTTTAATGCAGGAGAATTTTCTCTCCCCACTTTTGATTCACTGAACATTACTGCTGCATCTTTAAATGATGATGGATTGGATAACCATACTATACTGCTCTATTACTCAACTGCTGCTTCTAGTTTGGCTGTAACATTGATGCTAGCTATTTTTATTGTTTATATGGTCTCCAGAGACAACGTTTCATGCTCCATCTGTCTATAA'
			},
			{
				'strain':		'B/Wisconsin/01/2010',
				'isolate_id':	'EPI_ISL_76940',
				'date':			'2010-02-20',
				'region':		'NorthAmerica',
				'seq':			'ATGAAGGCAATAATTGTACTACTCATGGTAGTAACATCCAATGCAGATCGAATCTGCACTGGGATAACATCTTCAAACTCACCTCATGTGGTCAAAACAGCTACTCAAGGGGAGGTCAATGTGACTGGCGTGATACCACTGACAACAACACCAACAAAATCTTATTTTGCAAATCTCAAAGGAACAAGGACCAGAGGGAAACTATGCCCGGACTGTCTCAACTGTACAGATCTGGATGTGGCCTTGGGCAGGCCAATGTGTGTGGGGACCACACCTTCTGCTAAAGCTTCAATACTCCACGAGGTCAGACCTGTTACATCCGGGTGCTTTCCTATAATGCACGACAGAACAAAAATCAGGCAACTACCCAATCTTCTCAGAGGATATGAAAATATCAGGTTATCAACCCAAAACGTTATCGATGCAGAAAAAGCACCAGGAGGACCCTACAGACTTGGAACCTCAGGATCTTGCCCTAACGCTACCAGTAAAATCGGATTTTTTGCAACAATGGCTTGGGCTGTCCCAAAGGACAACTACAAAAATGCAACGAACCCACTAACAGTAGAAGTACCATACATTTGTACAGAAGGGGAAGACCAAATTACTGTTTGGGGGTTCCATTCAGATAACAAAACCCAAATGAAGAGCCTCTATGGAGACTCAAATCCTCAAAAGTTCACCTCATCTGCTAATGGAGTAACCACACATTATGTTTCTCAGATTGGCGACTTCCCAGATCAAACAGAAGACGGAGGACTACCACAAAGCGGCAGAATTGTTGTTGATTACATGATGCAAAAACCTGGGAAAACAGGAACAATTGTCTATCAAAGAGGTGTTTTGTTGCCTCAAAAGGTGTGGTGCGCGAGTGGCAGGAGCAAAGTAATAAAAGGGTCATTGCCTTTAATTGGTGAAGCAGATTGCCTTCATGAAAAATACGGTGGATTAAACAAAAGCAAGCCTTACTACACAGGAGAACATGCAAAAGCCATAGGAAATTGCCCAATATGGGTAAAAACACCTTTGAAGCTTGCCAATGGAACCAAATATAGACCTCCTGCAAAACTATTGAAGGAAAGGGGTTTCTTCGGAGCTATTGCTGGTTTCCTAGAAGGAGGATGGGAAGGAATGATTGCAGGTTGGCACGGATACACATCTCACGGAGCACATGGAGTGGCAGTGGCGGCAGACCTTAAGAGTACACAAGAAGCTATAAATAAGATAACAAAAAATCTCAATTCTTTGAGTGAGCTAGAAGTAAAGAACCTTCAAAGACTAAGTGGTGCCATGGATGAACTCCACAACGAAATACTCGAGCTGGATGAGAAAGTGGATGATCTCAGAGCTGACACTATAAGCTCACAAATAGAACTTGCAGTCTTGCTTTCCAACGAAGGAATAATAAACAGTGAAGACGAGCATCTATTGGCACTTGAGAGAAAACTAAAGAAAATGCTGGGTCCCTCTGCTGTAGACATAGGAAACGGATGCTTCGAAACCAAACACAAATGCAACCAGACCTGCTTAGACAGGATAGCTGCTGGCACCTTTAATGCAGGAGAATTTTCTCTCCCCACTTTTGATTCATTGAACATTACTGCTGCATCTTTAAATGATGATGGATTGGATAACCATACTATACTGCTCTATTACTCAACTGCTGCTTCTAGTTTGGCTGTAACATTAATGCTAGCTATTTTTATTGTTTATATGGTCTCCAGAGACAACGTTTCATGCTCCATCTGTCTATAA'
			},
			{
				'strain':		'B/Massachusetts/02/2012',
				'isolate_id':	'EPI_ISL_121434',
				'date':			'2012-03-13',
				'region':		'NorthAmerica',
				'seq':			'ATGAAGGCAATAATTGTACTACTAATGGTAGTAACATCCAATGCAGATCGAATCTGCACTGGGATAACATCTTCAAACTCACCTCATGTGGTCAAAACAGCTACTCAAGGGGAGGTCAATGTGACTGGTGTGATACCACTAACAACAACACCAACAAAATCTTATTTTGCAAATCTCAAAGGAACAAAGACCAGAGGGAAACTATGCCCAGACTGTCTCAACTGTACAGATCTGGATGTGGCCCTGGGCAGGCCAATGTGTGTGGGAACTACACCTTCTGCGAAAGCTTCAATACTTCACGAAGTCAGACCTGTTACATCCGGGTGCTTCCCTATAATGCACGACAGAACAAAAATCAGGCAACTAGCCAATCTTCTCAGAGGATATGAAAATATCAGGTTATCAACCCAAAACGTTATCGATGCAGAAAAGGCACCAGGAGGACCCTACAGACTTGGAACCTCAGGATCTTGCCCTAACGCTACCAGTAAAAGCGGATTTTTCGCAACAATGGCTTGGGCTGTCCCAAAGGACAACAACAAAAATGCAACGAACCCATTAACAGTAGAAGTACCATACATTTGTGCAGAAGGGGAAGACCAAATTACTGTTTGGGGGTTCCATTCAGATAACAAAACCCAAATGAAGAACCTCTATGGAGACTCAAATCCTCAAAAGTTCACCTCATCTGCTAATGGAGTAACCACACATTATGTTTCTCAGATTGGCGGCTTCCCAGATCAAACAGAAGACGGAGGACTACCACAAAGCGGCAGAATTGTCGTTGATTACATGATGCAAAAACCTGGGAAAACAGGAACAATTGTCTATCAAAGAGGTGTTTTGTTGCCTCAAAAGGTGTGGTGCGCGAGTGGCAGGAGCAAAGTAATAAAAGGGTCCTTGCCTTTAATTGGTGAAGCAGATTGCCTTCATGAAAAATACGGTGGATTAAACAAAAGCAAGCCTTACTACACAGGAGAACATGCAAAAGCCATAGGAAATTGCCCAATATGGGTGAAAACACCTTTGAAGCTTGCCAATGGAACCAAATATAGACCTCCTGCAAAACTATTAAAGGAAAGGGGTTTCTTCGGAGCTATTGCTGGTTTCCTAGAAGGAGGATGGGAAGGAATGATTGCAGGTTGGCACGGATACACATCTCACGGAGCACATGGAGTGGCAGTTGCTGCAGACCTTAAGAGCACACAAGAAGCTATAAACAAGATAACAAAAAATCTCAACTCTTTGAGTGAGCTAGAAGTAAAGAATCTTCAAAGGCTAAGTGGTGCCATGGATGAACTCCACAACGAAATACTCGAGCTGGATGAGAAAGTGGATGACCTCAGAGCTGACACTATAAGTTCACAAATAGAACTTGCAGTCTTGCTTTCCAACGAAGGAATAATAAACAGTGAAGACGAGCATCTATTGGCACTTGAGAGAAAACTAAAGAAAATGCTGGGTCCCTCTGCTGTAGACATAGGAAATGGATGCTTCGAAACCAAACACAAATGCAACCAGACCTGCTTAGACAGGATAGCTGCTGGCACCTTTAATGCAGGAGAGTTTTCTCTCCCCACTTTTGATTCATTGAACATTACTGCTGCATCTTTAAATGATGATGGATTGGATAACCATACTATACTGCTCTATTACTCAACTGCTGCTTCTAGTTTGGCTGTAACATTGATGCTAGCTATTTTTATTGTTTATATGGTCTCCAGAGACAACGTTTCATGCTCCATCTGTCTATAA'
			},
			{
				'strain':    	'B/PHUKET/3073/2013',
				'isolate_id':	'EPI_ISL_161843',
				'date':    		'2013-11-21',
				'region':   	'SoutheastAsia',
				'seq':			'ATGAAGGCAATAATTGTACTACTCATGGTAGTAACATCCAACGCAGATCGAATCTGCACTGGGATAACATCTTCAAACTCACCTCATGTGGTCAAAACAGCTACTCAAGGGGAGGTCAATGTGACTGGCGTGATACCACTGACAACAACACCAACAAAATCTTATTTTGCAAATCTCAAAGGAACAAGGACCAGAGGGAAACTATGCCCGGACTGTCTCAACTGTACAGATCTGGATGTGGCCTTGGGCAGGCCAATGTGTGTGGGGACCACACCTTCTGCTAAAGCTTCAATACTCCATGAGGTCAGACCTGTTACATCCGGGTGCTTTCCTATAATGCACGACAGAACAAAAATCAGGCAACTACCCAATCTTCTCAGAGGATATGAAAAGATCAGGTTATCAACCCAAAACGTTATCGATGCAGAAAAAGCACCAGGAGGACCCTACAGACTTGGAACCTCAGGATCTTGCCCTAACGCTACCAGTAAAATCGGATTTTTTGCAACAATGGCTTGGGCTGTCCCAAAGGACAACTACAAAAATGCAACGAACCCACTAACAGTGGAAGTACCATACATTTGTACAGAAGGGGAAGACCAAATTACTGTTTGGGGGTTCCATTCGGATAACAAAACCCAAATGAAGAGCCTCTATGGAGACTCAAATCCTCAAAAGTTCACCTCATCTGCTAATGGAGTAACCACGCATTATGTTTCTCAGATTGGCGACTTCCCAGATCAAACAGAAGACGGAGGACTACCACAAAGCGGCAGAATTGTTGTTGATTACATGATGCAAAAACCTGGGAAAACAGGAACAATTGTCTATCAAAGGGGTGTTTTGTTGCCTCAAAAGGTGTGGTGCGCGAGTGGCAGGAGCAAAGTAATAAAAGGGTCATTGCCTTTAATTGGTGAAGCAGATTGCCTTCATGAAGAATACGGTGGATTAAACAAAAGCAAGCCTTACTACACAGGAAAACATGCAAAAGCCATAGGAAATTGCCCAATATGGGTAAAAACACCTTTGAAGCTTGCCAATGGAACCAAATATAGACCTCCTGCAAAACTATTGAAGGAAAGGGGTTTCTTCGGAGCTATTGCTGGTTTCCTAGAAGGAGGATGGGAAGGAATGATTGCAGGTTGGCACGGATACACATCTCACGGAGCACATGGAGTGGCAGTGGCGGCAGACCTTAAGAGTACACAAGAAGCTATAAATAAGATAACAAAAAATCTCAATTCTTTGAGTGAACTAGAAGTAAAGAACCTTCAAAGACTAAGTGGTGCCATGGATGAACTCCACAACGAAATACTCGAGCTGGATGAAAAAGTGGATGATCTCAGAGCTGACACTATAAGCTCACAAATAGAACTTGCAGTCTTGCTTTCCAACGAAGGAATAATAAACAGTGAAGACGAGCATCTATTGGCACTTGAGAGAAAACTAAAGAAAATGCTGGGTCCCTCTGCTGTAGACATAGGAAACGGATGCTTCGAAACCAAACACAAATGCAACCAGACCTGCTTAGACAGGATAGCTGCTGGCACCTTTAATGCAGGAGAATTTTCTCTCCCCACTTTTGATTCATTGAACATTACTGCTGCATCTTTAAATGATGATGGATTGGATAACCATACTATACTGCTCTATTACTCAACTGCTGCTTCTAGTTTGGCTGTAACATTAATGCTAGCTATTTTTATTGTTTATATGGTCTCCAGAGACAACGTTTCATGCTCCATCTGTCTATAAAGAAGGTTAGGCCTTGTATTTTCCTTTATTGTAGTGCTTGTTTGCTTGTCATCATTACAAAGAAAC'
			}
		]
		tmp_outgroup = SeqIO.read('source-data/Yam_outgroup.gb', 'genbank')
		genome_annotation = tmp_outgroup.features
		self.cds = {x.qualifiers['gene'][0]:x for x in genome_annotation
				if 'gene' in x.qualifiers and x.type=='CDS' and 
				x.qualifiers['gene'][0] in ['SigPep', 'HA1', 'HA2']}
		self.outgroup = {
				'strain':'B/Singapore/11/94',
				'isolate_id':'EPI_ISL_20980',
				'date':'1994-05-10',
				'region':'China',
				'seq':str(tmp_outgroup.seq).upper()
			}

class BYam_clean(virus_clean):
	def __init__(self,**kwargs):
		virus_clean.__init__(self, **kwargs)

	def clean_outbreaks(self):
		"""Remove duplicate strains, where the geographic location, date of sampling and sequence are identical"""
		virus_hashes = set()
		new_viruses = []
		for v in self.viruses:
			try:
				geo = re.search(r'B/([^/]+)/', v.strain).group(1)
				if geo:
					vhash = (geo, v.date, str(v.seq))
					if vhash not in virus_hashes:
						new_viruses.append(v)
						virus_hashes.add(vhash)
			except:
				print "Error parsing geo info for", v.strain
		self.viruses = MultipleSeqAlignment(new_viruses)
		return new_viruses

	def clean_outliers(self):
		from seq_util import hamming_distance as distance
		"""Remove outlier viruses"""
		remove_viruses = []
		
		outlier_seqs = [
			"-----------ATGAAGGCCATAATTGTACTACTCATGGTAGTAACATCCAATGCAGATCGAATCTGCACTGGGATAACATCTTCAAACTCACCTCATGTGGTCAAAACAGCTACTCAAGGGGAGGTCAATGTGACTGGCGTGATACCACTGACAACAACACCAACAAAATCTTATTTTGCAAATCTCAAAGGAACAAGGACCAGAGGGAAACTATGTCCGGACTGTCTCAACTGTACAGATCTGGATGTGGCCTTGGGCAGGCCAATGTGTGTGGGGACCACACCTTCTGCTAAAGCTTCAATACTCCACGAAGTCAGACCTGTTACATCCGGGTGCTTTCCTATAATGCACGACAGAACAAAAATCAGGCAACTACCCAATCTTCTCAGAGGATATGAAAATATCAGGTTATCAACCCAAAACGTTATCGATGCAGAAAAAGCACCAGGAGGACCTTACAGACTTGGAACCTCAGGATCTTGCCCTAACGCTACCAGTAAAATCGGATTTTTCGCAACAATGGCTTGGGCTGTCCCAAAGGACAACTACAAAAATGCAACGAACCCACTAACAGTAGAAGTACCATACATTTGTGCAGAAGGGGAAGACCAAATTACTGTTTGGGGGTTCCATTCAGATAACAAAACCCAAATGAAGAACCTCTATGGAGACTCAAATCCTCAAAAGTTCACCTCATCTGCTAATGGAGTAACCACACATTATGTTTCTCAGATTGGCGACTTTCCAGATCAAACAGAAGACGGAGGACTACCACAAAGCGGCAGAATTGTTGTTGATTACATGGTGCAAAGACCTGGGAAAACAGGAACAATTGTCTATCAAAGAGGCGTTTTGTTGCCTCAAAAGGTGTGGTGCGCGAGTGGCAGGAGCAAAGTAATAAAAGGGTCATTGCCTTTAATTGGTGAAGCAGATTGCCTTCATGAAAAATACGGTGGATTAAACAAAAGCAAGCCTTACTACACAGGAGAACATGCAAAGGCCATAGGAAATTGCCCAATATGGGTGAAGACACCCTTGAAGCTGGCCAATGGAACCAAATATAGACCTCCTGCAAAACTATTAAAGGAAAGGGGTTTCTTCGGAGCTATTGCTGGTTTCTTAGAAGGAGGATGGGAAGGAATGATTGCAGGTTGGCACGGATACACGTCCCATGGGGCACATGGAGTAGCGGTGGCAGCAGACCTTAAGAGCACTCAAGAGGCCATAAACAAGATAACAAAAAATCTCAACTCTTTGAGTGAGCTGGAAGTAAAGAATCTTCAAAGACTAAGCGGTGCCATGGATGAACTCCACAACGAAATACTAGAACTAGACGAGAAAGTGGATGATCTCAGAGCTGATACAATAAGCTCACAAATAGAACTCGCAGTCCTGCTTTCCAATGAAGGAATAATAAACAGTGAAGATGAACATCTCTTGGCTCTTGAAAGAAAGCTGAAGAAAATGCTGGGCCCCTCTGCTGTAGAGATAGGGAATGGATGCTTTGAAACCAAACACAAGTGCAACCAGACCTGTCTCGACAGAATAGCTGCTGGTACCTTTGATGCAGGAGAATTTTCTCTCCCCACCTTTGATTCACTGAATATTACTGCTGCATCTTTAAATGACGATGGATTGGATAATCATACTATACTGCTTTACTACTCAACTGCTGCCTCCAGTTTGGCTGTAACACTGATGATAGCTATCTTTGTTGTTTATATGGTCTCCAGAGACAATGTTTCTTGCTCCATCTGTCTATAA--------------------------------------------------------------------------"
			]

		for outlier_seq in outlier_seqs:
			for v in self.viruses:
				dist = distance(Seq(outlier_seq), v)
				if (dist < 0.02):
					remove_viruses.append(v)
					if self.verbose>1:
						print "\tremoving", v.strain				

		self.viruses = MultipleSeqAlignment([v for v in self.viruses if v not in remove_viruses])

	def clean(self):
		self.clean_generic()
		self.clean_outbreaks()
		print "Number of viruses after outbreak filtering:",len(self.viruses)
		self.clean_outliers()
		print "Number of viruses after outlier filtering:",len(self.viruses)		

class BYam_process(process, BYam_filter, BYam_clean, BYam_refine, HI_tree):
	"""docstring for BYam_process, BYam_filter"""
	def __init__(self,verbose = 0, force_include = None, 
				force_include_all = False, max_global= True, **kwargs):
		self.force_include = force_include
		self.force_include_all = force_include_all
		self.max_global = max_global
		process.__init__(self, **kwargs)
		BYam_filter.__init__(self,**kwargs)
		BYam_clean.__init__(self,**kwargs)
		BYam_refine.__init__(self,**kwargs)
		HI_tree.__init__(self,**kwargs)
		self.verbose = verbose

	def run(self, steps, viruses_per_month=50, raxml_time_limit = 1.0, lam_HI=.5, lam_avi=1, lam_pot=.1):
		if 'filter' in steps:
			print "--- Virus filtering at " + time.strftime("%H:%M:%S") + " ---"
			self.filter()
			if self.force_include is not None and os.path.isfile(self.force_include):
				with open(self.force_include) as infile:
					forced_strains = [fix_name(line.strip()).upper() for line in infile]
			else:
				forced_strains = []
			self.subsample(viruses_per_month, 
				prioritize=forced_strains, all_priority=self.force_include_all, 
				region_specific = self.max_global)
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
		if 'HI' in steps:
			print "--- Adding HI titers to the tree " + time.strftime("%H:%M:%S") + " ---"
			self.determine_variable_positions()
			self.map_HI(training_fraction=1.0, method = 'nnl1reg', 
					lam_HI=lam_HI, lam_avi=lam_avi, lam_pot=lam_pot, map_to_tree=True)
			self.map_HI(training_fraction=1.0, method = 'nnl1reg', force_redo=True,
					lam_HI=lam_HI, lam_avi=lam_avi, lam_pot=lam_pot, map_to_tree=False)
			self.dump()
		if 'export' in steps:
			self.add_titers()
			self.temporal_regional_statistics()
			# exporting to json, including the BYam specific fields
			self.export_to_auspice(tree_fields = [
				'ep', 'ne', 'rb', 'aa_muts','accession','isolate_id', 'lab','db', 'country',
				'dHI', 'cHI', 'mean_HI_titers','HI_titers','HI_titers_raw', 'serum', 'HI_info',
				'avidity_tree','avidity_mut', 'potency_mut', 'potency_tree', 'mean_potency_mut', 'mean_potency_tree'], 
				annotations = ['2', '3', '3a'])
			self.generate_indexHTML()

		if 'HIvalidate' in steps:
			print "--- generating validation figures " + time.strftime("%H:%M:%S") + " ---"
			import matplotlib.pyplot as plt
			htmlpath = '../auspice/'
			if self.virus_type is not None: 
				htmlpath+=self.virus_type+'/'
			if self.resolution is not None: 
				htmlpath+=self.resolution+'/'

			self.check_symmetry(plot=True)
			plt.savefig(htmlpath+'HI_symmetry.png')

			self.map_HI(training_fraction=0.9, method='nnl1reg',lam_HI=lam_HI, lam_avi=lam_avi, 
						lam_pot = lam_pot, force_redo=True, map_to_tree=False)
			self.validate(plot=True)
			plt.savefig(htmlpath+'HI_prediction.png')

if __name__=="__main__":
	all_steps = ['filter', 'align', 'clean', 'tree', 'ancestral', 'refine', 'frequencies','HI', 'export']
	from process import parser
	params = parser.parse_args()

	lt = time.localtime()
	num_date = round(lt.tm_year+(lt.tm_yday-1.0)/365.0,2)
	params.time_interval = (num_date-params.years_back, num_date) 
	if params.interval is not None and len(params.interval)==2 and params.interval[0]<params.interval[1]:
		params.time_interval = (params.interval[0], params.interval[1])
	dt= params.time_interval[1]-params.time_interval[0]
	params.pivots_per_year = 12.0 if dt<5 else 6.0 if dt<10 else 3.0
	steps = all_steps[all_steps.index(params.start):(all_steps.index(params.stop)+1)] + ["HIvalidate"]
	if params.skip is not None:
		for tmp_step in params.skip:
			if tmp_step in steps:
				print "skipping",tmp_step
				steps.remove(tmp_step)

	# add all arguments to virus_config (possibly overriding)
	virus_config.update(params.__dict__)
	# pass all these arguments to the processor: will be passed down as kwargs through all classes
	myBYam = BYam_process(**virus_config) 
	if params.test:
		myBYam.load()
	else:
		myBYam.run(steps,viruses_per_month = virus_config['viruses_per_month'], 
			raxml_time_limit = virus_config['raxml_time_limit'])
