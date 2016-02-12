from __future__ import print_function
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import json
import cPickle as pickle
from collections import defaultdict

tree = pickle.load(open('../../augur/data/H3N2_1y_tree.pkl', 'r'))

#mutation = (160, 'T','K', 'R')
mutation = (53, 'D','N')
conditional = [(159, 'Y')]

count_passage = {mut:defaultdict(int) for mut in mutation[1:]}
for leaf in tree.leaf_iter():
	if hasattr(leaf,'passage'):
		if all([leaf.aa_seq['HA1'][pos-1]==aa for pos, aa in conditional]):
			mut = leaf.aa_seq['HA1'][mutation[0]-1]
			if mut in count_passage:
				if 'MDCK' in leaf.passage.upper():
					passage_info = 'MDCK'
				elif 'SIAT' in leaf.passage.upper():
					passage_info = 'SIAT'
				elif any([x in leaf.passage for x in ['linical', 'riginal','irect']]):
					passage_info = 'direct'
				else:
					passage_info=leaf.passage

				count_passage[mut][passage_info] += 1
	else:
		print(leaf.strain, 'without passage info')

for mut in mutation[1:]:
	print('########################')
	print('# ', mut)
	print('########################')
	tmp = sorted(count_passage[mut].items(), key=lambda x:-x[1])
	for cat, count in tmp:
		print(cat, count)
