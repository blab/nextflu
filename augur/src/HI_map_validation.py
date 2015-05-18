import matplotlib.pyplot as plt
from tree_titer import plot_tree, plot_dHI_distribution
import cPickle, argparse

def validation_figures(params):
	virus_config.update(params.__dict__)
	# pass all these arguments to the processor: will be passed down as kwargs through all classes
	myflu = flu_process(**virus_config) 
	myflu.load()
	fig_prefix = 'figures/'+params.prefix.split('/')[-1]

	####  FIT VALIDATION  #######################################################
	myflu.map_HI_to_tree(training_fraction=params.training, method = 'nnl1reg', 
		lam_HI=params.reg, lam_pot=params.pot, lam_avi=params.avi, subset_strains = params.train_strains)
	myflu.validate(plot=True)
	plt.savefig(fig_prefix+'HI_scatter.pdf')

	####  effects and mutations  #######################################################
	plot_dHI_distribution(myflu.tree)
	plt.savefig(fig_prefix+'dHI_distribution.pdf')

	####  cHI colored tree  #######################################################
	plot_tree(myflu.tree)
	plt.savefig(fig_prefix+'cHI_tree.pdf')

	####  VIRUS EFFECTS   #######################################################
	plt.figure()
	plt.title('histogram of inferred virus effects')
	plt.hist(myflu.virus_effect.values())
	plt.xlabel('virus avidities')
	plt.savefig(fig_prefix+'HI_avidity_histogram.pdf')
	print "virus effects:", np.mean(myflu.virus_effect.values()), np.std(myflu.virus_effect.values())

	####  SERUM POTENCIES  #######################################################
	with open(fig_prefix+'HI_potencies.txt','w') as outfile:
		for serum, val in myflu.serum_potency.iteritems():
			outfile.write(serum[0]+'\t'+serum[1]+'\t'+str(round(val,4))+'\n')
	print "potencies:", np.mean(myflu.serum_potency.values()), np.std(myflu.serum_potency.values())

	####  DISTANCE ASYMMETRIES #######################################################
	reciprocal_measurements = []
	reciprocal_measurements_titers = []
	for (testvir, serum) in myflu.HI_normalized:
		tmp_recip = [v for v in myflu.HI_normalized if serum[0]==v[0] and testvir==v[1][0]]
		for v in tmp_recip:
			val_fwd = myflu.HI_normalized[(testvir,serum)]
			val_bwd = myflu.HI_normalized[v]
			diff_uncorrected = val_fwd - val_bwd
			diff_corrected = (val_fwd - myflu.serum_potency[serum] - myflu.virus_effect[testvir])\
							-(val_bwd - myflu.serum_potency[v[1]] - myflu.virus_effect[serum[0]])
			val_bwd = myflu.HI_normalized[v]
			reciprocal_measurements.append([testvir, serum, diff_uncorrected, diff_corrected])
			reciprocal_measurements_titers.append([testvir, serum, val_fwd, val_bwd, 
			                                      (val_fwd - myflu.serum_potency[serum] - myflu.virus_effect[testvir]),
                      							  (val_bwd - myflu.serum_potency[v[1]] - myflu.virus_effect[serum[0]]),
												  ])

	plt.figure()
	plt.title('asymmetry in reciprocal titers')
	plt.hist([x[2] for x in reciprocal_measurements],alpha=0.7, label="uncorrected", normed=True)
	plt.hist([x[3] for x in reciprocal_measurements],alpha=0.7, label="corrected", normed=True)
	plt.xlabel('distance asymmetry')
	plt.legend()
	plt.savefig(fig_prefix+'HI_titer_asymmetry.pdf')

	####  Ultrametricity #######################################################
	symmetrized = {(v,s[0]): (val_fwd, val_bwd, cval_fwd, cval_bwd) for v,s,val_fwd, val_bwd,cval_fwd, cval_bwd in reciprocal_measurements_titers}
	all_reciprocal = set([v[0] for v in reciprocal_measurements_titers])
	ultra_deviation = [[],[]]
	ultra_norm = [[],[]]
	from random import sample
	from itertools import product
	for trial in range(10000):
		four = sample(all_reciprocal, 4)
		distances = {}
		for i1, v1 in enumerate(four):
			for v2 in four:
				if (v1,v2) in symmetrized:
					if v1 != v2:
						distances[(v1,v2)] = symmetrized[(v1,v2)]
				else:
					distances[(v1,v2)] = np.nan
		if np.nan in distances.values():
			continue
		else:
			for ci in [0,2]:
				for d12 in distances[(four[0], four[1])][ci:(ci+2)]:
					for d13 in distances[(four[0], four[2])][ci:(ci+2)]:
						for d14 in distances[(four[0], four[3])][ci:(ci+2)]:
							for d23 in distances[(four[1], four[2])][ci:(ci+2)]:
								for d24 in distances[(four[1], four[3])][ci:(ci+2)]:
									for d34 in distances[(four[2], four[3])][ci:(ci+2)]:
										tmp = sorted([d12 + d34, d13 + d24, d14 + d23])
										ultra_deviation[ci/2].append(tmp[-1]-tmp[-2])
										ultra_norm[ci/2].append(tmp[-1]-tmp[0])
										#print tmp
	plt.figure()
	plt.title('deviations from ultra metricity')
	plt.hist(np.array(ultra_deviation[0])/np.mean(ultra_norm[0]),label = "uncorrected", alpha=0.7,normed=True)
	plt.hist(np.array(ultra_deviation[1])/np.mean(ultra_norm[1]),label = "corrected", alpha=0.7,normed=True)
	plt.xlabel('deviation')
	plt.legend()
	plt.savefig(fig_prefix+'HI_titer_ultrametricity.pdf')

	#### titer effects ###############################################################
	dHI_list = []
	for node in myflu.tree.postorder_node_iter():
		dHI_list.append((node.dHI, node.mutations, node))
	dHI_list.sort()
	return dHI_list

def scan_regularization(params, grid):
	virus_config.update(params.__dict__)
	# pass all these arguments to the processor: will be passed down as kwargs through all classes
	myH3N2 = flu_process(**virus_config) 
	myflu.load()
	fig_prefix = 'figures/'+params.prefix.split('/')[-1]

	####  looping over different combinations of regularizers  ########################
	accuracy = np.zeros((len(grid), len(grid), len(grid), 5))
	for hi, lam_HI in enumerate(grid):
		for pi, lam_pot in enumerate(grid):
			for ai, lam_avi in enumerate(grid):
				myflu.map_HI_to_tree(training_fraction=params.training, method = 'nnl1reg', 
				lam_HI=lam_HI, lam_pot=lam_pot, lam_avi=lam_avi, subset_strains = params.train_strains)
				myflu.validate(plot=False)
				####  calculated asymmetries 
				reciprocal_measurements = []
				for (testvir, serum) in myflu.HI_normalized:
					tmp_recip = [v for v in myflu.HI_normalized if serum[0]==v[0] and testvir==v[1][0]]
					for v in tmp_recip:
						val_fwd = myflu.HI_normalized[(testvir,serum)]
						val_bwd = myflu.HI_normalized[v]
						diff_uncorrected = val_fwd - val_bwd
						diff_corrected = (val_fwd - myflu.serum_potency[serum] - myflu.virus_effect[testvir])\
										-(val_bwd - myflu.serum_potency[v[1]] - myflu.virus_effect[serum[0]])
						reciprocal_measurements.append([testvir, serum, diff_uncorrected, diff_corrected])
				accuracy[hi,pi, ai]=myflu.rms_error, myflu.abs_error, myflu.slope, myflu.intercept, np.std([x[3] for x in reciprocal_measurements])
				print lam_HI, lam_pot, lam_avi, accuracy[hi,pi,ai]

	return accuracy

if __name__=="__main__":
	plt.ion()
	parser = argparse.ArgumentParser(description='Process virus sequences, build tree, and prepare of web visualization')
	parser.add_argument('--prefix', type = str, default = 'data/', help='path+prefix of file dumps')
	parser.add_argument('--flutype', type = str, default = 'H3N2', help='flu strain')
	parser.add_argument('--training', type = float, default = 0.8, help='fraction of data used for training')
	parser.add_argument('--train_strains', default = False, action = 'store_true', help='subset measurements or strains to train')
	parser.add_argument('--reg', type = float, default = 1.0, help='regularization parameter')
	parser.add_argument('--avi', type = float, default = 1.0, help='regularization parameter')
	parser.add_argument('--pot', type = float, default = 1.0, help='regularization parameter')
	parser.add_argument('--resolution', type = str,  help ="label for the resolution")	
	parser.add_argument('--min_aamuts', type = str, default = '0', help='minimal number of aminoacid mutations to include branch or epi for epitope or rbs for receptor binding site')
	params = parser.parse_args()
	if params.flutype=='H3N2':
		from H3N2_process import *
		flu_process = H3N2_process
	elif params.flutype=='H1N1pdm':
		from H1N1pdm_process import *
		flu_process = H1N1pdm_process
	elif params.flutype=='H1N1':
		from H1N1_process import *
		flu_process = H1N1_process
	elif params.flutype=='Vic':
		from Vic_process import *
		flu_process = BVic_process
	elif params.flutype=='Yam':
		from Yam_process import *
		flu_process = BYam_process
	try:
		params.min_aamuts = int(params.min_aamuts)
	except:
		pass
	params.__dict__['HI_fname']='source-data/'+params.flutype+'_HI_titers.txt'	

	dHI_list = validation_figures(params)
	#grid = [0.1, 0.3, 1, 3, 10]
	#accuracy = scan_regularization(params, grid)
#
#	#for ii in range(accuracy.shape[-1]):
#	#	plt.figure()
#	#	for hi, lam_HI in enumerate(grid):
#	#		plt.subplot(len(grid)//2+len(grid)%2, 2, hi+1)
#	#		plt.title('lam_HI='+str(lam_HI))
	#		plt.imshow(accuracy[hi,:,:,ii], interpolation='nearest')

