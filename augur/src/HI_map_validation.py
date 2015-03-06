from H3N2_process import *
import matplotlib.pyplot as plt
from tree_titer import plot_tree, plot_dHI_distribution
import cPickle


if __name__=="__main__":
	parser = argparse.ArgumentParser(description='Process virus sequences, build tree, and prepare of web visualization')
	parser.add_argument('--prefix', type = str, default = 'data/', help='path+prefix of file dumps')
	parser.add_argument('--training', type = float, default = 0.8, help='fraction of data used for training')
	parser.add_argument('--reg', type = float, default = 1.0, help='regularization parameter')
	parser.add_argument('--avi', type = float, default = 1.0, help='regularization parameter')
	parser.add_argument('--pot', type = float, default = 1.0, help='regularization parameter')

	params = parser.parse_args()
	virus_config.update(params.__dict__)
	# pass all these arguments to the processor: will be passed down as kwargs through all classes
	myH3N2 = H3N2_process(**virus_config) 
	myH3N2.load()
	fig_prefix = 'figures/'+params.prefix.split('/')[-1]

	####  FIT VALIDATION  #######################################################
	myH3N2.map_HI_to_tree(training_fraction=params.training, method = 'nnl1reg', 
		lam_HI=params.reg, lam_pot=params.pot, lam_avi=params.avi)
	myH3N2.validate(plot=True)
	plt.savefig(fig_prefix+'HI_scatter.pdf')

	####  effects and mutations  #######################################################
	plot_dHI_distribution(myH3N2.tree)
	plt.savefig(fig_prefix+'dHI_distribution.pdf')

	####  cHI colored tree  #######################################################
	plot_tree(myH3N2.tree)
	plt.savefig(fig_prefix+'cHI_tree.pdf')

	####  VIRUS EFFECTS   #######################################################
	plt.figure()
	plt.title('histogram of inferred virus effects')
	plt.hist(myH3N2.virus_effect.values())
	plt.xlabel('virus avidities')
	plt.savefig(fig_prefix+'HI_avidity_histogram.pdf')

	####  SERUM POTENCIES  #######################################################
	with open(fig_prefix+'HI_potencies.txt','w') as outfile:
		for serum, val in myH3N2.serum_potency.iteritems():
			outfile.write(serum+'\t'+str(round(val,4))+'\n')


	####  DISTANCE ASYMMETRIES #######################################################
	reciprocal_measurements = []
	for (testvir, serum) in myH3N2.HI_normalized:
		if (serum, testvir) in myH3N2.HI_normalized:
			val_fwd = myH3N2.HI_normalized[(testvir,serum)]
			val_bwd = myH3N2.HI_normalized[(serum, testvir)]
			diff_uncorrected = val_fwd - val_bwd
			diff_corrected = (val_fwd - myH3N2.serum_potency[serum] - myH3N2.virus_effect[testvir])\
							-(val_bwd - myH3N2.serum_potency[testvir] - myH3N2.virus_effect[serum])
			reciprocal_measurements.append([testvir, serum, diff_uncorrected, diff_corrected])

	plt.figure()
	plt.title('asymmetry in reciprocal titers')
	plt.hist([x[2] for x in reciprocal_measurements],alpha=0.7, label="uncorrected", normed=True)
	plt.hist([x[3] for x in reciprocal_measurements],alpha=0.7, label="corrected", normed=True)
	plt.xlabel('distance asymmetry')
	plt.legend()
	plt.savefig(fig_prefix+'HI_titer_asymmetry.pdf')
