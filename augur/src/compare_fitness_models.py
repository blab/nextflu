from H3N2_process import *
import matplotlib.pyplot as plt
import cPickle

if __name__=="__main__":
	parser = argparse.ArgumentParser(description='Process virus sequences, build tree, and prepare of web visualization')
	parser.add_argument('--prefix', type = str, default = 'data/', help='path+prefix of file dumps')
	parser.add_argument('--seasons', default = [1995, 2015], help='year range to evalutate predictors for', nargs='+')
	params = parser.parse_args()
	virus_config.update(params.__dict__)
	virus_config["time_interval"] = map(float, params.seasons) 
	# pass all these arguments to the processor: will be passed down as kwargs through all classes
	myH3N2 = H3N2_process(**virus_config) 
	myH3N2.load()

	fig_prefix = 'figures/'+params.prefix.split('/')[-1]
	### EVALUTE INDIVIDUAL PREDICTORS  ##########################################
	avg_distance = None
	pred_errors = {}
	pred_suscept = {}
	pred_params = {}
	plt.figure()
	plt.plot(virus_config["time_interval"], [1,1], c='k')
	for pred_set in [['lb'], ['ep'], ['tol'], ['ne_star']]: #, ['HI']]:
		myH3N2.predictors = pred_set
		myH3N2.predict(niter=0)
		if avg_distance is None:
			avg_distance = np.array([myH3N2.allele_frequency_distance(myH3N2.season_af[s], myH3N2.season_af[t])
							for s,t in myH3N2.fit_test_season_pairs])
		pred_label = "/".join(pred_set)
		pred_errors[pred_label] = np.array(myH3N2.seasonal_errors)
		pred_suscept[pred_label] = np.array(myH3N2.seasonal_susceptibility)
		pred_params[pred_label] = np.array(myH3N2.params)
		plt.plot([x[0][0] for x in myH3N2.fit_test_season_pairs], pred_errors[pred_label]/avg_distance, 
				label = pred_label+":"+", ".join(map(str,map(lambda x:round(x,3), pred_params[pred_label]))))
		plt.plot([x[0][0] for x in myH3N2.fit_test_season_pairs], pred_suscept[pred_label]/avg_distance)
	plt.legend()
	plt.xlabel('year')
	plt.ylabel('distance to next (rel. avg)')
	plt.savefig(fig_prefix+'single_prediction.pdf')

	### EVALUTE PAIRS OF PREDICTORS  ##########################################
	plt.figure()
	plt.plot(virus_config["time_interval"], [1,1], c='k')
	for pred_set in [['lb', 'tol']]: #, ['lb', 'HI'], ['tol', 'HI']]:
		myH3N2.predictors = pred_set
		myH3N2.predict(niter=0)
		if avg_distance is None:
			avg_distance = np.array([myH3N2.allele_frequency_distance(myH3N2.season_af[s], myH3N2.season_af[t])
							for s,t in myH3N2.fit_test_season_pairs])
		pred_label = "/".join(pred_set)
		pred_errors[pred_label] = np.array(myH3N2.seasonal_errors)
		pred_params[pred_label] = np.array(myH3N2.params)
		plt.plot([x[0][0] for x in myH3N2.fit_test_season_pairs], pred_errors[pred_label]/avg_distance, 
				label = pred_label+":"+", ".join(map(str,map(lambda x:round(x,3), pred_params[pred_label]))))

	plt.legend()
	plt.xlabel('year')
	plt.ylabel('distance to next (rel. avg)')
	plt.savefig(fig_prefix+'combination_prediction.pdf')

	with open(fig_prefix+'fit_results.pkl', 'w') as pfile:
		cPickle.dump((pred_params, pred_errors), pfile)

