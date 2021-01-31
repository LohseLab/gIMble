import numpy as np
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use("Agg")
import scipy.stats as stats

def chisquare(observed, expected, p_th=0.05, recombination=False):
	#expected counts need to be larger than 5 for chisquare
	#combining bins with less than 5 counts into one bin.
	obs = np.reshape(observed, -1)
	exp = np.reshape(expected, -1)
	if not recombination:
		assert np.all(obs[exp == 0] == 0), "chisquare error: at least one cell with expected frequency 0 has observed values"
	#bin all values with counts smaller than 5
	binning_idxs = np.digitize(exp, np.array([5]), right=True)
	exp_smaller = np.sum(exp[binning_idxs==0])
	obs_smaller = np.sum(obs[binning_idxs==0])
	exp_binned = np.append(exp[binning_idxs==1], exp_smaller)
	obs_binned = np.append(obs[binning_idxs==1], obs_smaller)
	not_zeros = exp_binned > 0
	if sum(not_zeros) > 1:
		chisq = stats.chisquare(obs_binned[not_zeros], exp_binned[not_zeros])
		print("chisquare value:", chisq)
		assert chisq.pvalue > p_th

def bonferroni(p_threshold, pvals):
	assert min(pvals) > p_threshold / len(pvals)

def scatter_loglog(observed, expected, min_value, name):
	obs = np.reshape(observed, -1)
	exp = np.reshape(expected, -1)
	not_zeros = exp > 0
	assert sum(not_zeros) > 1, "expected probabilities are all 0."
	subs_obs = obs[not_zeros]
	subs_exp = exp[not_zeros]		
	log_obs = np.full(subs_obs.shape, min_value)
	np.log(subs_obs, out=log_obs, where=subs_obs>0)
	log_exp=np.log(subs_exp)
	fig, ax = plt.subplots(figsize=(10,10))
	ax.scatter(x=-log_exp, y=-log_obs)
	plt.axhline(y=min_value, color='g', linestyle='dashed')
	plt.axvline(x=min_value, color='g', linestyle='dashed')
	ax.plot([0,min_value+5],[0,min_value+5])
	ax.set_xlim((0,min_value+1))
	ax.set_ylim((0,min_value+1));
	ax.set_xlabel('lnP_gimble')
	ax.set_ylabel('lnP_sims')
	ax.set_title(name)
	ax.figure.savefig(f'tests/output/scatter_{name}.png', dpi=300)
	plt.clf()