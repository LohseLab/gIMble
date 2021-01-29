import numpy as np
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use("Agg")
import scipy.stats as stats

def chisquare(observed, expected, p_th=0.05):
	obs = np.reshape(observed, -1)
	exp = np.reshape(expected, -1)
	assert np.all(obs[exp == 0] == 0)
	not_zeros = exp > 0
	if sum(not_zeros) > 1:
		chisq = stats.chisquare(obs[not_zeros], exp[not_zeros])
		assert chisq.pvalue > p_th

def bonferroni(p_threshold, pvals):
	assert min(pvals) > p_threshold / len(pvals)

def scatter_xy(observed, expected, min_value, name):
	obs = np.reshape(observed, -1)
	exp = np.reshape(expected, -1)
	not_zeros = exp > 0
	assert sum(not_zeros) > 1, "expected probabilities are all 0."
	subs_obs = obs[not_zeros]
	subs_exp = exp[not_zeros]		
	log_obs = np.full(subs_obs.shape, min_value)
	np.log(subs_obs, out=log_obs)
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