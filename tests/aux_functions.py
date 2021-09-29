import numpy as np
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use("Agg")
import scipy.stats
import os
import collections
import pathlib
import lib.simulate

def get_test_config_file(task=None, model=None, label=None):
    if all([task, model, label]):
        return str(pathlib.Path.cwd() / pathlib.Path('tests/input/gimble.%s.%s.%s.config.ini' % (task, model, label)))
    return None

def get_numpy_array(model, label):
    dir_str = str(pathlib.Path.cwd() / pathlib.Path('tests/input/%s.%s.npy' % (model, label)))
    return np.squeeze(np.load(dir_str))

def sim_combi_gt_matrix(p=2, s=2, alleles=[0, 1]):
    '''
    generates numpy-genotype-array all possible combinations of elements in `alleles`
    WARNING: only works for (p*s)<30 
    ''' 
    arrays = p*s*[np.array(alleles)]
    la = len(arrays)
    dtype = np.result_type(*arrays)
    arr = np.empty([len(a) for a in arrays] + [la], dtype=dtype)
    for i, a in enumerate(np.ix_(*arrays)):
        arr[...,i] = a
    return arr.reshape(-1, s, p)

def sim_random_gt_matrix(p=2, s=2, n=10, alleles=[0, 1]):
    '''generates random numpy-genotype-array with `n` genotypes''' 
    rng = np.random.default_rng(12345)
    return rng.choice(alleles, size=(n,s,p))

def sim_sites(n):
    '''draws n integers from (0, n)'''
    rng = np.random.default_rng()
    return np.unique(np.sort(rng.integers(low=0, high=n, size=n)))

def sim_block_variation(b=10, kmax=3):
    '''simulates blocks variation array'''
    rng = np.random.default_rng(12345)
    out =  rng.choice(range(kmax+1), size=(b, 4))
    # FGVs need to be fixed
    rows_to_fix = (out[:,2] > 0) & (out[:,3] > 0)
    cols_to_fix = rng.choice([2,3], size=rows_to_fix[rows_to_fix==True].shape[0])
    out[rows_to_fix, cols_to_fix] = 0
    return out 

def sim_window_variation(b=10, w=10, kmax=3):
    '''simulates windows variation array'''
    rng = np.random.default_rng(12345)
    out = rng.choice(range(kmax+1), size=(w, b, 4))
    # FGVs need to be fixed
    rows_to_fix = (out[:,:,2] > 0) & (out[:,:,3] > 0)
    cols_to_fix = rng.choice([2,3], size=rows_to_fix[rows_to_fix==True].shape[0])
    out[rows_to_fix, cols_to_fix] = 0
    return out 

def bsfs_to_2d(bsfs):
    '''needed for testing that tally-arrays and bsfs are identical'''
    if not np.any(bsfs):
        return None
    non_zero_idxs = np.nonzero(bsfs)
    if bsfs.ndim == 4: # blocks
        return np.concatenate([bsfs[non_zero_idxs].reshape(non_zero_idxs[0].shape[0], 1), np.array(non_zero_idxs, dtype=np.uint64).T], axis=1)
    elif bsfs.ndim == 5: # windows
        non_zero_idxs_array = np.array(non_zero_idxs, dtype=np.uint64).T
        first = non_zero_idxs_array[:,0].reshape(non_zero_idxs[0].shape[0], 1)
        second = bsfs[non_zero_idxs].reshape(non_zero_idxs[0].shape[0], 1)
        third = non_zero_idxs_array[:,1:]
        return np.concatenate([first, second, third], axis=1)
    else:
        raise ValueError('bsfs_to_2d: bsfs.ndim must be 4 (blocks) or 5 (windows)')

def chisquare(observed, expected, p_th=0.05, recombination=False, all_sims=False):
    #expected counts need to be larger than 5 for chisquare
    #combining bins with less than 5 counts into one bin.
    obs = np.reshape(observed, -1)
    exp = np.reshape(expected, -1)
    if not (recombination or all_sims):
        assert np.all(obs[exp == 0] == 0), "chisquare error: at least one cell with expected frequency 0 has observed values"
    #bin all values with counts smaller than 5
    binning_idxs = np.digitize(exp, np.array([5]), right=True)
    exp_binned = exp[binning_idxs==1]
    obs_binned = obs[binning_idxs==1]
    not_zeros = exp_binned > 0
    if sum(not_zeros)<1:
        assert False #expected probabilities are all 0
    else:
        chisq = scipy.stats.chisquare(obs_binned[not_zeros], exp_binned[not_zeros])
        print("chisquare value:", chisq)
        print("exp", exp_binned)
        print("obs", obs_binned)
        assert chisq.pvalue > p_th

def chisquare_contingency(observed1, observed2, p_th=0.05):
    obs1 = np.reshape(observed1, -1)
    obs2 = np.reshape(observed2, -1)
    all_obs = np.vstack((obs1, obs2))
    assert not np.all(all_obs==0), 'all counts are 0'
    binning_idxs = np.digitize(all_obs,np.array([5]), right=True)
    binning_across_cols = (binning_idxs==0).any(axis=0)
    obs_all_smaller = np.sum(all_obs[:,binning_across_cols==True],axis=1)
    if np.all(obs_all_smaller==0):
        obs_all_binned = all_obs[:,binning_across_cols==False]
    else:
        obs_all_binned = np.hstack((all_obs[:,binning_across_cols==False], obs_all_smaller.reshape((2,1))))
    chi2, p, dof, ex = scipy.stats.chi2_contingency(obs_all_binned)
    print("chisquare value:", p)
    print("obs_all_binned")
    assert p>p_th

def bonferroni(p_threshold, pvals):
    assert min(pvals) > p_threshold / len(pvals)

def scatter_loglog(observed, expected, min_value, name, num_replicates, xlabel='gimble', ylabel='sims'):
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
    ax.set_ylim((0,min_value+1))
    ax.set_xlabel(f'lnP_{xlabel}')
    ax.set_ylabel(f'lnP_{ylabel}')
    #add confidence interval
    #pick single point
    p = 0.01
    ref_point = -np.log(p)
    ci = np.array(scipy.stats.binom.interval(1-0.05/subs_exp.size, num_replicates, p))
    ci = -np.log(ci/num_replicates)
    ax.plot([0, ref_point],[0,ci[0]], 'k--')
    ax.plot([0, ref_point],[0,ci[1]], 'k--')
    ax.set_title(name)
    if not os.path.isdir('tests/output'):
        os.mkdir('tests/output')
    ax.figure.savefig(f'tests/output/scatter_{name}.png', dpi=300)
    plt.clf()

def sim_ETPs(global_info, sim_configs, n, blocks, chunks, replicates):
    threads = 1
    add_global_info = {
                        'blocks':blocks, 
                        'chunks':chunks, 
                        'sample_pop_sizes':[n for _ in global_info['sample_pop_ids']],
                        'replicates':replicates
                        }
    all_interpop_comparisons = lib.simulate.all_interpopulation_comparisons(*add_global_info['sample_pop_sizes'])
    sim_global_info = {**global_info, **add_global_info}
    simmed_ETPs = lib.simulate.run_sims(sim_configs, sim_global_info, all_interpop_comparisons, chunks, threads, disable_tqdm=True)
    return simmed_ETPs

def make_single_demography(config):
    print(isinstance(config, dict))
    print(config['parameters_grid_points'])
    graph_list = lib.gimble.config_to_demes_graph(config, None)
    #return graph_list[0]
    return next(lib.simulate.make_demographies(graph_list))

def downsample(A, N):
    distribution = [i for i, j in enumerate(A) for _ in range(j)]
    sample = collections.Counter(random.sample(distribution, N))
    return [sample[i] for i in range(len(A))]