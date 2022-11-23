from docopt import docopt
import sys, os
import numpy as np
import msprime
import allel
import zarr
import concurrent.futures
import contextlib
from tqdm import tqdm
import itertools
import lib.gimble
import pandas as pd
from functools import partial
import collections
import multiprocessing

# def run_sims(config, threads=1, discrete=True, disable_tqdm=False): 
#     samples = {
#                'A': config['simulate']['sample_size_A'], 
#                'B': config['simulate']['sample_size_B']
#                }
#     mutation_model = 'infinite_alleles'
#     recombination_rates = config['simulate']['recombination_rate']
#     demographies = config['demographies']
#     seeds = config['seeds']
#     if isinstance(recombination_rates, float):
#         recombination_rates = itertools.repeat(recombination_rates)
#     for idx, (demography, rec_rate, seed) in enumerate(tqdm(
#         zip(demographies, recombination_rates, seeds),
#         desc='Overall simulation progress',
#         ncols=100, 
#         unit_scale=True, 
#         total=config['parameters_grid_points'], 
#         disable=disable_tqdm
#         )
#     ):
#         if threads > 1:
#             result_list = run_sim_parallel(
#                 seed, 
#                 rec_rate, 
#                 samples, 
#                 demography, 
#                 config['simulate']['ploidy'], 
#                 config['simulate']['block_length'], 
#                 config['simulate']['comparisons'], 
#                 config['max_k'],
#                 config['mu']['mu'],
#                 mutation_model,
#                 discrete,
#                 config['simulate']['blocks_per_replicate'], 
#                 disable_tqdm,
#                 idx,
#                 threads
#                 )
#         else:
#             result_list = run_sim_serial(
#                 seed, 
#                 rec_rate, 
#                 samples, 
#                 demography, 
#                 config['simulate']['ploidy'], 
#                 config['simulate']['block_length'], 
#                 config['simulate']['comparisons'], 
#                 config['max_k'], 
#                 config['mu']['mu'],
#                 mutation_model,
#                 discrete,
#                 config['simulate']['blocks_per_replicate'], 
#                 disable_tqdm,
#                 idx
#                 )
#         chunks = config['simulate']['blocks_per_replicate'].size
#         result_list = _combine_chunks(result_list, chunks)
#         yield np.array(result_list)

# def sim_worker(seed, blocks, recombination_rate, samples, demography, ploidy, block_length, comparisons, max_k, mutation_rate, mutation_model, discrete):
#     ancestry_seed, mutation_seed = seed
#     sequence_length = block_length * blocks
#     #run simulation:
#     ts = msprime.sim_ancestry(
#         samples=samples, 
#         demography = demography, 
#         ploidy=ploidy,
#         sequence_length=sequence_length,
#         discrete_genome=discrete,
#         recombination_rate=recombination_rate,
#         random_seed=ancestry_seed
#         )
#     ts = msprime.sim_mutations(
#         ts, 
#         rate=mutation_rate, 
#         random_seed=mutation_seed, 
#         discrete_genome=discrete, 
#         model=mutation_model
#         )
#     num_samples = sum(samples.values())
#     #make bsfs
#     if discrete:
#         positions = np.array([site.position for site in ts.sites()], dtype=np.int64)
#     else:
#         positions, sequence_length = infinite_sites(ts, blocks, sequence_length)
#     genotype_matrix = get_genotypes(ts, ploidy, num_samples)
#     bsfs = generate_bsfs(genotype_matrix, positions, comparisons, max_k, blocks, sequence_length)
#     return bsfs

def get_sim_args(config):
    print("[+] Preparing simulations...")
    # MODIFIERS : windows/recombination_rates, demographies (if grid)
    # CONSTANTS : blocks, samples, ploidy, block_length, comparisons, max_k, mutation_rate, mutation_model, discrete
    # VARIABLES : seeds, recombination_rates, demographies (if grid)
    # sim_key: the simulation

    # Query tally: 
    # - window_idx, replicate_idx, parameters of demography, seeds, ...
    # Query gridsearch (best model, slicing):
    # - window_idx, replicate_idx, recombination, as in data
    constant_params = {
        'blocks': config['simulate']['blocks'],
        'samples': {'A': config['simulate']['sample_size_A'], 'B': config['simulate']['sample_size_B']},
        'ploidy': config['simulate']['ploidy'],
        'block_length': config['simulate']['block_length'],
        'comparisons': config['simulate']['comparisons'],
        'max_k': config['max_k'],
        'mutation_rate': config['mu']['mu'],
        'mutation_model': 'infinite_alleles',
        'discrete_genome': config['simulate']['discrete_genome']
        }
    simulate_jobs = []
    for replicate_idx in range(config['simulate']['replicates']):
        ancestry_seeds = config['simulate']['ancestry_seeds_by_replicate'][replicate_idx]
        mutation_seeds = config['simulate']['mutation_seeds_by_replicate'][replicate_idx]
        demographies = config['simulate']['demographies']
        recombination_rates = config['simulate']['recombination_rate']
        for window_idx, demography, recombination_rate, ancestry_seed, mutation_seed in zip(
            range(config['simulate']['windows']), 
            demographies,
            recombination_rates, 
            ancestry_seeds, 
            mutation_seeds):
            #idx += 1
            simulate_jobs.append([
                window_idx, 
                replicate_idx,
                ancestry_seed,
                mutation_seed,
                demography.demography,
                recombination_rate,
                constant_params['blocks'],
                constant_params['samples'],
                constant_params['ploidy'],
                constant_params['block_length'],
                constant_params['comparisons'],
                constant_params['max_k'], 
                constant_params['mutation_rate'], 
                constant_params['mutation_model'], 
                constant_params['discrete_genome']])
    return simulate_jobs
    
def simulate_call(simulate_job):
    window_idx, replicate_idx, ancestry_seed, mutation_seed, demography, recombination_rate, blocks, samples, ploidy, block_length, comparisons, max_k, mutation_rate, mutation_model, discrete_genome = simulate_job
    sequence_length = block_length * blocks
    #run simulation:
    ts = msprime.sim_ancestry(
        samples=samples, 
        demography = demography, 
        ploidy=ploidy,
        sequence_length=sequence_length,
        discrete_genome=discrete_genome,
        recombination_rate=recombination_rate,
        random_seed=ancestry_seed
        )
    ts = msprime.sim_mutations(
        ts, 
        rate=mutation_rate, 
        random_seed=mutation_seed, 
        discrete_genome=discrete_genome, 
        model=mutation_model
        )
    num_samples = sum(samples.values())
    #make bsfs
    if discrete_genome:
        positions = np.array([site.position for site in ts.sites()], dtype=np.int64)
    else:
        positions, sequence_length = infinite_sites(ts, blocks, sequence_length)
    genotype_matrix = get_genotypes(ts, ploidy, num_samples)
    bsfs = generate_bsfs(genotype_matrix, positions, comparisons, max_k, blocks, sequence_length)
    return {'window_idx': window_idx, 'replicate_idx': replicate_idx, 'bsfs': bsfs}

# def simulate_function(nlopt_values, grad, gfEvaluatorObj, dataset, dataset_idx, config, nlopt_iterations, verbose):
    # global NLOPT_LOG_QUEUE
    # nlopt_traceback = None
    # nlopt_iterations[dataset_idx] += 1 
    # scaled_values_by_symbol, scaled_values_by_parameter, unscaled_values_by_parameter, block_length = scale_nlopt_values(nlopt_values, config)
    # try:
    #     ETPs = gfEvaluatorObj.evaluate_gf(scaled_values_by_symbol, scaled_values_by_symbol[sage.all.SR.var('theta')])
    #     likelihood = calculate_composite_likelihood(ETPs, dataset)
    # except Exception as exception:
    #     nlopt_log_iteration_tuple = get_nlopt_log_iteration_tuple(dataset_idx, nlopt_iterations[dataset_idx], block_length, 'N/A', scaled_values_by_parameter, unscaled_values_by_parameter)
    #     nlopt_traceback = traceback.format_exc(exception)
    #     NLOPT_LOG_QUEUE.put((nlopt_log_iteration_tuple, nlopt_traceback))
    # nlopt_log_iteration_tuple = get_nlopt_log_iteration_tuple(dataset_idx, nlopt_iterations[dataset_idx], block_length, likelihood, scaled_values_by_parameter, unscaled_values_by_parameter)
    # NLOPT_LOG_QUEUE.put((nlopt_log_iteration_tuple, nlopt_traceback))
    # #if verbose:
    # #    elapsed = lib.gimble.format_time(timer() - start_time)
    # #    #process_idx = multiprocessing.current_process()._identity
    # #    #print_nlopt_line(dataset_idx, nlopt_iterations[dataset_idx], likelihood, unscaled_values_by_parameter, elapsed, process_idx)
    # #    print_nlopt_line(dataset_idx, nlopt_iterations[dataset_idx], likelihood, unscaled_values_by_parameter, elapsed)
    # return likelihood

SIMULATE_QUEUE = multiprocessing.Queue()
SIMULATE_ITERATIONS_MANAGER = multiprocessing.Manager()

@contextlib.contextmanager
def poolcontext(*args, **kwargs):
    pool = multiprocessing.Pool(*args, **kwargs)
    yield pool
    pool.terminate()

# def new_simulate_collections(config):
#     simulate_jobs = get_sim_args(config)
#     simulate_windows_by_replicate = collections.defaultdict(list) # unsorted list
#     if config['num_cores'] <= 1:
#         for simulate_job in tqdm(simulate_jobs):
#             simulate_window = simulate_call(simulate_job)
#             #print("simulate_window", simulate_window)
#             simulate_windows_by_replicate[simulate_window['replicate_idx']].append(simulate_window)
#     else:
#         with poolcontext(processes=config['num_cores']) as pool:
#             for simulate_window in tqdm(pool.imap_unordered(simulate_call, simulate_jobs), total=len(simulate_jobs)):
#                 simulate_windows_by_replicate[simulate_window['replicate_idx']].append(simulate_window)
    
#     # sanitize results
#     tally_by_replicate_idx = {}
#     for replicate_idx, windows in simulate_windows_by_replicate.items():
#         tallies = [tally['bsfs'] for tally in sorted(windows, key=lambda k: (k['window_idx']))]

#         tally_by_replicate_idx[replicate_idx] = np.array(tallies)
#     return tally_by_replicate_idx

def new_simulate(config):
    simulate_jobs = get_sim_args(config)
    simulate_windows_by_replicate = collections.defaultdict(list) # unsorted list
    tallies = np.zeros(tuple([config['simulate']['replicates'], config['simulate']['windows']] + list(config['max_k'] + 2)))
    # save tallies as numpy arrays directly!!!
    if config['num_cores'] <= 1:
        for simulate_job in tqdm(simulate_jobs):
            simulate_window = simulate_call(simulate_job)
            #print("simulate_window", simulate_window)
            tallies[simulate_window['replicate_idx']][simulate_window['window_idx']] = simulate_window['bsfs']
    else:
        with poolcontext(processes=config['num_cores']) as pool:
            for simulate_window in tqdm(pool.imap_unordered(simulate_call, simulate_jobs), total=len(simulate_jobs)):
                tallies[simulate_window['replicate_idx']][simulate_window['window_idx']] = simulate_window['bsfs']
    return tallies

# def run_sim_parallel(seeds, recombination_rate, samples, demography, ploidy, block_length, comparisons, max_k, mutation_rate, mutation_model, discrete, blocks_per_replicate, disable_tqdm, idx, threads):
#     with concurrent.futures.ProcessPoolExecutor(max_workers=threads) as pool:
#         run_sims_specified = partial(
#         sim_worker,
#         recombination_rate=recombination_rate,
#         samples=samples,
#         demography=demography,
#         ploidy=ploidy,
#         block_length=block_length,
#         comparisons=comparisons,
#         max_k=max_k,
#         mutation_rate=mutation_rate,
#         mutation_model=mutation_model,
#         discrete=discrete
#     )
    
#         result_list = list(
#                         tqdm(
#                             pool.map(
#                                 run_sims_specified, 
#                                 seeds, 
#                                 itertools.cycle(blocks_per_replicate)
#                                 ),
#                             desc=f'running parameter combination {idx}',
#                             ncols=100, 
#                             unit_scale=True, 
#                             total=len(seeds),
#                             disable=disable_tqdm)
#                         )
#     return result_list
       
# def run_sim_serial(seeds, recombination_rate, samples, demography, ploidy, block_length, comparisons, max_k, mutation_rate, mutation_model, discrete, blocks_per_replicate, disable_tqdm, idx):
#     # [SIMULATION]
#     shape = np.insert(max_k+2, 0, seeds.shape[0]) 
#     result_list = np.zeros(shape, dtype=lib.gimble._return_np_type(np.sum(blocks_per_replicate)))
#     for sub_idx, (seed, block_per_replicate) in enumerate(
#                         tqdm(
#                             zip(
#                                 seeds, 
#                                 itertools.cycle(blocks_per_replicate)
#                                 ), 
#                             desc=f'running parameter combination {idx}',
#                             ncols=100, 
#                             unit_scale=True, 
#                             disable=disable_tqdm,
#                             total=len(seeds)
#                             )
#                         ):
#         result_list[sub_idx] = sim_worker(
#                 seed=seed,
#                 blocks=block_per_replicate,
#                 recombination_rate=recombination_rate,
#                 samples=samples,
#                 demography=demography,
#                 ploidy=ploidy,
#                 block_length=block_length,
#                 comparisons=comparisons,
#                 max_k=max_k,
#                 mutation_rate=mutation_rate,
#                 mutation_model=mutation_model,
#                 discrete=discrete
#             )
#     return result_list

def generate_bsfs(genotype_matrix, positions, comparisons, max_k, blocks, total_length):
    sa_genotype_array = allel.GenotypeArray(genotype_matrix)
    num_comparisons = len(comparisons)
    result = np.zeros((num_comparisons, blocks, len(max_k)), dtype=np.int64)
    # generate all comparisons

    for idx, pair in enumerate(comparisons):
        block_sites = np.arange(total_length, dtype=np.int64).reshape(blocks, -1)
        new_positions_variant_bool = np.isin(
            positions, block_sites, assume_unique=True
            )
        subset_genotype_array = sa_genotype_array.subset(new_positions_variant_bool, pair) #all variants are included
        *redundant, variation = lib.gimble.blocks_to_arrays(block_sites, subset_genotype_array, positions)
        result[idx] = variation
    
    result = result.reshape(-1, result.shape[-1])
    # count mutuples (clipping at k_max, if supplied)
    mutuples, counts = np.unique(np.clip(result, 0, max_k+1), return_counts=True, axis=0)
    # define out based on max values for each column
    dtype = lib.gimble._return_np_type(counts)
    out = np.zeros(tuple(max_k + 2), dtype)
    # assign values
    out[tuple(mutuples.T)] = counts
    return out

def infinite_sites(ts, blocks, total_length): 
    positions = np.array([int(site.position) for site in ts.sites()])
    new_positions = lib.gimble.fix_pos_array(positions)
    if ts.num_sites>0 and new_positions[-1]>=total_length:
        blocklength = int(np.ceil(new_positions[-1]/blocks))
        total_length = blocks*blocklength
    if ts.num_sites==0:
        new_positions = [0,]
    return (new_positions.astype(np.int64), total_length)

def _combine_chunks(result_list, chunks):
    if chunks>1:
        return np.array([np.add.reduce(m) for m in np.split(result_list,list(range(0,len(result_list),chunks))[1:])])
    else:
        return np.array(result_list)

def get_genotypes(ts, ploidy, num_samples):
    if ts.num_mutations == 0:
        return np.zeros((1,num_samples, ploidy), dtype=np.uint8)
    shape = (ts.num_sites, num_samples, ploidy)
    return np.reshape(ts.genotype_matrix(), shape)

def all_interpopulation_comparisons(*popsizes):
    popA, popB, *rest = popsizes
    if len(rest)>0:
        raise ValueError("More than 2 population sizes were provided to simulate. We cannot cope with that just yet.")
    return list(itertools.product(range(popA), range(popA, popA + popB)))

def make_demographies(config):
    #return (msprime.Demography.from_demes(graph) for graph in lib.gimble.config_to_demes_graph(config)) # generator
    return [msprime.Demography.from_demes(graph) for graph in lib.gimble.config_to_demes_graph(config)] # list !!!