import sys
import numpy as np
import msprime
import allel
import itertools
import collections
import lib.gimble

# needs
#lib.gimble.blocks_to_arrays
#lib.gimble._return_np_type
#lib.gimble.fix_pos_array

def get_sim_args_by_replicate_idx(kwargs):
    print("[+] Preparing simulations...")
    # MODIFIERS : windows/recombination_rates, demographies (if grid)
    # CONSTANTS : blocks, samples, ploidy, block_length, comparisons, max_k, mutation_rate, mutation_model, discrete
    # VARIABLES : seeds, recombination_rates, demographies (if grid)
    # sim_key: the simulation
    global CONSTANT_PARAMS
    CONSTANT_PARAMS = {
        'blocks': kwargs['blocks'],
        'samples': {'A': kwargs['samples_A'], 'B': kwargs['samples_B']},
        'ploidy': 2,
        'block_length': kwargs['block_length'],
        'comparisons': list(itertools.product(range(kwargs['samples_A']),
                range(kwargs['samples_A'], kwargs['samples_A'] + kwargs['samples_B']))),
        'max_k': kwargs['kmax'],
        'mutation_rate': kwargs['mu'],
        'mutation_model': 'infinite_alleles',
        'discrete_genome': (not kwargs['continuous_genome'])
        }
    simulate_jobs_by_replicate_idx = collections.defaultdict(list)
    recombination_rates = kwargs['recombination_rate']
    demographies = [demography.get_demography() for demography in kwargs['demographies']]
    # print([demography.get_parameter_dict() for demography in kwargs['demographies']])
    for replicate_idx in range(kwargs['replicates']):
        ancestry_seeds = kwargs['ancestry_seeds_by_replicate'][replicate_idx]
        mutation_seeds = kwargs['mutation_seeds_by_replicate'][replicate_idx]
        for window_idx, demography, recombination_rate, ancestry_seed, mutation_seed in zip(
            range(kwargs['windows']), 
            demographies,
            recombination_rates, 
            ancestry_seeds, 
            mutation_seeds):
            simulate_jobs_by_replicate_idx[replicate_idx].append([
                window_idx, 
                replicate_idx,
                ancestry_seed,
                mutation_seed,
                demography,
                recombination_rate,
                ])
    return simulate_jobs_by_replicate_idx


def simulate_call(simulate_job):
    '''simulate call for 1 window'''
    window_idx, replicate_idx, ancestry_seed, mutation_seed, demography, recombination_rate = simulate_job
    blocks = CONSTANT_PARAMS['blocks']
    samples = CONSTANT_PARAMS['samples']
    ploidy = CONSTANT_PARAMS['ploidy']
    block_length = CONSTANT_PARAMS['block_length']
    comparisons = CONSTANT_PARAMS['comparisons']
    max_k = CONSTANT_PARAMS['max_k'] 
    mutation_rate = CONSTANT_PARAMS['mutation_rate'] 
    mutation_model = CONSTANT_PARAMS['mutation_model'] 
    discrete_genome = CONSTANT_PARAMS['discrete_genome']
    sequence_length = CONSTANT_PARAMS['block_length'] * CONSTANT_PARAMS['blocks']
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

def generate_bsfs(genotype_matrix, positions, comparisons, max_k, blocks, total_length):
    sa_genotype_array = allel.GenotypeArray(genotype_matrix)
    num_comparisons = len(comparisons)
    result = np.zeros((num_comparisons, blocks, len(max_k)), dtype=np.int64)
    # generate all comparisons
    for idx, pair in enumerate(comparisons):
        block_sites = np.arange(total_length, dtype=np.int64).reshape(blocks, -1)
        subset_genotype_array = sa_genotype_array.subset(sel0=None, sel1=pair) 
        *_, variation = lib.gimble.blocks_to_arrays(block_sites, subset_genotype_array, positions)
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

def get_genotypes(ts, ploidy, num_samples):
    if ts.num_mutations == 0:
        return np.zeros((1, num_samples, ploidy), dtype=np.uint8)
    shape = (ts.num_sites, num_samples, ploidy)
    return np.reshape(ts.genotype_matrix(), shape)