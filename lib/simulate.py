from docopt import docopt
import sys
import numpy as np
import msprime
import allel
import multiprocessing
import contextlib
from tqdm import tqdm
import itertools
import lib.gimble
from functools import partial

def dict_product(d):
    return (dict(zip(d, x)) for x in itertools.product(*d.values()))

def expand_params(d):
    for key, value in d.items():
        if len(value)>1:
            assert len(value)==3, 'MIN, MAX and STEPSIZE need to be specified'
            d[key]=np.arange(value[0], value[1], value[2])

def run_sim(parameterObj):
    threads = parameterObj.threads
    ploidy = parameterObj._config["ploidy"]
    params = parameterObj._config["parameters"]
    blocks = parameterObj._config["blocks"]
    blocklength = parameterObj._config["blocklength"]
    windows = parameterObj._config["windows"]
    replicates = parameterObj._config['replicates']
    if windows:
        print(f"[+] simulating {replicates} replicate(s) of {windows} window(s) of {blocks} block(s) each")
        total_length = blocklength*blocks*windows
    else:
        print(f"[+] simulating {replicates} replicate(s) of {blocks} independent block(s)")
        total_length=blocklength

    expand_params(params)
    sim_configs = dict_product(params)
    msprime_configs = (make_sim_configs(config, ploidy) for config in sim_configs)
    print("[+] running simulations")
    with multiprocessing.Pool(processes=threads) as pool:
        #param_configs = partial()
        pool.map(run_ind_sim, msprime_configs)

def make_sim_configs(params, ploidy):
    sample_size_A = params['sample_size_A']
    sample_size_B = params['sample_size_B']   
    C_A = params['C_A']
    C_B = params['C_B']
    total_samples = sample_size_A + sample_size_B
    demographic_events = None
    population_configurations = [
        msprime.PopulationConfiguration(
            sample_size=sample_size_A * ploidy, initial_size=C_A
        ),
        msprime.PopulationConfiguration(
            sample_size=sample_size_B * ploidy, initial_size=C_B
        ),
    ]
    migration_matrix = np.zeros((3, 3))
    if 'M_A_B' in params:
        # migration A to B backwards
        migration_matrix[1, 0] = params['M_A_B']
    if 'M_B_A' in params:
        # migration B to A
        migration_matrix[0, 1] = params['M_B_A']
    if 'C_A_B' in params:
        C_AB = params['C_A_B'] if params['C_A_B'] else C_A+C_B
        population_configurations += [
            msprime.PopulationConfiguration(sample_size=0, initial_size=C_AB),
        ]
        # demographic events: specify in the order they occur backwards in time
        demographic_events = [
            msprime.MassMigration(
                time=params['T'], source=0, destination=2, proportion=1.0
            ),
            msprime.MassMigration(
                time=params['T'], source=1, destination=2, proportion=1.0
            ),
            msprime.MigrationRateChange(params["T"], 0),
        ]

    return (population_configurations, demographic_events, migration_matrix) 

def run_ind_sim(msprime_config):
    population_configurations, demographic_events, migration_matrix = msprime_config
    mutation_rate = 0.0001
    recombination_rate = 0.0
    length = 1000
    ploidy=2
    total_samples = 2 

    ts = msprime.simulate(
        length=length,
        recombination_rate=recombination_rate,
        population_configurations=population_configurations,
        demographic_events=demographic_events,
        migration_matrix=migration_matrix,
        mutation_rate=mutation_rate,
    )

    """
    #with msprime 1.0
    ts = run_ind_sim(
        population_configurations=population_configurations,
        demographic_events=demographic_events,
        migration_matrix=migration_matrix,
        length=blocklength,
        mutation_rate=params["theta"],
        recombination_rate=0.0,
    )
    tsm = msprime.mutate(ts, rate=mutation_rate, discrete=True)
    positions = np.array([site.position for site in tsm.sites()])
    """
    #with infinite sites = pre-msprime 1.0
    positions = np.array([int(site.position) for site in ts.sites()])
    new_positions = fix_pos_array(positions)
    print("[+] finished simulations")
    genotype_matrix = get_genotypes(ts, ploidy, total_samples)
    #sa_genotype_array = allel.GenotypeArray(genotype_matrix)
    print("[+] generated genotype matrix")
    #block_sites = np.arange(blocks*blocklength).reshape(blocks,blocklength)
    #block_sites_variant_bool = np.isin(block_sites, new_positions, assume_unique=True)
    #block_sites = lib.gimble.genotype_to_mutype_array(sa_genotype_array, block_sites_variant_bool, block_sites, debug=True)
    #multiallelic, missing, monomorphic, variation = lib.gimble.block_sites_to_variation_arrays(block_sites)

def get_genotypes(ts, ploidy, num_samples):
    shape = (ts.num_mutations, num_samples, ploidy)

    return np.reshape(ts.genotype_matrix(), shape)

def fix_pos_array(pos_array):
    #print('pos_array_0', pos_array)
    # get boolean array for first and subsequent duplicates (True) (required sorted) 
    idxs = np.insert((np.diff(pos_array)==0).astype(np.bool), 0, False)
    if np.any(idxs): 
        # if there are duplicates, get new values by incrementing by one
        new_values = pos_array[idxs]+1
        # get non-duplicate values
        uniq_values = pos_array[~idxs]
        # insert new_values in non-duplicated values (required sorted)
        new_idxs = np.searchsorted(uniq_values, new_values)
        # recursive call
        return fix_pos_array(np.sort(np.insert(uniq_values, new_idxs, new_values)))
    else:
        # if there are no duplicated values
        #print('pos_array_1', pos_array)
        return pos_array
        