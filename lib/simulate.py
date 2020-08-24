from docopt import docopt
import sys, os
import numpy as np
import msprime
import allel
import zarr
import multiprocessing
import contextlib
from tqdm import tqdm
from tqdm.auto import trange
import itertools
import lib.gimble
import pandas as pd
from functools import partial
import collections

def run_sim(parameterObj, gimbleStore):
    threads = parameterObj.threads
    global_info = parameterObj.config["simulations"]
    ploidy = global_info["ploidy"]
    blocks = global_info["blocks"]
    blocklength = global_info["blocklength"]
    replicates = global_info["replicates"]
    sim_configs = parameterObj.grid
    global_info['sample_pop_ids'] = parameterObj.config['populations']['sample_pop_ids']
    A,B = global_info['sample_pop_ids']
    global_info['reference_pop'] = parameterObj.config['populations']['reference_pop']
    global_info['mu'] = parameterObj.config['mu']['mu']
    
    msprime_configs = (make_sim_configs(config, global_info) for config in sim_configs)
    all_interpop_comparisons = all_interpopulation_comparisons(
        global_info[f"sample_size_{A}"], global_info[f"sample_size_{B}"]
    )
    
    print(f"[+] simulating {replicates} replicate(s) of {blocks} block(s) for {len(sim_configs)} parameter combinations")
    #with tqdm(total=replicates*len(sim_configs), desc="[%] running sims ", ncols=100, unit_scale=True) as pbar:
    for idx, (config, zarr_attrs) in enumerate(tqdm(zip(msprime_configs, sim_configs),desc='Overall simulation progress',ncols=100, unit_scale=True, total=len(sim_configs))):
        seeds = np.random.randint(1, 2 ** 32, replicates)
        result_list = []
        if threads > 1:
            with multiprocessing.Pool(processes=threads) as pool:

                run_sims_specified = partial(
                run_ind_sim,
                msprime_config=config,
                ploidy=ploidy,
                blocks=blocks,
                blocklength=blocklength,
                comparisons=all_interpop_comparisons,
            )
    
                result_list = list(tqdm(pool.imap(run_sims_specified, seeds),desc=f'running parameter combination {idx}',ncols=100, unit_scale=True, total=replicates))
        else:
            for seed in tqdm(seeds,desc=f'running parameter combination {idx}',ncols=100, unit_scale=True):
                result_list.append(
                    run_ind_sim(
                        seed=seed,
                        msprime_config=config,
                        ploidy=ploidy,
                        blocks=blocks,
                        blocklength=blocklength,
                        comparisons=all_interpop_comparisons
                    )
                )
            
        name = f"parameter_combination_{idx}"
        run_count = gimbleStore.data['sims'].attrs['run_count']
        g = gimbleStore.data[f'sims/run_{run_count}'].create_dataset(name, data=np.array(result_list), overwrite=True)
        g.attrs.put(zarr_attrs)
        g.attrs['seeds']=tuple([int(s) for s in seeds])
            
def make_sim_configs(params, global_info):
    A, B = global_info["sample_pop_ids"]
    sample_size_A = global_info[f"sample_size_{A}"]
    sample_size_B = global_info[f"sample_size_{B}"]
    num_samples = sample_size_A + sample_size_B
    C_A = params[f"Ne_{A}"]
    C_B = params[f"Ne_{B}"]
    if f"Ne_{A}_{B}" in params:
        C_AB = params[f"Ne_{A}_{B}"]
    elif f"Ne_{B}_{A}" in params:
        C_AB = params[f"Ne_{B}_{A}"]
    else: 
        C_AB = params[f"Ne_{global_info['reference_pop']}"]
    mu = global_info["mu"]
    rec_rate = params["recombination"]

    population_configurations = [
        msprime.PopulationConfiguration(
            sample_size=sample_size_A * global_info["ploidy"], initial_size=C_A
        ),
        msprime.PopulationConfiguration(
            sample_size=sample_size_B * global_info["ploidy"], initial_size=C_B
        ),
        msprime.PopulationConfiguration(
            sample_size=0, initial_size=C_AB
        )
    ]

    migration_matrix = np.zeros((3, 3))  # migration rate needs to be divided by 4Ne
    #migration matirx: M[i,j]=k k is the fraction of population i consisting of migrants
    # from population j, FORWARDS in time.
    #here migration is defined backwards in time
    if f"me_{A}_{B}" in params:
        # migration A to B backwards, forwards in time, migration from B to A
        migration_matrix[0, 1] = params[f"me_{A}_{B}"] #this needs to be verified
    if f"me_{B}_{A}" in params:
        # migration B to A, forwards in time, migration from A to B
        migration_matrix[1, 0] = params[f"me_{B}_{A}"]
    
    # demographic events: specify in the order they occur backwards in time
    demographic_events = []
    if params["T"]:
        demographic_events = [
            msprime.MassMigration(
                time=params["T"], source=0, destination=2, proportion=1.0
            ),
            msprime.MassMigration(
                time=params["T"], source=1, destination=2, proportion=1.0
            ),
            msprime.MigrationRateChange(params["T"], 0),
        ]

    return (
        population_configurations,
        demographic_events,
        migration_matrix,
        mu,
        num_samples,
        rec_rate,
    )


def run_ind_sim(
    seed,
    msprime_config,
    ploidy,
    blocks,
    blocklength,
    comparisons
):
    (
        population_configurations,
        demographic_events,
        migration_matrix,
        mu,
        num_samples,
        rec_rate
    ) = msprime_config
    total_length = blocks * blocklength
    ts = msprime.simulate(
        length=total_length,
        recombination_rate=rec_rate,
        population_configurations=population_configurations,
        demographic_events=demographic_events,
        migration_matrix=migration_matrix,
        mutation_rate=mu,
        random_seed=seed, #error was when 3582573439
    )
    
    """
    #with msprime 1.0 -> finite sites mutations
    ts = run_ind_sim(
        population_configurations=population_configurations,
        demographic_events=demographic_events,
        migration_matrix=migration_matrix,
        length=blocklength,
        mutation_rate=params["mu"],
        recombination_rate=0.0,
    )
    tsm = msprime.mutate(ts, rate=mutation_rate, discrete=True)
    positions = np.array([site.position for site in tsm.sites()])
    """
    # with infinite sites = pre-msprime 1.0
    positions = np.array([int(site.position) for site in ts.sites()])
    #print(f"[+] {ts.num_sites} mutation(s) along the simulated sequence")
    new_positions = lib.gimble.fix_pos_array(positions)
    if ts.num_sites>0 and any(p>=total_length for p in new_positions):
        blocklength = new_positions[-1]
        total_length = blocks*blocklength
    if ts.num_sites==0:
        new_positions = [0,]
    genotype_matrix = get_genotypes(ts, ploidy, num_samples)
    sa_genotype_array = allel.GenotypeArray(genotype_matrix)
    # always the same for all pairwise comparisons
    #print("[+] generated genotype matrix")
    # generate all comparisons
    num_comparisons = len(comparisons)
    #result = np.zeros((num_comparisons, blocks, blocklength), dtype="int8")
    result = np.zeros((num_comparisons, blocks, 4), dtype="int64") #get number of mutypes
    for idx, pair in enumerate(comparisons):
        block_sites = np.arange(total_length).reshape(blocks, blocklength)
        # slice genotype array
        #subset_genotype_array = sa_genotype_array.subset(sel1=pair)
        block_sites_variant_bool = np.isin(
            block_sites, new_positions, assume_unique=True
        )
        new_positions_variant_bool = np.isin(
            new_positions, block_sites, assume_unique=True
        )
        subset_genotype_array = sa_genotype_array.subset(new_positions_variant_bool, pair)
        #result[idx] = lib.gimble.genotype_to_mutype_array(
        #    subset_genotype_array, block_sites_variant_bool, block_sites, debug=False
        #)
        block_sites = lib.gimble.genotype_to_mutype_array(
            subset_genotype_array, block_sites_variant_bool, block_sites, debug=False
        )
        multiallelic, missing, monomorphic, variation = lib.gimble.block_sites_to_variation_arrays(block_sites)
        result[idx] = variation
    return result


def get_genotypes(ts, ploidy, num_samples):
    if ts.num_mutations == 0:
        return np.zeros((1,num_samples, ploidy), dtype='int8')
    shape = (ts.num_mutations, num_samples, ploidy)
    return np.reshape(ts.genotype_matrix(), shape)

def dict_product(d):
    if len(d)>0:
        return [dict(zip(d, x)) for x in itertools.product(*d.values())]

def all_interpopulation_comparisons(popA, popB):
    return list(itertools.product(range(popA), range(popA, popA + popB)))
