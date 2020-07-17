from docopt import docopt
import sys, os
import numpy as np
import msprime
import allel
import zarr
import multiprocessing
import contextlib
from tqdm import tqdm
import itertools
import lib.gimble
from functools import partial
import pandas as pd
import collections


def run_sim(parameterObj):
    threads = parameterObj.threads
    ploidy = parameterObj._config["ploidy"]
    params = parameterObj._config["parameters"]
    blocks = parameterObj._config["blocks"]
    blocklength = parameterObj._config["blocklength"]
    replicates = parameterObj._config["replicates"]
    zarr_store = parameterObj.zstore
    store = lib.gimble.load_store(parameterObj)
    
    expand_params(params)
    sim_configs = dict_product(params)
    msprime_configs = (make_sim_configs(config, ploidy) for config in sim_configs)
    all_interpop_comparisons = all_interpopulation_comparisons(
        params["sample_size_A"][0], params["sample_size_B"][0]
    )
    print(f"[+] simulating {replicates} replicate(s) of {blocks} block(s) for {len(sim_configs)} parameter combinations")
    store.data.require_group('sims')
    #check how many simulation runs we have already done
    #count number of groups in this level and add one
    if len(store.data['sims'])==0:
        runcount = 0
    else:
        runcount = max([int(group[0].split('_')[-1]) for group in list(store.data['sims'].groups())])
    new_run_number = runcount + 1
    root = store.data.create_group(f'sims/run_{new_run_number}')
    with tqdm(total=replicates*len(sim_configs), desc="[%] running sims ", ncols=100, unit_scale=True) as pbar:
        for idx, (config, zarr_attrs) in enumerate(zip(msprime_configs, sim_configs)):
            seeds = np.random.randint(1, 2 ** 32, replicates)

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
                    result_list = pool.map(run_sims_specified, seeds)
            else:
                result_list = []
                for seed in seeds:
                    result_list.append(
                        run_ind_sim(
                            seed=seed,
                            msprime_config=config,
                            ploidy=ploidy,
                            blocks=blocks,
                            blocklength=blocklength,
                            comparisons=all_interpop_comparisons,
                        )
                    )

            name = f"parameter_combination_{idx}"
            g = root.create_group(name)
            g.attrs.put(zarr_attrs)
            for idx2, (d, s) in enumerate(zip(result_list, seeds)):
                g.create_dataset(f"replicate_{idx2}", data=d, overwrite=True)
                g[f"replicate_{idx2}"].attrs["seed"] = str(s)
                pbar.update(1)

def make_sim_configs(params, ploidy):
    sample_size_A = params["sample_size_A"]
    sample_size_B = params["sample_size_B"]
    num_samples = sample_size_A + sample_size_B
    C_A = params["C_A"]
    C_B = params["C_B"]
    mutation_rate = params["theta"]
    rec_rate = params["recombination"]

    demographic_events = None
    population_configurations = [
        msprime.PopulationConfiguration(
            sample_size=sample_size_A * ploidy, initial_size=C_A
        ),
        msprime.PopulationConfiguration(
            sample_size=sample_size_B * ploidy, initial_size=C_B
        ),
    ]
    migration_matrix = np.zeros((3, 3))  # migration rate needs to be divided by 4Ne
    if "M_A_B" in params:
        # migration A to B backwards
        migration_matrix[1, 0] = params["M_A_B"]
    if "M_B_A" in params:
        # migration B to A
        migration_matrix[0, 1] = params["M_B_A"]
    if "C_A_B" in params:
        C_AB = params["C_A_B"] if params["C_A_B"] else C_A + C_B
        population_configurations += [
            msprime.PopulationConfiguration(sample_size=0, initial_size=C_AB),
        ]
        # demographic events: specify in the order they occur backwards in time
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
        mutation_rate,
        num_samples,
        rec_rate,
    )


def run_ind_sim(
    seed,
    msprime_config,
    ploidy,
    blocks,
    blocklength,
    comparisons,
):
    (
        population_configurations,
        demographic_events,
        migration_matrix,
        theta,
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
        mutation_rate=theta,
        random_seed=seed, #error was when 3582573439
    )

    """
    #with msprime 1.0 -> finite sites mutations
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
    # with infinite sites = pre-msprime 1.0
    positions = np.array([int(site.position) for site in ts.sites()])
    #print(f"[+] {ts.num_sites} mutation(s) along the simulated sequence")
    new_positions = lib.gimble.fix_pos_array(positions)
    if ts.num_sites>0 and any(p>=total_length for p in new_positions):
        blocklength = new_positions[-1]
        total_length = blocks*blocklength
    genotype_matrix = get_genotypes(ts, ploidy, num_samples)
    sa_genotype_array = allel.GenotypeArray(genotype_matrix)
    # always the same for all pairwise comparisons
    #print("[+] generated genotype matrix")
    # generate all comparisons
    num_comparisons = len(comparisons)
    result = np.zeros((num_comparisons, blocks, blocklength), dtype="int8")
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
        result[idx] = lib.gimble.genotype_to_mutype_array(
            subset_genotype_array, block_sites_variant_bool, block_sites, debug=False
        )
        multiallelic, missing, monomorphic, variation = lib.gimble.block_sites_to_variation_arrays(block_sites)
        print(variation)
    return result


def get_genotypes(ts, ploidy, num_samples):
    shape = (ts.num_mutations, num_samples, ploidy)
    return np.reshape(ts.genotype_matrix(), shape)

def dict_product(d):
    return [dict(zip(d, x)) for x in itertools.product(*d.values())]


def expand_params(d):
    for key, value in d.items():
        if len(value) > 1 and key!="recombination":
            assert len(value) >= 3, "MIN, MAX and STEPSIZE need to be specified"
            sim_range = np.arange(value[0], value[1]+value[2], value[2], dtype=float)
            if len(value)==4:
                if not any(np.isin(sim_range, value[3])):
                    print(f"[-] Specified range for {key} does not contain specified grid center value")  
            d[key] = sim_range

def all_interpopulation_comparisons(popA, popB):
    return list(itertools.product(range(popA), range(popA, popA + popB)))
