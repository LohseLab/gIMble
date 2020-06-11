import pathlib
import oyaml as yaml
import collections
from timeit import default_timer as timer
from docopt import docopt
import sys
from lib.gimble import RunObj
import numpy as np
import msprime


def run_sim(parameterObj):
    ploidy = parameterObj._config["ploidy"]
    params = parameterObj._config["parameters"]
    blocks = parameterObj._config["blocks"]
    blocklength = parameterObj._config["blocklength"]
    windows = parameterObj._config["windows"]
    sample_size_A = params["sample_size_A"]
    inital_size_A = params["C_A"]
    sample_size_B = params["sample_size_B"]
    inital_size_B = params["C_B"]
    total_samples = sample_size_A + sample_size_B
    inital_size_AB = 1
    demographic_events = None
    population_configurations = [
        msprime.PopulationConfiguration(
            sample_size=sample_size_A * ploidy, initial_size=inital_size_A
        ),
        msprime.PopulationConfiguration(
            sample_size=sample_size_B * ploidy, initial_size=inital_size_B
        ),
    ]
    migration_matrix = np.zeros((3, 3))
    if "M_A_B" in params:
        # migration A to B backwards
        migration_matrix[1, 0] = params["M_A_B"]
    if "M_B_A" in params:
        # migration B to A
        migration_matrix[0, 1] = params["M_B_A"]
    if "C_A_B" in params:
        initial_size_AB = params["C_A_B"]
        population_configurations += [
            msprime.PopulationConfiguration(sample_size=0, initial_size=inital_size_AB),
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

    print("[+] running simulations")

    ts = run_ind_sim(
        population_configurations=population_configurations,
        demographic_events=demographic_events,
        migration_matrix=migration_matrix,
        length=blocklength,
        mutation_rate=params["theta"],
        recombination_rate=0.0,
    )
    print("[+] finished simulations")
    genotype_matrix = get_genotypes(ts, ploidy, total_samples)
    print("[+] generated genotype matrix")

def run_ind_sim(
    population_configurations,
    demographic_events,
    migration_matrix,
    length,
    mutation_rate,
    recombination_rate=0.0,
):
    ts = msprime.simulate(
        length=length,
        recombination_rate=recombination_rate,
        population_configurations=population_configurations,
        demographic_events=demographic_events,
        migration_matrix=migration_matrix,
        mutation_rate=mutation_rate,
    )
    return ts

def get_genotypes(ts, ploidy, num_samples):
    shape = (ts.num_mutations, num_samples, ploidy)

    return np.reshape(ts.genotype_matrix(), shape)
        