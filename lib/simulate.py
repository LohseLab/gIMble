import collections
import itertools
from tqdm import tqdm
from sys import exit
import pathlib
import allel
import numpy as np
import pandas as pd
import shutil
import zarr
import os
import logging
import re
import msprime

class ParameterObj(object):
    def __init__(self, args):
        print(args)
        self.module = 'simulate'
        self.pop_ids = args['--pop_ids'].split(',')
        self.samples = [int(element) for element in args['--samples'].split(',')]
        assert len(self.samples)==len(self.pop_ids), "number of populations in samples by population does not match specified number of populations"
        
        self.pop_Ne = [float(element) for element in args['--pop_sizes'].split(',')]
        
        self.Ne_dict = {pop:Ne for pop, Ne in zip(self.pop_ids, self.pop_Ne)}
        
        self.num_blocks = int(args['--blocks'])
        self.block_length = int(args['--block_length'])
        
        self.num_windows = int(args['--windows'])
        self.window_size = int(args['--window_size'])

        self.m_AB, self.m_BA = 0.0, 0.0

        self.migration = args['--migration_string']
        self.migration_rate = float(args['--migration_rate'])
        if self.migration:
            assert self.migration_rate>0, 'Specified migration but no migration rate'
        #this is not pretty
            if self.migration in ['A>B', 'B<A']:
                self.m_BA = self.migration_rate
            else:
                self.m_AB = self.migration_rate   

        self.mutation_rate = float(args['--mutation_rate'])

        #what forms could this take: only "(A,B)"?
        self.join_string = args['--join_string']
        self.ancestral_Ne = float(args['--ancestral_Ne'])

        self.splitT = float(args['--split_time'])
        self.ploidy = int(args['--ploidy'])
        self.recombination_rate = 0.0
        self.parse_parameters()

        #print(self.__dict__)

    def parse_parameters(self):

        if self.join_string and self.ancestral_Ne == 0:
            ancestral = self.join_string[self.join_string.find("(")+1:self.join_string.find(")")].split(',')
            assert len(ancestral[0])>0, 'use brackets for the join_string'
            self.ancestral_Ne = sum(pop_Ne)
            
        elif self.ancestral_Ne > 0 and not self.join_string:
            pass
        else:
            raise AssertionError ("either provide join_string or ancestral_Ne")


    def simulate(self):
        output = get_genotypes(self, run_sim(self))

def run_sim(parameterObj):
    migration_matrix = [[0 for _ in range(3)] for _ in range(3)]
    migration_matrix[0][1] = parameterObj.m_AB
    migration_matrix[1][0] = parameterObj.m_BA
    #M_jk is fraction of ind in j that consists of migrants from k

    population_configurations = [msprime.PopulationConfiguration(sample_size=s*parameterObj.ploidy, initial_size=n) 
        for n, s in zip(parameterObj.pop_Ne,parameterObj.samples)] + [msprime.PopulationConfiguration(sample_size=0, initial_size=parameterObj.ancestral_Ne)]
    
    #demographic events: specify in the order they occur backwards in time
    demographic_events = [
    msprime.MassMigration(time=parameterObj.splitT, source=0, destination=2, proportion=1.0),
    msprime.MassMigration(time=parameterObj.splitT, source=1, destination=2, proportion=1.0),
    msprime.MigrationRateChange(parameterObj.splitT, 0),
    ]

    replicates = msprime.simulate(
        num_replicates=parameterObj.num_blocks,
        length = parameterObj.block_length, 
        recombination_rate = parameterObj.recombination_rate,
        population_configurations = population_configurations,
        demographic_events = demographic_events,
        migration_matrix = migration_matrix,
        mutation_rate = parameterObj.mutation_rate)
    
    return replicates

def get_genotypes(parameterObj, replicates):
    num_replicates = parameterObj.num_blocks
    ploidy = parameterObj.ploidy
    num_samples = sum(parameterObj.samples)

    temp = []
    result = []
    num_mutations = 0
    for i, ts in enumerate(replicates):
        shape = (ts.num_mutations, num_samples, ploidy)
        temp.append(np.reshape(ts.genotype_matrix(), shape))
        num_mutations = max(ts.num_mutations, num_mutations)
    for genotype_matrix in temp:
        padding = np.zeros((num_mutations, num_samples, ploidy))
        padding[:genotype_matrix.shape[0],:,:] = genotype_matrix
        result.append(padding)
    #np.save('blocks', result)
    print("saved blocks.npy file")
    return result
