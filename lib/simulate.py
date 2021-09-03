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
import pandas as pd
from functools import partial
import collections

def run_sims(demographies, recombination_rates, global_info, replicates, seeds, threads=1, discrete=True, disable_tqdm=False): 
	samples = {
               'A': global_info['simulate']['sample_size_A'], 
               'B': global_info['simulate']['sample_size_B']
               }
	mutation_model = 'infinite_alleles'
	if isinstance(recombination_rates, float):
		recombination_rates = itertools.repeat(recombination_rates)
	for idx, (demography,rec_rate, seed) in enumerate(tqdm(
		zip(demographies,recombination_rates, seeds),
		desc='Overall simulation progress',
		ncols=100, 
		unit_scale=True, 
		total=global_info['parameters_grid_points'], 
		disable=disable_tqdm
		)
	):
		if threads > 1:
			result_list = run_sim_parallel(
				seed, 
				rec_rate, 
				samples, 
				demography, 
				global_info['simulate']['ploidy'], 
				global_info['simulate']['block_length']*global_info['simulate']['blocks'], 
				global_info['simulate']['comparisons'], 
				global_info['max_k'],
				global_info['mu']['mu'],
				mutation_model,
				discrete,
				global_info['simulate']['blocks'], 
				disable_tqdm,
				idx,
				threads
				)
		else:
			result_list = run_sim_serial(
				seed, 
				rec_rate, 
				samples, 
				demography, 
				global_info['simulate']['ploidy'], 
				global_info['simulate']['block_length']*global_info['simulate']['blocks'], 
				global_info['simulate']['comparisons'], 
				global_info['max_k'], 
				global_info['mu']['mu'],
				mutation_model,
				discrete,
				global_info['simulate']['blocks'], 
				disable_tqdm,
				idx
				)
		result_list = _combine_chunks(result_list, global_info['simulate']['chunks'])
		yield np.array(result_list)

def sim_worker(seed, recombination_rate, samples, demography, ploidy, sequence_length, blocks, comparisons, max_k, mutation_rate, mutation_model, discrete):
	ancestry_seed, mutation_seed = seed
	#run simulation:
	ts = msprime.sim_ancestry(
		samples=samples, 
    	demography = demography, 
    	ploidy=ploidy,
    	sequence_length=sequence_length,
    	discrete_genome=discrete,
    	recombination_rate=recombination_rate,
    	random_seed=ancestry_seed
    	)
	ts = msprime.sim_mutations(
		ts, 
		rate=mutation_rate, 
		random_seed=mutation_seed, 
		discrete_genome=discrete, 
		model=mutation_model
		)
	num_samples = sum(samples.values())
	#make bsfs:
	if discrete:
		positions = np.array([site.position for site in ts.sites()])
	else:
		positions, sequence_length = infinite_sites(ts, blocks, sequence_length)
	genotype_matrix = get_genotypes(ts, ploidy, num_samples)
	bsfs = generate_bsfs(genotype_matrix, positions, comparisons, max_k, blocks, sequence_length)
	return bsfs

def run_sim_parallel(seeds, recombination_rate, samples, demography, ploidy, sequence_length, comparisons, max_k, mutation_rate, mutation_model, discrete, blocks, disable_tqdm, idx, threads):
	with multiprocessing.Pool(processes=threads) as pool:
		run_sims_specified = partial(
		sim_worker,
		recombination_rate=recombination_rate,
		samples=samples,
		demography=demography,
		ploidy=ploidy,
		sequence_length=sequence_length,
		blocks=blocks,
		comparisons=comparisons,
		max_k=max_k,
		mutation_rate=mutation_rate,
		mutation_model=mutation_model,
		discrete=discrete
	)
	
		result_list = list(tqdm(pool.imap(run_sims_specified, seeds),desc=f'running parameter combination {idx}',ncols=100, unit_scale=True, total=len(seeds),disable=disable_tqdm))
	return result_list
	   
def run_sim_serial(seeds, recombination_rate, samples, demography, ploidy, sequence_length, comparisons, max_k, mutation_rate, mutation_model, discrete, blocks, disable_tqdm, idx):
	# [SIMULATION]
	shape = np.insert(max_k+2, 0, seeds.shape[0]) 
	result_list = np.zeros(shape, dtype=lib.gimble._return_np_type(blocks))
	for sub_idx, seed in enumerate(tqdm(seeds, desc=f'running parameter combination {idx}',ncols=100, unit_scale=True, disable=disable_tqdm)):
		result_list[sub_idx] = sim_worker(
				seed=seed,
				recombination_rate=recombination_rate,
				samples=samples,
				demography=demography,
				ploidy=ploidy,
				sequence_length=sequence_length,
				blocks=blocks,
				comparisons=comparisons,
				max_k=max_k,
				mutation_rate=mutation_rate,
				mutation_model=mutation_model,
				discrete=discrete
			)
	return result_list

def generate_bsfs(genotype_matrix, positions, comparisons, max_k, blocks, total_length):
	sa_genotype_array = allel.GenotypeArray(genotype_matrix)
	num_comparisons = len(comparisons)
	result = np.zeros((num_comparisons, blocks, len(max_k)), dtype=np.int64)
	# generate all comparisons

	for idx, pair in enumerate(comparisons):
		block_sites = np.arange(total_length).reshape(blocks, -1)
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
	return (new_positions, total_length)

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
	return (msprime.Demography.from_demes(graph) for graph in lib.gimble.config_to_demes_graph(config))