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

def run_sims(sim_configs, global_info, all_interpop_comparisons, chunks=1, threads=1, store=None):
	"""
	Arguments:
		sim_configs {list} -- [containing dicts with keys: Ne_x, me_x_x, T, recombination]
		global_info {dict} -- [keys: mu, ploidy, sample_pop_ids, sample_pop_sizes, blocklength,
									blocks, k_max, chunks, replicates, reference_pop]
		all_interpop_comparisons {list} -- 
	
	Keyword Arguments:
		chunks {int} -- [number of chunks in which sequence is being split up prior to simulating] (default: {1})
		threads {int} -- [number of cores used to parallelize over replicates] (default: {1})
		store {[object]} -- [zarr store] (default: {None})
	
	Returns:
		[type] -- [description]
	"""
	msprime_configs = (make_sim_configs(config, global_info) for config in sim_configs)
	all_results=[]
	for idx, (single_config, zarr_attrs) in enumerate(tqdm(zip(msprime_configs, sim_configs),desc='Overall simulation progress',ncols=100, unit_scale=True, total=len(sim_configs))):
		seeds = np.random.randint(1, 2 ** 32, global_info["replicates"])
		if threads > 1:
			result_list = run_sim_parallel(single_config, global_info, seeds, all_interpop_comparisons, idx, threads)
		else:
			result_list = run_sim_serial(single_config, global_info, seeds, all_interpop_comparisons, idx)
		result_list = _combine_chunks(result_list, chunks)
		if store is not None:
			name = f"parameter_combination_{idx}"
			store.create_dataset(name, data=result_list, overwrite=True)
			store.attrs.put(zarr_attrs)
			store.attrs['seeds']=tuple([int(s) for s in seeds])
		else:
			all_results.append(result_list)
	return np.array(all_results, dtype=np.int64) #check do we need int64?

def simulate_parameterObj(sim_configs, parameterObj, gimbleStore):
	global_info = compile_global_info(parameterObj)
	all_interpop_comparisons = all_interpopulation_comparisons(*global_info['sample_pop_sizes'])
	print(f'[+] simulating {int(global_info["replicates"]//global_info["chunks"])} replicate(s) of {int(global_info["blocks"]*global_info["chunks"])} block(s) for {len(sim_configs)} parameter combinations')
	if parameterObj.label:
		group_name=parameterObj.label
	else:
		run_count = gimbleStore._return_group_last_integer('sims')
		group_name = f"run_{run_count}"
	gimbleStore.data.require_group(f'sims/{group_name}')
	gimbleStore.data[f'sims/{group_name}'].attrs.put(global_info)
	gimbleStore.data[f'sims/{group_name}'].attrs['fixed_param_grid'] = parameterObj.fixed_param_grid
	run_sims(sim_configs, global_info, all_interpop_comparisons, global_info["chunks"], parameterObj.threads, gimbleStore.data[f'sims/{group_name}'])

def compile_global_info(parameterObj):
	#global_info_list =['mu', 'ploidy', 'sample_pop_ids', 'sample_pop_sizes', 'blocklength', 
	#                      'blocks','k_max', 'chunks', 'replicates']
	global_info = parameterObj.config['simulations'].copy()    
	global_info['mu'] = parameterObj.config['mu']['mu']
	global_info['blocklength'] = parameterObj.config['mu']['blocklength']
	global_info['k_max'] = list(parameterObj.config['k_max'].values())
	global_info['sample_pop_ids'] = sorted(parameterObj.config['populations']['sample_pop_ids'])
	global_info['sample_pop_sizes'] = [global_info[f'sample_size_{pop_id}'] for pop_id in global_info['sample_pop_ids']]
	global_info['reference_pop'] = parameterObj.config['populations']['reference_pop']
	#check if number of chunks is valid
	if global_info['chunks']>1:
		if global_info['blocks']==1:
			print(f"[-] Can't split 1 block into {global_info['chunks']} chunks. Simulation will continue without chunking.")
			global_info['chunks']=1
		else:
			global_info['blocks']//=global_info['chunks']
			global_info['replicates']*=global_info['chunks']
	return global_info

def run_sim_parallel(config, global_info, seeds, all_interpop_comparisons, idx, threads):

	with multiprocessing.Pool(processes=threads) as pool:
		run_sims_specified = partial(
		run_ind_sim,
		msprime_config=config,
		ploidy=global_info["ploidy"],
		blocks=global_info["blocks"],
		blocklength=global_info["blocklength"],
		comparisons=all_interpop_comparisons,
		k_max=global_info["k_max"]
	)
	
		result_list = list(tqdm(pool.imap(run_sims_specified, seeds),desc=f'running parameter combination {idx}',ncols=100, unit_scale=True, total=len(seeds)))
	return result_list
	   
def run_sim_serial(config, global_info, seeds, all_interpop_comparisons, idx):
	result_list = []
	for seed in tqdm(seeds,desc=f'running parameter combination {idx}',ncols=100, unit_scale=True):
		result_list.append(
			run_ind_sim(
				seed=seed,
				msprime_config=config,
				ploidy=global_info["ploidy"],
				blocks=global_info["blocks"],
				blocklength=global_info["blocklength"],
				comparisons=all_interpop_comparisons,
				k_max=global_info["k_max"]
			)
		)
	return result_list
			
def make_sim_configs(params, global_info):
	A, B = global_info["sample_pop_ids"]
	sample_size_A, sample_size_B = global_info['sample_pop_sizes']
	num_samples = sum(global_info['sample_pop_sizes'])
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
	comparisons,
	k_max
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
		random_seed=seed,
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
	result = np.zeros((num_comparisons, blocks, len(k_max)), dtype="int64") #get number of mutypes
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
	#flatten
	result = result.reshape(-1, result.shape[-1])
	#result = np.hstack(result).reshape(num_comparisons*blocks, len(k_max))
	#apply get_bsfs to this
	# count mutuples (clipping at k_max, if supplied)
	max_k = np.array(k_max) + 1 if k_max else None
	mutuples, counts = np.unique(np.clip(result, 0, max_k), return_counts=True, axis=0)
	# define out based on max values for each column
	out = np.zeros(tuple(max_k + 1), np.int64)
	# assign values
	out[tuple(mutuples.T)] = counts
	return out

def process_single_simulation_run():
	pass

def _combine_chunks(result_list, chunks):
	if chunks>1:
		return np.array([np.add.reduce(m) for m in np.split(result_list,list(range(0,len(result_list),chunks))[1:])])
	else:
		return np.array(result_list)

def get_genotypes(ts, ploidy, num_samples):
	if ts.num_mutations == 0:
		return np.zeros((1,num_samples, ploidy), dtype='int8')
	shape = (ts.num_mutations, num_samples, ploidy)
	return np.reshape(ts.genotype_matrix(), shape)

def all_interpopulation_comparisons(*popsizes):
	popA, popB, *rest = popsizes
	if len(rest)>0:
		raise ValueError("More than 2 population sizes were provided to simulate. We cannot cope with that just yet.")
	return list(itertools.product(range(popA), range(popA, popA + popB)))