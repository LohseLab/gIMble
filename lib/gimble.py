import itertools
from tqdm import tqdm
from lib.functions import plot_mutuple_barchart
import allel
import numpy as np
import pandas as pd
import shutil
import zarr
import os
import string
import logging
import collections
import sys
import warnings
import pathlib
import matplotlib.pyplot as plt
# np.set_printoptions(threshold=sys.maxsize)

# Formatting functions

'''
ZARR
 - seq
    - meta 
    - chroms

 - sim
    - meta
        - 
        - simulaton = {'A' : [1], 'B': [2,3,4,5,6,7,8, ... 100] }
    - 1
        - GTS
    - 2
    - 3
    - 4
    - 5

block_starts, seq = "seq/blocks/%s/starts" % (seq)
'''

def format_bases(bases):
    return "%s b" % format(bases, ',d')

def format_percentage(fraction, precision=2):
    return "{:.{}%}".format(fraction, precision)

def format_fraction(fraction, precision=2):
    return "{:.{}}".format(fraction, precision)

def format_count(count):
    return "%s" % str(format(count, ',d'))

def get_n50_from_lengths(lengths):
    length_sorted = sorted(lengths, reverse=True)
    cum_sum = np.cumsum(length_sorted)
    half = int(sum(lengths)/2)
    cum_sum_2 =min(cum_sum[cum_sum >= half])
    n50_idx = np.where(cum_sum == cum_sum_2)
    return length_sorted[int(n50_idx[0][0])]
    
def fix_pos_array(pos_array):
    '''
    De-duplicates array by shifting values forward until there aren't any collisions
    '''
    # get boolean array for first and subsequent duplicates (True) (required sorted) 
    idxs = np.insert((np.diff(pos_array)==0).astype(np.bool), 0, False)
    if np.any(idxs): 
        # if there are duplicates, get new values by incrementing by one
        new_values = pos_array[idxs]+1
        # get non-duplicate values
        uniq_values = pos_array[~idxs]
        # insert new_values in non-duplicated values (required sorted)
        new_idxs = np.searchsorted(uniq_values, new_values)
        # recursive call
        return fix_pos_array(np.sort(np.insert(uniq_values, new_idxs, new_values)))
    # if there are no duplicated values
    return pos_array

############## Only needed once we have demand for multidimensional pairing function
# def multidimensional_box_pairing(lengths: List[int], indexes: List[int]) -> int:
#     n = len(lengths)
#     index = 0
#     dimension_product = 1
#     for dimension in reversed(range(n)):
#         index += indexes[dimension] * dimension_product
#         dimension_product *= lengths[dimension]
#     return index
# def multidimensional_szudzik_pairing(*indexes: int) -> int:
#     n = len(indexes)
#     if n == 0:
#         return 0
#     shell = max(indexes)
#     def recursive_index(dim: int):
#         slice_dims = n - dim - 1
#         subshell_count = (shell + 1) ** slice_dims - shell ** slice_dims
#         index_i = indexes[dim]
#         if index_i == shell:
#             return subshell_count * shell + multidimensional_box_pairing([shell + 1] * slice_dims, indexes[dim + 1:])
#         else:
#             return subshell_count * index_i + recursive_index(dim + 1)
#     return shell ** n + recursive_index(0)

def szudzik_pairing(folded_minor_allele_counts):
    # adapted from: https://drhagen.com/blog/superior-pairing-function/
    if isinstance(folded_minor_allele_counts, np.ndarray):
        # assumes folded_minor_allele_counts is array with shape (n,2)
        return np.where(
            (folded_minor_allele_counts[:,0] >= folded_minor_allele_counts[:,1]),
            np.square(folded_minor_allele_counts[:,0]) + folded_minor_allele_counts[:,0] + folded_minor_allele_counts[:,1],
            folded_minor_allele_counts[:,0] + np.square(folded_minor_allele_counts[:,1])
            )
    elif isinstance(folded_minor_allele_counts, tuple):
        # assumes folded_minor_allele_counts is tuple of size 2
        a, b = folded_minor_allele_counts
        if a >= b:
            return (a**2) + a + b 
        return a + (b**2)
    else:
        pass

def get_coverage_counts(coverages, idxs, num_blocks):
    # not used at the moment
    num_sample_sets = idxs[-1] + 1 # this is correct, don't worry about it ...
    temp = coverages + (num_sample_sets * np.arange(coverages.shape[0]))[:, None]
    blocks_per_sample_set = np.bincount(temp.ravel(), minlength=(num_sample_sets * coverages.shape[0])).reshape(-1, num_sample_sets)
    # remove columns/sample_sets that only contain zeroes and normalise
    return blocks_per_sample_set[:, ~(blocks_per_sample_set == 0).all(0)] / (num_blocks / len(idxs)) / num_blocks 

def block_sites_to_variation_arrays(block_sites, max_type_count=4):
    block_count = block_sites.shape[0]
    max_type_count += 3 # ideally this should be the maximum amount of mutypes + 2 + 1 
    temp_sites = block_sites + (max_type_count * np.arange(block_count).reshape(block_count, 1))
    # return multiallelic, missing, monomorphic, variation
    return np.hsplit(np.bincount(temp_sites.ravel(), minlength=(block_count * max_type_count)).reshape(-1, max_type_count), [1, 2, 3])

def calculate_popgen_from_array(mutype_array, sites):
    # print('# Mutypes: 0=MULTI, 1=MISS, 2=MONO, 3=HetB, 4=HetA, 5=HetAB, 6=Fixed')
    pi_1 = float("%.8f" % np.divide(np.sum(mutype_array[:, 1]) + np.sum(mutype_array[:, 2]), sites)) # average heterozygosity
    pi_2 = float("%.8f" % np.divide(np.sum(mutype_array[:, 0]) + np.sum(mutype_array[:, 2]), sites)) # average heterozygosity
    d_xy = float("%.8f" % np.divide(np.divide(np.sum(mutype_array[:, 0]) + np.sum(mutype_array[:, 1]) + np.sum(mutype_array[:, 2]), 2.0) + np.sum(mutype_array[:, 3]), sites))
    mean_pi = (pi_1 + pi_2) / 2.0
    total_pi = (d_xy + mean_pi) / 2.0 # special case of pairwise Fst
    f_st = np.nan
    if (total_pi):
        f_st = float("%.8f" % ((total_pi - mean_pi) / total_pi)) # special case of pairwise Fst
    fgv = len(mutype_array[(mutype_array[:, 2] > 0) & (mutype_array[:, 3] > 0)])
    return (pi_1, pi_2, d_xy, f_st, fgv)

def genotype_to_mutype_array(sa_genotype_array, idx_block_sites_in_pos, block_sites, debug=True):
    if debug == True:
        print("# Blocksites", block_sites)
    np_genotype_array = np.array(sa_genotype_array)
    #print('np_genotype_array', np_genotype_array)
    #print('block_sites', block_sites.shape)
    #print(block_sites)
    #print('idx_block_sites_in_pos', idx_block_sites_in_pos.shape)
    #print(idx_block_sites_in_pos)
    #print('sa_genotype_array', sa_genotype_array.shape)
    #print(sa_genotype_array)
    #print('np_genotype_array', np_genotype_array.shape)
    #print(np_genotype_array)
    np_allele_count_array = np.ma.masked_equal(sa_genotype_array.count_alleles(), 0, copy=False)    
    #print('np_allele_count_array', np_allele_count_array)
    allele_map = np.ones((np_allele_count_array.shape), dtype='int8') * np.arange(np_allele_count_array.shape[-1])
    idx_max_global_allele_count = np.nanargmax(np_allele_count_array, axis=1)
    idx_min_global_allele_count = np.nanargmin(np_allele_count_array, axis=1)
    has_major_allele = (idx_max_global_allele_count != idx_min_global_allele_count)
    idx_min_prime_allele = np.amin(np_genotype_array[:,0], axis=1)
    idx_min_global_allele = np.amin(np.amin(np_genotype_array, axis=1), axis=1)
    idx_max_global_allele = np.amax(np.amax(np_genotype_array, axis=1), axis=1)
    idx_major_allele = np.where(
        has_major_allele, 
        idx_max_global_allele_count, 
        idx_min_prime_allele)
    idx_minor_allele = np.where(
        has_major_allele, 
        idx_min_global_allele_count, 
        np.where((
            idx_min_global_allele == idx_min_prime_allele),
            np.max((idx_min_global_allele, idx_max_global_allele), axis=0), 
            np.min((idx_min_global_allele, idx_max_global_allele), axis=0)))
    # for each genotype (np.arange(allele_map.shape[0])), set minor allele to 1 (1st do minor, so that overwritten if monomorphic)
    allele_map[np.arange(allele_map.shape[0]), idx_minor_allele] = 1 
    # for each genotype (np.arange(allele_map.shape[0])), set major allele to 0
    allele_map[np.arange(allele_map.shape[0]), idx_major_allele] = 0
    folded_minor_allele_counts = sa_genotype_array.map_alleles(allele_map).to_n_alt(fill=-1)
    folded_minor_allele_counts[np.any(sa_genotype_array.is_missing(), axis=1)] = np.ones(2) * -1        # -1, -1 for missing => -1
    folded_minor_allele_counts[(np_allele_count_array.count(axis=1) > 2)] = np.ones(2) * (-1, -2)       # -1, -2 for multiallelic => -2
    block_sites_pos = block_sites.flatten()
    block_sites[idx_block_sites_in_pos] = szudzik_pairing(folded_minor_allele_counts) + 2               # add 2 so that not negative for bincount
    block_sites[~idx_block_sites_in_pos] = 2                                                            # monomorphic = 2 (0 = multiallelic, 1 = missing)
    if debug == True:
        pos_df = pd.DataFrame(block_sites_pos[idx_block_sites_in_pos.flatten()], dtype='int', columns=['pos'])
        genotypes_df = pd.DataFrame(np_genotype_array.reshape(np_genotype_array.shape[0], 4), dtype='i4', columns=['a1', 'a2', 'b1', 'b2'])        
        block_sites_df = pos_df.join(genotypes_df)
        folded_minor_allele_count_df = pd.DataFrame(folded_minor_allele_counts, dtype='int8', columns=['fmAC_a', 'fmAC_b'])
        block_sites_df = block_sites_df.join(folded_minor_allele_count_df)
        variants = pd.DataFrame(block_sites[idx_block_sites_in_pos], dtype='int', columns=['SVar'])
        block_sites_df = block_sites_df.join(variants)
        print('# Mutypes: 0=MULTI, 1=MISS, 2=MONO, 3=HetB, 4=HetA, 5=HetAB, 6=Fixed')
        print(block_sites_df)
    return block_sites

def create_ranges(aranges):
    # does the transformation form 0-based (BED) to 1-based (VCF) coordinate system 
    # https://stackoverflow.com/a/47126435
    l = (aranges[:, 1] - aranges[:, 0])
    clens = l.cumsum()
    ids = np.ones(clens[-1], dtype=int)
    ids[0] = aranges[0, 0]
    ids[clens[:-1]] = aranges[1:, 0] - aranges[:-1, 1] + 1
    return ids.cumsum()

def cut_windows(mutype_array, idxs, start_array, end_array, sample_set_array, num_blocks=10, num_steps=3):
    coordinate_sorted_idx = np.argsort(end_array)
    mutype_array_sorted = mutype_array.take(coordinate_sorted_idx, axis=0)
    window_idxs = np.arange(mutype_array_sorted.shape[0] - num_blocks + 1)[::num_steps, None] + np.arange(num_blocks)
    window_mutypes = mutype_array_sorted.take(window_idxs, axis=0)
    #print("window_mutypes", window_mutypes[:])
    block_starts = start_array.take(coordinate_sorted_idx, axis=0).take(window_idxs, axis=0)
    #print("block_starts", block_starts[:])
    window_starts = np.min(start_array.take(coordinate_sorted_idx, axis=0).take(window_idxs, axis=0), axis=1).T
    #print("window_starts", window_starts[:])
    block_ends = end_array.take(coordinate_sorted_idx, axis=0).take(window_idxs, axis=0)#[:,-1]
    window_ends = np.max(end_array.take(coordinate_sorted_idx, axis=0).take(window_idxs, axis=0), axis=1).T
    window_midpoints = (block_starts / 2) + (block_ends / 2)
    window_pos_mean = np.mean(window_midpoints, axis=1).T
    window_pos_median = np.median(window_midpoints, axis=1).T
    # coverages = get_coverage_counts(sample_set_array.take(window_idxs, axis=0), idxs, num_blocks)
    return window_mutypes, window_starts, window_ends, window_pos_mean, window_pos_median

def cut_blocks(interval_starts, interval_ends, block_length, block_span, block_gap_run):
    sites = create_ranges(np.array((interval_starts, interval_ends)).T)
    block_sites = np.concatenate([
        x[:block_length * (x.shape[0] // block_length)].reshape(-1, block_length) 
            for x in np.split(sites, np.where(np.diff(sites) > block_gap_run)[0] + 1)
        ]) 
    block_span_valid_mask = (((block_sites[:, -1] - block_sites[:, 0]) + 1) <= block_span)
    return block_sites[block_span_valid_mask] 

def create_store(parameterObj):
    store = Store()
    store._from_input(parameterObj)
    return store

def load_store(parameterObj):
    store = Store() 
    store._from_zarr(parameterObj)
    return store

class RunObj(object):
    '''Superclass for ParameterObjs'''
    def __init__(self, params):
        self._PATH = params['path']
        self._VERSION = params['version']
        self._MODULE = params['module']
        self._CWD = params['cwd']

    def __repr__(self):
        return("[+] VER := %s\n[+] CWD := %s\n[+] CMD := %s\n" % (
            self._VERSION, 
            self._CWD, 
            self._get_cmd()
            ))

    def _get_cmd(self):
        return "%s/gimble %s %s" % (
            self._PATH, 
            self._MODULE, 
            "".join(["--%s " % " ".join((k, str(v))) for k,v in self.__dict__.items() if not k.startswith("_")]))

    def _get_path(self, infile, path=False):
        if infile is None:
            return None
        path = pathlib.Path(infile).resolve()
        if not path.exists():
            sys.exit("[X] File not found: %r" % str(infile))
        if path:
            return path
        return str(path)

    def _get_int(self, string):
        try:
            return(int(string))
        except TypeError():
            sys.exit("[X] %r can't be converted to interger." % string)

    def _get_float(self, string):
        try:
            return(float(string))
        except TypeError():
            sys.exit("[X] %r can't be converted to float." % string)


SCHEMA = {
    'main';
        'attrs':
            'module': 'command'
    }

seqs,
    meta

sims,
    meta
META_SCHEMA = {
            'meta' : {
                'files': {'vcf_f': None, 'sample_f': None, 'genome_f': None, 'bed_f': None}, 
                # genome_file
                'sequences': {'sequence_ids': [], 'sequence_lengths': []},
                # sample_file
                'samples': {'sample_ids' = [], 'sample_ids_vcf' = [], 'sample_ids_bed' = [], 'population_ids' = [], 'sample_sets' = [], 'sample_sets_cartesian' = []},
                # bed_file
                'bed': {'intervals': None, 'span': None}
                # blocks
                'blocks': {'block_length' : None, 'block_span' : None, 'block_max_missing' : None, 'block_max_multiallelic' : None},
                # windows
                'windows': {'window_size': None, 'window_step': None}
                }
            }

# for each sequence
#'variants/pos'
#'variants/gt_matrix'
#'intervals/starts'
#'intervals/ends'
#'intervals/interval_matrix'
#'blocks/starts'
#'blocks/ends'
#'blocks/variation'
#'blocks/missing'
#'blocks/multiallelic'
#'windows/starts'
#'windows/ends'
#'windows/variation'
#'windows/pos_mean'
#'windows/pos_median'

class Store(object):
    def __init__(self, zarr_dir, create=False, overwrite=False):
        self.prefix = self._get_zarr_prefix(zarr_dir)
        self.path = self._get_zarr_path(zarr_dir)
        self.data = self._get_store_data(create, overwrite)
    
    def tree(self):
        return self.data.tree()
    
    def log_stage(self, parameterObj):
        '''only remembers last command per module'''
        self.data.attrs[parameterObj._MODULE] = parameterObj._get_cmd()
    
    def get_stage(self, stage):
        '''can only get last command per module'''
        return self.data.attrs[stage]
    
    def _get_store_data(self, create, overwrite):
        if create:
            if os.path.isdir(self.path):
                logging.info("[-] Gimble store %r already exists." % self.path)
                if not overwrite:
                    logging.info("[X] Please specify '-f' to overwrite.")
                    sys.exit(1)
                logging.info("[+] Deleting existing ZARR store %r" % self.path)
                shutil.rmtree(self.path)
            logging.info("[+] Creating ZARR store %r" % self.path)
            return zarr.open(str(self.path), mode='w')
        logging.info("[+] Loading ZARR store %r" % self.path)
        return zarr.open(str(self.path), mode='r+')
    
    def _get_zarr_prefix(self, zarr_dir):
        '''used as prefix for output files'''
        if isinstance(zarr_dir, pathlib.Path):
            # remove suffix 
            # ['blocks', 'windows', 'inference'] 
            return str(zarr_dir.resolve().stem)
        # ['setup']
        return zarr_dir
    
    def _get_zarr_path(self, zarr_dir):
        if isinstance(zarr_dir, pathlib.Path):
            # ['blocks', 'windows', 'inference'] 
            return str(zarr_dir) 
        # add suffix 
        # ['setup']
        return "%s.z" % zarr_dir
    
    def setup(self, parameterObj):
        logging.info("[#] Initialising store...")
        self._init_meta(overwrite=True)
        logging.info("[#] Processing genome file %r..." % parameterObj.genome_file)
        self._parse_genome_file(parameterObj.genome_file)
        logging.info("[#] Processing sample file %r..." % parameterObj.sample_file)
        self._parse_sample_file(parameterObj.sample_file, 2)
        logging.info("[#] Processing vcf file %r..." % parameterObj.vcf_file)
        self._parse_vcf_file(parameterObj.vcf_file)
        logging.info("[#] Processing vcf file %r..." % parameterObj.bed_file)
        self._parse_bed_file(parameterObj.bed_file)
        self.add_stage(parameterObj)
    
    def _init_meta(self, overwrite=False):
        # seqs
        defaults_by_group = {
            'seqs/meta/files': {'vcf_f': '', 'sample_f': '', 'genome_f': '', 'bed_f': ''},
            'seqs/meta/sequences' : {'sequence_ids': [], 'sequence_lengths': [], 'n50' = 0},
            'seqs/meta/samples' : {'samples' : [], 'populations' : [], 'sample_sets' : [], 'sample_sets_cartesian' : []},
            'seqs/meta/variants/': {
                'variant_counts': [], 'vcf_samples': [], 'gt_idx_by_sample': {}, 
                'count_hom_ref_by_sample': {}, 'count_hom_alt_by_sample': {}, 'count_het_by_sample': {}, 'count_missing_by_sample': {}}
            'seqs/meta/intervals' : {'count': 0, 'span': 0, 'bed_samples': []},
            'seqs/meta/blocks' : {'block_length' : 0, 'block_span' : 0, 'block_max_missing' : 0, 'block_max_multiallelic' : 0, 'span_in_blocks': 0},
            'seqs/meta/windows' : {'window_size': 0, 'window_step': 0, 'window_count': 0, 'span_in_windows': 0}
        }
        for group, default in defaults_by_group.items():
            self.data.require_group(group, overwrite=overwrite)
            self.data[group].attrs.put(default)    
        # sims
        # ...

    def _parse_genome_file(self, genome_file):
        df = pd.read_csv(
            genome_file, sep="\t", usecols=[0,1], names=['sequence_id', 'sequence_length'], 
            header=None, dtype={'sequence_id': str, 'sequence_length': 'Int64'})
        print("[TBT] Test for bad TSV")
        print(df) 
        if df.isnull().values.any():
            sys.exit("[X] Bad TSV file %r." % genome_file)
        self.data.attrs['seqs/meta/sequences/']['ids'] = df['sequence_id'].to_list()
        self.data.attrs['seqs/meta/sequences/']['lengths'] = df['sequence_length'].to_list()
        self.data.attrs['seqs/meta/sequences/']['n50'] = get_n50_from_lengths(self.data['seqs/meta/sequences/lengths'])
        for sequence_id in self.data.attrs['sequence_ids']:
            self.data.create_group('seqs/%s/' % sequence_id)
        self.data.attrs['seqs/meta/files/']['genome_f'] = str(genome_file)

    def _parse_sample_file(self, sample_file, pairedness):
        df = pd.read_csv(
            sample_file, sep=",", usecols=[0,1], names=['samples', 'populations'],
            header=None, dtype={'samples': str, 'populations': str})
        print("[TBT] Test for bad CSV")
        print(df)
        if df.isnull().values.any():
            sys.exit("[X] Bad CSV file %r." % sample_file)
        self.data.attrs['seqs/meta/samples/']['samples'] = df['sequence_id'].to_list()
        self.data.attrs['seqs/meta/samples/']['populations'] = df['population_ids'].to_list()
        self.data.attrs['seqs/meta/samples/']['population_ids'] = sorted(set(df['population_id'].to_list()))
        population_by_sample = {sample: population for sample, population in zip(
            self.data.attrs['seqs/meta/samples/']['samples'], 
            self.data.attrs['seqs/meta/samples/']['populations'])}
        self.data.attrs['seqs/meta/samples/']['population_by_sample'] = population_by_sample
        self.data.attrs['seqs/meta/samples/']['sample_sets'] = [
            tuple(sorted(x, key=(population_by_sample.get if population_by_sample[x[0]] != population_by_sample[x[1]] else None))) 
            for x in itertools.combinations(population_by_sample.keys(), 2)]
        self.data.attrs['seqs/meta/samples/']['sample_sets_cartesian'] = [
            True if 
                self.data.attrs['seqs/meta/samples/population_by_sample'][sample_set[0]] == self.data.attrs['seqs/meta/samples/population_by_sample'][sample_set[1]] 
            else False for sample_set in self.data.attrs['seqs/meta/samples/sample_sets']]
        self.data.attrs['seqs/meta/files/']['sample_f'] = str(sample_file)

    def _parse_vcf_file(self, vcf_file):
        sequences = self.data.attrs['seqs/meta/sequences/']['id']
        lengths = self.data.attrs['seqs/meta/sequences/']['lengths']
        samples = self.data.attrs['seqs/meta/samples/']['samples']
        variant_counts = []
        count_hom_ref_by_sample = {}
        count_hom_alt_by_sample = {}
        count_het_by_sample = {}
        count_missing_by_sample = {}
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            gt_key, pos_key, sample_key = 'calldata/GT', 'variants/POS', 'samples'
            vcf_samples = allel.read_vcf(vcf_f, fields=[sample_key])[sample_key]
            query_samples = [vcf_sample for vcf_sample in vcf_samples if vcf_sample in set(samples)] # order as they appear genotypes
            self.data.attrs['seqs/meta/variants/']['gt_samples'] = query_samples
            self.data.attrs['seqs/meta/variants/']['gt_idx_by_sample'] = {query_sample: idx for idx, query_sample in enumerate(query_samples)}
            for idx, (sequence, length) in tqdm(enumerate(zip(sequences, lengths)), total=len(sequences), desc="[%] Reading variants...", ncols=100):
                vcf_data = allel.read_vcf(vcf_file, region=seq_id, samples=samples_query, fields=[gt_key, pos_key])
                # genotypes
                gt_matrix = vcf_data[gt_key]
                self.data.create_dataset("seqs/%s/variants/gt_matrix'" % seq_id, data=gt_matrix)
                variant_counts.append(gt_matrix.shape[0])
                sa_genotype_matrix = allel.GenotypeArray(gt_matrix)
                count_array_hom_ref = sa_genotype_matrix.count_hom_ref(axis=0)
                count_array_hom_alt = sa_genotype_matrix.count_hom_alt(axis=0)
                count_array_het = sa_genotype_matrix.count_het(axis=0)
                count_array_missing = sa_genotype_matrix.count_missing(axis=0)
                for sample, vcf_idx in self.data.attrs['seqs/meta/variants/']['gt_idx_by_sample'].items():
                    count_hom_ref_by_sample[sample] = hom_ref_by_sample.get(sample, 0) + count_array_hom_ref[idx]
                    count_hom_alt_by_sample[sample] = hom_alt_by_sample.get(sample, 0) + count_array_hom_alt[idx]
                    count_het_by_sample[sample] = het_by_sample.get(sample, 0) + count_array_het[idx]
                    count_missing_by_sample[sample] = missing_by_sample.get(sample, 0) + count_array_missing[idx]
                # positions
                pos_array = vcf_data[pos_key] - 1 # port to BED (0-based) coordinates
                unique_pos, counts_pos = np.unique(pos_array, return_counts=True)
                duplicates = unique_pos[counts_pos > 1]
                if duplicates.any():
                    print("\n[-] Sequence %r: %s VCF records with non-unique positions found. Rescuing records by shifting position... (abort if this is not desired)" % (seq_id, len(duplicates)))
                    pos_array = fix_pos_array(pos_array)
                self.data.create_dataset("seqs/%s/variants/pos" % seq_id, data=pos_array)
        self.data.attrs['seqs/meta/variants/']['count_hom_ref_by_sample'] = count_hom_ref_by_sample
        self.data.attrs['seqs/meta/variants/']['count_hom_alt_by_sample'] = count_hom_alt_by_sample
        self.data.attrs['seqs/meta/variants/']['count_het_by_sample'] = count_het_by_sample
        self.data.attrs['seqs/meta/variants/']['count_missing_by_sample'] = count_missing_by_sample
        self.data.attrs['seqs/meta/variants/']['variant_counts'] = variant_counts
        self.data.attrs['seqs/meta/files/']['vcf_f'] = str(sample_file)

    def info(self, verbose = False):
        print("[=] ==========================================")
        SPACING = 16
        print(f"[+] [{'GStore'.center(SPACING)}] {self.path}")
        sequence_lengths = self.data.attrs['seqs/meta/sequences/']['lengths']
        sequence_span = format_bases(sum(sequence_lengths))
        sequence_count = format_count(len(sequence_lengths))
        sequence_n50 = format_bases(get_n50_from_lengths(sequence_lengths))
        print(f"[+] [{'Genome'.center(SPACING)}] {sequence_span} in {sequence_count} sequence(s) (N50={sequence_n50})")
        samples = self.data.attrs['seqs/meta/samples/']['samples']
        populations = self.data.attrs['seqs/meta/samples/']['populations']
        population_ids = self.data.attrs['seqs/meta/samples/']['population_ids']
        samples_count = len(samples)
        population_count = len(population_ids)
        print(f"[+] [{'Samples'.center(SPACING)}] {samples_count} in {population_count} populations")
        for idx, (population_id, population_char) in enumerate(zip(population_ids, string.ascii_uppercase)):
            population_samples = [sample for sample, population in zip(sample, populations) if population == population_id]
            population_sample_count = len(population_samples)
            print(f"[+] [{'Population {population_char}'.center(SPACING)}] {population_id} has {population_sample_count} samples")
            if verbose:
                population_sample_string = ", ".join(population_samples)
                print(f"[+] [{'Samples in {population_char}'.center(SPACING)}] {population_sample_string}")
        variant_counts = self.data.attrs['seqs/meta/variants/']['variant_counts']
        variant_samples = self.data.attrs['seqs/meta/variants/']['vcf_samples']
        variants_total = sum(variant_counts)
        print(f"[+] [{'Variants'.center(SPACE)}] {variants_total} VCF records for  populations")
        # TABLE
        logging.info("[+] Found %s samples from %s populations. Generated %s sets of %s samples" % (
            len(self.data.attrs['sample_ids']), 
            len(self.data.attrs['population_ids']),
            len(self.data.attrs['sample_sets']),
            pairedness
            ))
        #logging.info("[+] Found %s sequences of a total length of %s b..." % (
        #    len(self.data['seqs/meta/sequences/ids']), 
        #    sum(self.data['seqs/meta/sequences/lengths'])))

    



    def _parse_bed_file(self, parameterObj.bed_file):
        print("[+] Parsing BED file...")
        try:
            bed_df = pd.read_csv(parameterObj.bed_file, sep="\t", usecols=[0,1,2], names=['sequence_id', 'start', 'end'], 
                                    dtype={'sequence_id': str, 'start': 'Int64', 'end': 'Int64'}).sort_values(['sequence_id', 'start'], ascending=[True, True])
        except ValueError:
            sys.exit("[X] BED file %r does not contain the following the columns: 'sequence_id', 'start', 'end'" % (parameterObj.bed_file))
        count_sequences = bed_df['sequence_id'].nunique()
        count_intervals = len(bed_df.index)
        print("[+] Found %s BED intervals on %s sequences (%s intervals/sequence)..." % (
                                            lib.gimble.format_count(count_intervals), 
                                            lib.gimble.format_count(count_sequences), 
                                            lib.gimble.format_fraction(count_intervals / count_sequences)))
        bed_df['distance'] = np.where((bed_df['sequence_id'] == bed_df['sequence_id'].shift(-1)), (bed_df['start'].shift(-1) - bed_df['end']) + 1, np.nan)
        bed_df['length'] = (bed_df['end'] - bed_df['start'])
        distance_counter = collections.Counter(list(bed_df['distance'].dropna(how="any", inplace=False)))
        length_counter = collections.Counter(list(bed_df['length']))
        distance_f = parameterObj.bed_file.parent / (parameterObj.bed_file.stem + '.distance.png')
        plot_loglog(distance_counter, 'Distance to downstream BED interval', distance_f)
        length_f = parameterObj.bed_file.parent / (parameterObj.bed_file.stem + '.length.png')
        plot_loglog(length_counter, 'Length of BED interval', length_f)



class Store(object):

    def _from_zarr(self, parameterObj):
        '''should be made more explicit'''
        self.path = str(parameterObj.zstore)
        self.prefix = str(parameterObj.zstore).rstrip(".z")
        self.data = zarr.open(self.path, mode='r+')

    def _from_input(self, parameterObj):
        self.prefix = parameterObj.outprefix
        self.path = self._get_path(parameterObj.outprefix)
        self.data = zarr.open(self.path, mode='w')

        self._parse_genome_file(parameterObj.genome_file)
        self._parse_sample_file(
            str(parameterObj.sample_file), 
            parameterObj._pairedness
            )
        self._parse_vcf_file(str(parameterObj.vcf_file))
        self.data.attrs['bed_f'] = str(parameterObj.bed_file)
        self.add_stage(parameterObj)
        
    def get_base_count(self, sequence_id=None, kind='total'):
        if kind == 'total':
            if sequence_id is None:
                return sum(self.data.attrs['sequence_length'])
            else:
                if isinstance(sequence_id, list):
                    sequence_ids = set(sequence_id)
                elif isinstance(sequence_id, str):
                    sequence_ids = set([sequence_id])
                else:
                    return None
                return sum([s_length for s_id, s_length in zip(
                    self.data.attrs['sequence_length'], 
                    self.data.attrs['sequence_ids']) if s_id in sequence_ids])

    def get_pop_ids_by_sample_id(self):
        return self.data.attrs['pop_ids_by_sample_id']

    def info(self, verbose=False):
        pop_ids_by_sample_id = self.get_pop_ids_by_sample_id()
        pop_ids = self.data.attrs['pop_ids']
        print("[=] ==========================================")
        print("[+] [DataStore   ] %s" % self.path)
        print("[+] [Genome      ] %s in %s sequence(s)" % (format_bases(self.get_base_count(kind='total')), len(self.data.attrs['sequence_ids'])))
        print("[+] [Samples     ] %s in %s populations" % (len(pop_ids_by_sample_id), len(pop_ids)))
        
        if verbose:
            print("%s" % "\n".join(["[+] \t %s [%s]" % (sample_id, pop_id) for sample_id, pop_id in pop_ids_by_sample_id.items()]))
        print("[+] [Sample-sets ] %s" % (len(self.data.attrs['idx_cartesian_sample_sets'])))
        
        #self.data.attrs['block_length']
        #self.data.attrs['block_span']
        #print("[+] [  Variants  ] : %s variants" % (sum(self.data.attrs['variants'])))
        #interval_counts = []
        #interval_sites = []
        #block_sites = []
        #for sequence_id in self.data.attrs['sequence_ids']:
        #    for sample_set_idx in self.data.attrs['sequence_ids']
        #    interval_counts.append(self.data.attrs['block_span']) 
        #    interval_sites.append(self.data.attrs['block_span'])
        #    block_sites.append()

        #print("[+] [  Intervals ] : %s in %sb" % (self.data.attrs['pairedness'], len(self.data.attrs['idx_cartesian_sample_sets'])))
        #print("[+] [  Blocks    ] l=%s s=%s : %s" % (self.data.attrs['pairedness'], len(self.data.attrs['idx_cartesian_sample_sets'])))

        #self.data.attrs['window_size']
        #self.data.attrs['window_step']
#
        #for sequence_id in self.data.attrs['sequence_ids']:
        #    window = self.data["%s/windows/window_id" % sequence_id]
        #self.data.create_dataset("%s/windows/window_id")
        #self.data.create_dataset("%s/windows/midpoint_mean")
        #self.data.create_dataset("%s/windows/midpoint_median")
        #self.data.create_dataset("%s/%s/interval_sites")
        #self.data.create_dataset("%s/%s/block_sites")
        #self.data.create_dataset("%s/%s/blocks/starts")
        #self.data.create_dataset("%s/%s/blocks/ends")
        #self.data.create_dataset("%s/%s/blocks/variation")
        #self.data.create_dataset("%s/%s/blocks/multiallelic")
        #self.data.create_dataset("%s/%s/blocks/missing")
        #for sequence_id in self.data.attrs['sequence_ids']:

        #print("[+] [Blocks] = %s: %s" % (self.data.attrs['pairedness'], len(self.data.attrs['idx_cartesian_sample_sets'])))
        #print("[+] [Windows] = %s: %s" % (self.data.attrs['pairedness'], len(self.data.attrs['idx_cartesian_sample_sets'])))
        #print("[+] [Grids] = %s: %s" % (self.data.attrs['pairedness'], len(self.data.attrs['idx_cartesian_sample_sets'])))
        print("[=] ==========================================")

    def tree(self):
        return self.data.tree()

    def attrs(self):
        print([self.data.attrs['sample_sets'][idx] for idx in self.data.attrs['idx_cartesian_sample_sets']])
        return "\n".join(
            ["\t".join([k, str(len(v)), str(type(v)), str(v)]) if isinstance(v, list) else "\t".join([k, str(v), str(type(v)), str(v)]) for k, v in self.data.attrs.asdict().items()])


    def _parse_vcf_file(self, vcf_file):
        sample_ids = self.data.attrs['sample_ids']
        sequence_ids = self.data.attrs['sequence_ids']
        variant_counts = []
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            for seq_id in tqdm(sequence_ids, total=len(sequence_ids), desc="[%] Parsing input files ", ncols=100):
                zarr_key = ''
                #print(vcf_file)
                data_by_key = allel.read_vcf(vcf_file, region=seq_id, samples=sample_ids, fields=['samples', 'calldata/GT', 'variants/POS'])
                if data_by_key:
                    for key, data in data_by_key.items():
                        if key == 'samples':
                            self.data.attrs['sample_ids_vcf'] = list(data)
                            #print("vcf samples", list(data))
                            self.data.attrs['sample_ids_to_vcf_idx'] = {sample_id: idx for idx, sample_id in enumerate(data)}
                        elif key == 'calldata/GT':
                            zarr_key = "%s/gt" % seq_id
                            self.data.create_dataset(zarr_key, data=data) 
                        elif key == 'variants/POS':
                            zarr_key = "%s/pos" % seq_id
                            data = np.array(data)
                            variant_counts.append(data.shape[0])
                            #print('np.array(data).shape', np.array(data).shape)
                            unique_pos, counts_pos = np.unique(data, return_counts=True)
                            duplicates = unique_pos[counts_pos > 1]
                            if duplicates.any():
                                print("\n[-] Sequence %r: %s VCF records with non-unique positions found. Rescuing records by shifting position... (abort if this is not desired)" % (seq_id, len(duplicates)))
                                data = fix_pos_array(data)
                            self.data.create_dataset(zarr_key, data=data - 1) # port to BED (0-based) coordinates
                        else:
                            logging.error("[X] Unknown key %r" % key)
                            sys.exit()
            self.data.attrs['variants'] = variant_counts
            self.data.attrs['vcf_f'] = str(vcf_file)

    def has_data(self, datatype):
        if datatype in self.data.attrs:
            return True
        return False

    def make_blocks(self, parameterObj, debug=False):
        if self.has_data('blocks'):
            if not parameterObj.overwrite:
                sys.exit("[X] Store %r already contains blocks.\n[X] These blocks => %r\n[X] Please specify '--force' to overwrite." % (self.path, self.get_stage_cmd('blocks')))
            print('[-] Store %r already contains blocks. But these will be overwritten...' % (self.path))
        self.data.attrs['block_length'] = parameterObj.block_length
        self.data.attrs['block_gap_run'] = parameterObj.block_gap_run
        self.data.attrs['block_span'] = parameterObj.block_span
        self.data.attrs['block_max_missing'] = parameterObj.block_max_missing
        self.data.attrs['block_max_multiallelic'] = parameterObj.block_max_multiallelic
        self.data.attrs['mutypes_count'] = 4 # should be calculated from possible folded genotypes and pairedness
        sample_ids = self.data.attrs['sample_ids'] # from BED file
        sequence_ids = self.data.attrs['sequence_ids']
        sample_sets = self.data.attrs['sample_sets']
        #print(self.attrs())
        logging.info("[#] Processing BED file %r ..." % self.data.attrs['bed_f'])
        df = pd.read_csv(self.data.attrs['bed_f'], sep="\t", usecols=[0, 1, 2, 4], names=['sequence_id', 'start', 'end', 'samples'], dtype={'sequence_id': str, 'start': np.int, 'end': np.int, 'samples': str})
        # remove sequence_ids that are not in sequence_names_array, sort, reset index
        intervals_df = df[df['sequence_id'].isin(sequence_ids)].sort_values(['sequence_id', 'start'], ascending=[True, True]).reset_index(drop=True)
        # get length column
        intervals_df['length'] = intervals_df['end'] - intervals_df['start'] 
        # get coverage matrix and drop columns of samples that are not in sample_ids_array
        intervals_df = pd.concat([intervals_df, intervals_df.samples.str.get_dummies(sep=',').filter(sample_ids)], axis=1).drop(columns=['samples'])
        # remove those intervals including less than two sample_ids (query sample_id columns : intervals_df[intervals_df.columns.intersection(sample_ids)])
        intervals_df = intervals_df.loc[(intervals_df[intervals_df.columns.intersection(sample_ids)].sum(axis=1) > 1)]
        #interval_sites_by_sample_set_idx = collections.defaultdict(list)
        #block_sites_by_sample_set_idx = collections.defaultdict(list)
        with tqdm(total=(len(sequence_ids) * len(sample_sets)), desc="[%] Calculating bSFSs ", ncols=100, unit_scale=True) as pbar:        
            for seq_id in sequence_ids:        
                _intervals_df = intervals_df[intervals_df['sequence_id'] == seq_id]
                # To Do: put a check in so that it only progresses if there ARE intervals in the BED file for that sequence_id ... 
                seq_id_key = "%s/pos" % seq_id
                _pos = np.array([])
                if seq_id_key in self.data: 
                    _pos = self.data["%s/pos" % seq_id] # zarr.core.array
                    sa_genotype_array = allel.GenotypeArray(self.data["%s/gt" % seq_id])
                for sample_set_idx, sample_set in enumerate(sample_sets):
                    if sample_set_idx in self.data.attrs['idx_cartesian_sample_sets']:
                        try:
                            sample_set_intervals_df = _intervals_df[_intervals_df[sample_set].all(axis='columns')]
                        except KeyError:
                            sys.exit("[X] Sample set %r not found in BED file" % sample_set)
                        # Cut blocks based on intervals and block-algoritm parameters
                        block_sites = cut_blocks(
                            np.array(sample_set_intervals_df.start), 
                            np.array(sample_set_intervals_df.end), 
                            parameterObj.block_length, 
                            parameterObj.block_span, 
                            parameterObj.block_gap_run
                            )
                        block_starts = np.array(block_sites[:,0])           # has to be np.array() !
                        block_ends = np.array(block_sites[:,-1] + 1)        # has to be np.array() !
                        interval_space = sample_set_intervals_df['length'].sum()
                        block_space = np.sum(((block_sites[:,-1] - block_sites[:,0]) + 1)) # block_space == span !
                        if debug:
                            print("#", seq_id, sample_set_idx, sample_set)
                            print("# Block_sites 1", block_sites.shape)
                            print(block_sites)
                        # Variation
                        if seq_id_key in self.data: 
                            # positions inside block_sites that are variant     
                            idx_block_sites_in_pos = np.isin(block_sites, _pos, assume_unique=True) # will crash if non-unique pos
                            # variant positions in _pos that are within block_sites 
                            idx_pos_in_block_sites = np.isin(_pos, block_sites, assume_unique=True) # will crash if non-unique pos
                            # make room for reading in array
                            sample_set_vcf_idxs = [self.data.attrs['sample_ids_to_vcf_idx'][sample_id] for sample_id in sample_set] # list of indices
                            sa_sample_set_genotype_array = sa_genotype_array.subset(idx_pos_in_block_sites, sample_set_vcf_idxs)
                            block_sites = genotype_to_mutype_array(sa_sample_set_genotype_array, idx_block_sites_in_pos, block_sites, debug)
                        else:
                            block_sites[:] = 2 # if no variants, all invariant
                        multiallelic, missing, monomorphic, variation = block_sites_to_variation_arrays(block_sites, self.data.attrs['mutypes_count'])
                        pi_1, pi_2, d_xy, f_st, fgv = calculate_popgen_from_array(variation, (self.data.attrs['block_length'] * variation.shape[0])) 
                        if debug:
                            print("# Block_sites 2")
                            print(block_sites)
                            print('# Variation: 0=HetB, 1=HetA, 2=HetAB, 3=Fixed')
                            print(variation)
                            print("[+] Pi_%s = %s; Pi_%s = %s; D_xy = %s; F_st = %s; FGV = %s" % (self.data.attrs['pop_ids'][0], pi_1, self.data.attrs['pop_ids'][1], pi_2, d_xy, f_st, fgv)) 
                        # To Do: put explicit datatypes when saving stuff in ZARR (it might otherwise estimate it wrongly!)
                        self.data.create_dataset("%s/%s/interval_sites" % (seq_id, sample_set_idx), data=np.array([interval_space]), overwrite=True)
                        self.data.create_dataset("%s/%s/block_sites" % (seq_id, sample_set_idx), data=np.array([block_space]), overwrite=True)
                        self.data.create_dataset("%s/%s/blocks/starts" % (seq_id, sample_set_idx), data=block_starts, overwrite=True)
                        self.data.create_dataset("%s/%s/blocks/ends" % (seq_id, sample_set_idx), data=block_ends, overwrite=True)
                        self.data.create_dataset("%s/%s/blocks/variation" % (seq_id, sample_set_idx), data=variation, overwrite=True)
                        self.data.create_dataset("%s/%s/blocks/multiallelic" % (seq_id, sample_set_idx), data=multiallelic.flatten(), overwrite=True)
                        self.data.create_dataset("%s/%s/blocks/missing" % (seq_id, sample_set_idx), data=missing.flatten(), overwrite=True)
                    pbar.update(1)
        self.add_stage(parameterObj)

    def make_windows(self, parameterObj):
        block_status = self.has_data('blocks')
        if not block_status:
            sys.exit("[X] No blocks found. Please make blocks first.")
        if self.has_data('windows'):
            if not parameterObj.overwrite:
                sys.exit("[X] Store %r already contains windows.\n[X] These windows => %r\n[X] Please specify '--force' to overwrite." % (self.path, self.get_stage_cmd('windows')))
            print('[-] Store %r already contains windows. But these will be overwritten...' % (self.path))
        sample_sets_idxs = self.data.attrs['idx_cartesian_sample_sets']
        with tqdm(total=(len(self.data.attrs['sequence_ids']) * len(sample_sets_idxs)), desc="[%] Making windows ", ncols=100, unit_scale=True) as pbar: 
            for seq_id in self.data.attrs['sequence_ids']: 
                mutypes, sample_set_covs, starts, ends = [], [], [], []
                for sample_set_idx in sample_sets_idxs:
                    # first determine valid block mask
                    missing = np.array(self.data["%s/%s/blocks/missing" % (seq_id, sample_set_idx)])
                    multiallelic = np.array(self.data["%s/%s/blocks/multiallelic" % (seq_id, sample_set_idx)])
                    valid = np.less_equal(missing, self.data.attrs['block_max_missing']) & np.less_equal(multiallelic, self.data.attrs['block_max_multiallelic'])
                    mutype = np.array(self.data["%s/%s/blocks/variation" % (seq_id, sample_set_idx)])[valid]
                    mutypes.append(mutype)
                    start = np.array(self.data["%s/%s/blocks/starts" % (seq_id, sample_set_idx)])[valid]
                    starts.append(start)
                    end = np.array(self.data["%s/%s/blocks/ends" % (seq_id, sample_set_idx)])[valid]
                    ends.append(end)
                    sample_set_covs.append(np.full_like(end, sample_set_idx)) 
                    pbar.update()
                mutype_array = np.concatenate(mutypes, axis=0)
                #print('mutype_array', mutype_array.shape, mutype_array[:])
                start_array = np.concatenate(starts[:], axis=0)
                #print('starts', starts)
                end_array = np.concatenate(ends, axis=0)
                #print('end_array', end_array.shape, end_array[:])
                sample_set_array = np.concatenate(sample_set_covs, axis=0)
                window_variation, window_starts, window_ends, window_pos_mean, window_pos_median = cut_windows(mutype_array, sample_sets_idxs, start_array, end_array, sample_set_array, num_blocks=parameterObj.window_size, num_steps=parameterObj.window_step)
                self.data.create_dataset("%s/windows/variation" % seq_id, data=window_variation, overwrite=True)
                self.data.create_dataset("%s/windows/starts" % seq_id, data=window_starts, overwrite=True)
                self.data.create_dataset("%s/windows/ends" % seq_id, data=window_ends, overwrite=True)
                self.data.create_dataset("%s/windows/pos_mean" % seq_id, data=window_pos_mean, overwrite=True)
                self.data.create_dataset("%s/windows/pos_median" % seq_id, data=window_pos_median, overwrite=True)
        self.add_stage(parameterObj)
        self.data.attrs['window_size'] = parameterObj.window_size
        self.data.attrs['window_step'] = parameterObj.window_step

    def _get_path(self, outprefix):
        path = "%s.z" % outprefix
        if os.path.isdir(path):
            logging.info("[!] ZARR store %r exists. Deleting ..." % path)
            shutil.rmtree(path)
        logging.info("[+] Generating ZARR store %r" % path)
        return path

    def get_gts(self, sequence_ids, start, end, sample_ids):
        pass

    def dump_blocks(self, parameterObj):
        # applies the missing & multiallelic thresholds
        sample_sets_idxs = self.data.attrs['idx_cartesian_sample_sets']
        data_by_key_by_sample_set_idx = collections.defaultdict(lambda: collections.defaultdict(list))
        variation_global = []
        with tqdm(total=(len(self.data.attrs['sequence_ids']) * len(sample_sets_idxs)), desc="[%] Writing bSFSs ", ncols=100, unit_scale=True) as pbar: 
            for seq_id in self.data.attrs['sequence_ids']: 
                for sample_set_idx in sample_sets_idxs:
                    # aggregate interval/block sites as determined from BED file and blocking algorithm
                    data_by_key_by_sample_set_idx[sample_set_idx]['interval_sites'].append(np.array(self.data["%s/%s/interval_sites" % (seq_id, sample_set_idx)]))
                    data_by_key_by_sample_set_idx[sample_set_idx]['block_sites'].append(np.array(self.data["%s/%s/block_sites" % (seq_id, sample_set_idx)]))
                    # missing & multiallelic thresholds determine validity of blocks ...
                    missing = np.array(self.data["%s/%s/blocks/missing" % (seq_id, sample_set_idx)])
                    multiallelic = np.array(self.data["%s/%s/blocks/multiallelic" % (seq_id, sample_set_idx)])
                    valid = np.less_equal(missing, parameterObj.block_max_missing) & np.less_equal(multiallelic, parameterObj.block_max_multiallelic)
                    data_by_key_by_sample_set_idx[sample_set_idx]['block_sites_valid'].append(np.array([valid[valid == True].shape[0] * self.data.attrs['block_length']]))
                    # aggregate variation/location data of valid blocks 
                    data_by_key_by_sample_set_idx[sample_set_idx]['missing'].append(missing[valid])
                    data_by_key_by_sample_set_idx[sample_set_idx]['multiallelic'].append(multiallelic[valid])
                    variation = np.array(self.data["%s/%s/blocks/variation" % (seq_id, sample_set_idx)])[valid]
                    data_by_key_by_sample_set_idx[sample_set_idx]['variation'].append(variation)
                    variation_global.append(variation)
                    data_by_key_by_sample_set_idx[sample_set_idx]['starts'].append(np.array(self.data["%s/%s/blocks/starts" % (seq_id, sample_set_idx)])[valid])
                    data_by_key_by_sample_set_idx[sample_set_idx]['ends'].append(np.array(self.data["%s/%s/blocks/ends" % (seq_id, sample_set_idx)])[valid])
                    #sample_set.append(np.full_like(end[valid], sample_set_idx)) 
                    pbar.update()
        variation_global_array = np.concatenate(variation_global, axis=0)
        # popgen
        variation_global = []
        metrics_rows = []
        # is order (pi_1, pi_2, d_xy, f_st, fgv) correct?
        for sample_set_idx in data_by_key_by_sample_set_idx:
            sample_set_ids = self.data.attrs['sample_sets'][sample_set_idx]
            #print(data_by_key_by_sample_set_idx)
            block_sites = np.sum(np.concatenate(data_by_key_by_sample_set_idx[sample_set_idx]['block_sites'], axis=0))
            interval_sites = np.sum(np.concatenate(data_by_key_by_sample_set_idx[sample_set_idx]['interval_sites'], axis=0))
            block_sites_valid = np.sum(np.concatenate(data_by_key_by_sample_set_idx[sample_set_idx]['block_sites_valid'], axis=0))
            variation_array = np.concatenate(data_by_key_by_sample_set_idx[sample_set_idx]['variation'], axis=0)
            missing_count = np.sum(np.concatenate(data_by_key_by_sample_set_idx[sample_set_idx]['missing'], axis=0))
            multiallelic_count = np.sum(np.concatenate(data_by_key_by_sample_set_idx[sample_set_idx]['multiallelic'], axis=0))
            hetB_count, hetA_count, hetAB_count, fixed_count = np.sum(variation_array, axis=0)
            #pi_1, pi_2, d_xy, f_st, fgv = calculate_popgen_from_array(variation_array, (self.data.attrs['block_length'] * variation_array.shape[0]))    
            pi_1, pi_2, d_xy, f_st, fgv = calculate_popgen_from_array(variation_array, block_sites_valid)    
            metrics_rows.append([
                sample_set_ids[0], 
                sample_set_ids[1],     
                block_sites,
                interval_sites,
                block_sites_valid,
                np.divide(block_sites_valid, self.data.attrs['block_length']),
                fgv,
                missing_count,
                multiallelic_count,
                hetA_count, 
                hetB_count, 
                hetAB_count, 
                fixed_count,
                pi_1, 
                pi_2, 
                d_xy, 
                f_st
                ])
        # output metrics 
        header = [
            self.data.attrs['pop_ids'][0], 
            self.data.attrs['pop_ids'][1], 
            'block_sites', 
            'interval_sites', 
            'block_sites_valid', 
            'blocks', 
            'fgv', 
            'missing', 
            'multiallelic', 
            'hetA', 
            'hetB', 
            'hetAB', 
            'fixed', 
            'piA', 
            'piB', 
            'dxy', 
            'fst'
            ]
        pd.DataFrame(data=metrics_rows, columns=header, dtype='int64').to_hdf("%s.block_stats.h5" % self.prefix, 'bsfs', format='table')

        pi_1, pi_2, d_xy, f_st, fgv = calculate_popgen_from_array(variation_global_array, (self.data.attrs['block_length'] * variation_global_array.shape[0]))
        print("[+] Pi_%s = %s; Pi_%s = %s; D_xy = %s; F_st = %s; FGVs = %s / %s blocks (%s)" % (self.data.attrs['pop_ids'][0], pi_1, self.data.attrs['pop_ids'][1], pi_2, d_xy, f_st, fgv, variation_global_array.shape[0], format_percentage(fgv / variation_global_array.shape[0]))) 
        
        # mutuple barchart
        mutypes, counts = np.unique(variation_global_array, return_counts=True, axis=0)
        mutype_counter = collections.Counter({tuple(i):j for i,j in zip(mutypes, counts)})
        plot_mutuple_barchart('%s.mutuple_barchart.png' % self.prefix, mutype_counter)

        # mutuple tally
        bsfs = np.concatenate([counts[:, np.newaxis], mutypes], axis =-1)
        header = ['count'] + [x+1 for x in range(self.data.attrs['mutypes_count'])]
        pd.DataFrame(data=bsfs, columns=header, dtype='int64').to_hdf("%s.blocks.h5" % self.prefix, 'tally', format='table')

        # block coordinates (BED format)
        #header = ['block_id', 'start', 'end', 'sample_set', 'multiallelic', 'missing']
        #pd.DataFrame(data=bsfs, columns=header, dtype='int64').to_hdf("%s.blocks.h5" % self.prefix, 'bed', format='table')

    def dump_windows(self, parameterObj):
        window_info_rows = []
        window_mutuple_tally = []
        for sequence_id in tqdm(self.data.attrs['sequence_ids'], total=len(self.data.attrs['sequence_ids']), desc="[%] Generating output ", ncols=100):
            variations = self.data["%s/windows/variation" % sequence_id]
            #print(self.data["%s/windows/starts" % sequence_id][:])
            #print(self.data["%s/windows/pos_mean" % sequence_id][:])
            #print(self.data["%s/windows/pos_median" % sequence_id][:])
            window_ids = np.array(["_".join([sequence_id, _start, _end]) for (_start, _end) in zip(
                np.array(self.data["%s/windows/starts" % sequence_id]).astype(str), 
                np.array(self.data["%s/windows/ends" % sequence_id]).astype(str))])
            #window_ids = self.data["%s/windows/window_id" % sequence_id]
            midpoint_means = self.data["%s/windows/pos_mean" % sequence_id]
            midpoint_medians = self.data["%s/windows/pos_median" % sequence_id]
            for window_id, variation, midpoint_mean, midpoint_median in zip(window_ids, variations, midpoint_means, midpoint_medians):
                pi_1, pi_2, d_xy, f_st, fgv = calculate_popgen_from_array(variation, (self.data.attrs['block_length'] * variation.shape[0]))
                window_info_rows.append([window_id, sequence_id, midpoint_mean, midpoint_median, pi_1, pi_2, d_xy, f_st, fgv/variation.shape[0]])
                # mutuple barchart
                mutypes, counts = np.unique(variation, return_counts=True, axis=0)
                tally = np.concatenate([counts[:, np.newaxis], mutypes], axis =-1)
                windows = np.array([window_id] * tally.shape[0])
                window_mutuple_tally.append(np.concatenate([windows[:, np.newaxis], tally], axis =-1))
        
        window_bsfs_cols = ['window_id', 'count'] + [x+1 for x in range(self.data.attrs['mutypes_count'])]
        window_bsfs_df = pd.DataFrame(np.vstack(window_mutuple_tally), columns=window_bsfs_cols)
        print("[+] Made %s windows" % window_bsfs_df['window_id'].nunique()) 
        window_bsfs_f = "%s.window_bsfs.tsv" % self.prefix
        window_bsfs_df.to_csv(window_bsfs_f, sep='\t', index=False)
        print("[>] Created: %r" % str(window_bsfs_f))

        window_info_cols = ['window_id', 'sequence_id', 'midpoint_mean', 'midpoint_median', 'pi_%s' % self.data.attrs['pop_ids'][0], 'pi_%s' % self.data.attrs['pop_ids'][1], 'd_xy', 'f_st', 'fgv']
        window_info_df = pd.DataFrame(window_info_rows, columns=window_info_cols)
        window_info_f = "%s.window_info.tsv" % self.prefix
        window_info_df.to_csv(window_info_f, sep='\t', index=False)
        print("[>] Created: %r" % str(window_info_f))        
        self.plot_fst_genome_scan(window_info_df)
        self.plot_pi_genome_scan(window_info_df)
        #     plot_pi_scatter(window_df, '%s.pi_scatter.png' % parameterObj.dataset)
    
    def plot_pi_genome_scan(self, window_df):
        offset_by_sequence_id = {}
        offset = 0
        x_boundaries = []
        for sequence_id, sequence_length in zip(self.data.attrs['sequence_ids'], self.data.attrs['sequence_length']):
            offset_by_sequence_id[sequence_id] = offset
            x_boundaries.append(offset)
            offset += sequence_length
        x_boundaries.append(offset)
        #print([(sequenceObj.id, sequenceObj.length) for sequenceObj in sequenceObjs])
        #print(x_boundaries)
        fig = plt.figure(figsize=(18,4), dpi=200, frameon=True)
        #connecting dots
        ax = fig.add_subplot(111)  
        window_df['rel_pos'] = window_df['midpoint_median'] + window_df['sequence_id'].map(offset_by_sequence_id)
        window_df.sort_values(['rel_pos'], inplace=True)
        #print(window_df)
        pi_A_key = list(window_df.columns)[4]
        pi_B_key = list(window_df.columns)[5]
        ax.plot(window_df['rel_pos'], window_df[pi_A_key], color='orange', alpha=0.8, linestyle='-', linewidth=1, label=pi_A_key.replace('pi_', ''))
        ax.plot(window_df['rel_pos'], window_df[pi_B_key], color='dodgerblue', alpha=0.8, linestyle='-', linewidth=1, label=pi_B_key.replace('pi_', ''))
        y_lim = (min(window_df[pi_A_key].min(), window_df[pi_B_key].min()), max(window_df[pi_A_key].max(), window_df[pi_B_key].max()))
        ax.vlines(x_boundaries, y_lim[0], y_lim[1], colors=['lightgrey'], linestyles='dashed', linewidth=1)
        ax.set_ylim(y_lim)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.legend(numpoints=1)
        plt.ylabel('Pi')
        plt.xlabel("Genome coordinate")
        out_f = '%s.pi_genome_scan.png' % self.prefix
        plt.tight_layout()
        fig.savefig(out_f, format="png")
        print("[>] Created: %r" % str(out_f))
        plt.close(fig)

    def plot_fst_genome_scan(self, window_df):
        offset_by_sequence_id = {}
        offset = 0
        x_boundaries = []
        for sequence_id, sequence_length in zip(self.data.attrs['sequence_ids'], self.data.attrs['sequence_length']):
            offset_by_sequence_id[sequence_id] = offset
            x_boundaries.append(offset)
            offset += sequence_length
        x_boundaries.append(offset)
        fig = plt.figure(figsize=(18,4), dpi=200, frameon=True)
        #connecting dots
        ax = fig.add_subplot(111)  
        y_lim = (0.0, 1.0)
        window_df['rel_pos'] = window_df['midpoint_median'] + window_df['sequence_id'].map(offset_by_sequence_id)
        window_df.sort_values(['rel_pos'], inplace=True)
        ax.plot(window_df['rel_pos'], window_df['f_st'], color='lightgrey', alpha=0.8, linestyle='-', linewidth=1)
        scatter = ax.scatter(window_df['rel_pos'], window_df['f_st'], c=window_df['d_xy'], alpha=1.0, cmap='PiYG_r', edgecolors='white', marker='o', s=40, linewidth=0.2)
        cbar = fig.colorbar(scatter, ax=ax)
        cbar.ax.set_title('D_xy')
        ax.vlines(x_boundaries, 0.0, 1.0, colors=['lightgrey'], linestyles='dashed', linewidth=1)
        ax.set_ylim(y_lim)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        plt.ylabel('F_st')
        plt.xlabel("Genome coordinate")
        ax.autoscale_view(tight=None, scalex=True, scaley=True)
        out_f = '%s.fst_genome_scan.png' % self.prefix
        fig.savefig(out_f, format="png")
        plt.close(fig)
        print("[>] Created: %r" % str(out_f))