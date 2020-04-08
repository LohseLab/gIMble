#import collections
import itertools
from tqdm import tqdm
#from lib.functions import get_pop_ids_by_sample_ids_from_csv, plot_mutuple_barchart, format_bases, format_fraction, create_hdf5_store, tabulate_df, poolcontext, format_percentage, plot_distance_scatter, plot_fst_genome_scan, plot_pi_genome_scan, plot_pi_scatter, plot_sample_barchart
from lib.functions import plot_mutuple_barchart, format_percentage
from sys import exit
#import pathlib
import allel
import numpy as np
import pandas as pd
import shutil
import zarr
import os
import logging
import collections
import sys
import warnings

#from timeit import default_timer as timer
'''
[To do]
- write test that errors when no VCF index since it takes ages without

[ QC ]
- Genotype counts for each sample

[ setup ] 
- is it possible to run it without population file?
    - concept of samples has to come from VCF/BED

[ analyse ]
- parse population file and use it for query set

~/git/gIMble/gIMble setup -s ~/data/heliconius.populations.csv -v ~/data/heliconius.freebayes.autosomes.intergenic_and_nonrepeats.normalized.SnpGap_2.NonSNP.Balance.PASS.decomposed.vcf.gz -b ~/data/heliconius.master.normalised.slop_2.bed -o ~/data/heliconius_full_all -g ~/data/hmel2_5.chromosomes.autosomes.genomefile 
&& ~/git/gIMble/gIMble blocks -z heliconius_full_all.z -l 64 -r 5 -m 80 
&& ~/git/gIMble/gIMble windows -z heliconius_full_all.z -w 50000 -s 10000

        https://stupidpythonideas.blogspot.com/2014/01/grouping-into-runs-of-adjacent-values.html
        https://stackoverflow.com/questions/2154249/identify-groups-of-continuous-numbers-in-a-list
        https://stackoverflow.com/questions/4494404/find-large-number-of-consecutive-values-fulfilling-condition-in-a-numpy-array
        https://stackoverflow.com/questions/32357823/splitting-a-list-of-integers-by-breaks-in-the-data-when-the-break-size-is-variab

# # to deal with samples n!=2, it 
# # - needs to know the number of folded mutypes
# # - or needs to calculate the number of fMiC (folded minor allele counts)

# def genotypes(alleles=2, ploidy=2, samples=2):
#     return [list(x) for x in itertools.product([u for u in itertools.combinations_with_replacement(range(alleles), ploidy)], repeat=samples)]

# def genotype_count(alleles=2, ploidy=2, samples=2):
#     return len(genotypes(alleles, ploidy, samples))


'''
COLOURS = ['orange', 'dodgerblue']
FULL_MUTYPE_ORDER = ['hetA', 'fixed', 'hetB', 'hetAB', 'missing', 'multiallelic']
MUTYPE_ORDER = ['hetA', 'fixed', 'hetB', 'hetAB']
MUTYPE_OTHER = ['missing', 'multiallelic']
GT_ORDER = ['TOTAL', 'MISS', 'HOM', 'HET']

np.set_printoptions(threshold=sys.maxsize)

import matplotlib.pyplot as plt

def plot_histogram(x, out_f):
    fig, ax = plt.subplots(figsize=(14, 5))
    hist, bins = np.histogram(x, bins=50)
    width = 0.7 * (bins[1] - bins[0])
    center = (bins[:-1] + bins[1:]) / 2
    ax.bar(center, hist, align='center', width=width)
    fig.savefig('%s.png' % out_f, format="png")

######################## inference ################################

'''
Model file should yield:
    - n Paths with different amounts of Steps
        
'''


##############################################################
def szudzik_pairing(folded_minor_allele_counts):
    # folded_minor_allele_counts.shape -> (n, 2)
    # adapted from: https://drhagen.com/blog/superior-pairing-function/
    # for unpairing and multidimensional pairing functions, see: https://drhagen.com/blog/multidimensional-pairing-functions/
    return np.where(
        (folded_minor_allele_counts[:,0] >= folded_minor_allele_counts[:,1]),
        np.square(folded_minor_allele_counts[:,0]) + folded_minor_allele_counts[:,0] + folded_minor_allele_counts[:,1],
        folded_minor_allele_counts[:,0] + np.square(folded_minor_allele_counts[:,1])
        )

def get_coverage_counts(coverages, idxs, num_blocks):
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
    # print(df)
    #print('mutype_array', mutype_array.shape, mutype_array)
    #print('sites', sites)
    pi_1 = float("%.8f" % np.divide(np.sum(mutype_array[:, 0]) + np.sum(mutype_array[:, 2]), sites))
    pi_2 = float("%.8f" % np.divide(np.sum(mutype_array[:, 1]) + np.sum(mutype_array[:, 2]), sites))
    d_xy = float("%.8f" % np.divide(np.divide(np.sum(mutype_array[:, 0]) + np.sum(mutype_array[:, 1]) + np.sum(mutype_array[:, 2]), 2.0) + np.sum(mutype_array[:, 3]), sites))
    mean_pi = (pi_1 + pi_2) / 2.0
    total_pi = (d_xy + mean_pi) / 2.0 # special case of pairwise Fst
    f_st = np.nan
    if (total_pi):
        f_st = float("%.8f" % ((total_pi - mean_pi) / total_pi)) # special case of pairwise Fst
    fgv = len(mutype_array[(mutype_array[:, 2] > 0) & (mutype_array[:, 3] > 0)])
    return (pi_1, pi_2, d_xy, f_st, fgv)

def genotype_to_mutype_array(genotypeArray):
    gt_array = np.array(genotypeArray)
    allele_counts = np.ma.masked_equal(genotypeArray.count_alleles(), 0, copy=False)    
    allele_map = np.ones((allele_counts.shape), dtype='int8') * np.arange(allele_counts.shape[-1])
    idx_max_global_allele_count = np.nanargmax(allele_counts, axis=1)
    idx_min_global_allele_count = np.nanargmin(allele_counts, axis=1)
    has_major_allele = (idx_max_global_allele_count != idx_min_global_allele_count)
    idx_min_prime_allele = np.amin(gt_array[:,0], axis=1)
    idx_min_global_allele = np.amin(np.amin(gt_array, axis=1), axis=1)
    idx_max_global_allele = np.amax(np.amax(gt_array, axis=1), axis=1)
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
    allele_map[np.arange(allele_map.shape[0]), idx_minor_allele] = 1 # First "minor" allele, so that it gets overwritten if monomorphic
    allele_map[np.arange(allele_map.shape[0]), idx_major_allele] = 0
    folded_genotypes = genotypeArray.map_alleles(allele_map) #allel.GenotypeArray(genotypeArray.map_alleles(allele_map))
    folded_minor_allele_counts = folded_genotypes.to_n_alt(fill=-1)
    folded_minor_allele_counts[np.any(genotypeArray.is_missing(), axis=1)] = np.ones(2) * -1    # -1, -1 for missing => -1
    folded_minor_allele_counts[(allele_counts.count(axis=1) > 2)] = np.ones(2) * (-1, -2)       # -1, -2 for multiallelic => -2
    mutypes = szudzik_pairing(folded_minor_allele_counts) + 2                                   # add 2 so that not negative for bincount
    return mutypes
    
def consecutive(data, stepsize=1):
    return np.split(data, np.where(np.diff(data) != stepsize)[0] + 1)

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
    window_starts = start_array.take(coordinate_sorted_idx, axis=0).take(window_idxs, axis=0)
    window_min = np.min(start_array.take(coordinate_sorted_idx, axis=0).take(window_idxs, axis=0), axis=1).T
    window_ends = end_array.take(coordinate_sorted_idx, axis=0).take(window_idxs, axis=0)#[:,-1]
    window_max = np.max(end_array.take(coordinate_sorted_idx, axis=0).take(window_idxs, axis=0), axis=1).T
    window_midpoints = (window_starts / 2) + (window_ends / 2)
    window_mean = np.mean(window_midpoints, axis=1).T
    window_median = np.median(window_midpoints, axis=1).T
    # coverages = get_coverage_counts(sample_set_array.take(window_idxs, axis=0), idxs, num_blocks)
    return window_mutypes, window_min, window_max, window_mean, window_median

def cut_blocks(interval_starts, interval_ends, block_length, block_span, block_gap_run):
    sites = create_ranges(np.array((interval_starts, interval_ends)).T)
    block_sites = np.concatenate([x[:block_length * (x.shape[0] // block_length)].reshape(-1, block_length) for x in np.split(sites, np.where(np.diff(sites) > block_gap_run)[0] + 1)])
    return block_sites[(((block_sites[:, -1] - block_sites[:, 0]) + 1) <= block_span)] # return blocks where ((block_end - block_start) + 1) <= block_span

def create_store(parameterObj):
    store = Store() 
    store._from_input(parameterObj)
    return store

def load_store(parameterObj):
    store = Store() 
    store._from_zarr(parameterObj)
    return store

class Store(object):
    def __init__(self):
        self.path = None
        self.prefix = None
        self.data = None
        self.stages = {}

    def _from_zarr(self, parameterObj):
        self.path = parameterObj.zstore
        self.prefix = parameterObj.zstore.rstrip(".z")
        self.data = zarr.open(self.path, mode='r+')

    def _from_input(self, parameterObj):
        self.prefix = parameterObj.outprefix
        self.path = self._get_path(parameterObj.outprefix)
        self.data = zarr.open(self.path, mode='w')
        self._parse_genome_file(parameterObj.genome_file)
        self._parse_sample_file(
            parameterObj.sample_file, 
            parameterObj.pairedness
            )
        self._parse_vcf_file(parameterObj.vcf_file)
        self.data.attrs['bed_f'] = parameterObj.bed_file 
        
    def tree(self):
        return self.data.tree()

    def attrs(self):
        print([self.data.attrs['sample_sets'][idx] for idx in self.data.attrs['idx_cartesian_sample_sets']])
        return "\n".join(
            ["\t".join([k, str(len(v)), str(type(v)), str(v)]) if isinstance(v, list) else "\t".join([k, str(v), str(type(v)), str(v)]) for k, v in self.data.attrs.asdict().items()])

    def _parse_genome_file(self, genome_file):
        logging.info("[#] Processing Genome file %r..." % genome_file)
        df = pd.read_csv(genome_file, sep="\t", names=['sequence_id', 'sequence_length'], header=None)
        self.data.attrs['sequence_ids'] = df['sequence_id'].to_list() #df['sequence_id'].to_numpy(dtype=str)
        self.data.attrs['sequence_length'] = df['sequence_length'].to_list()
        logging.info("[+] Found %s sequences of a total length of %s b..." % (len(self.data.attrs['sequence_ids']), sum(self.data.attrs['sequence_length'])))
        for sequence_id in self.data.attrs['sequence_ids']:
            self.data.create_group('%s/' % sequence_id)
        self.data.attrs['genome_f'] = genome_file
        
    def _parse_sample_file(self, sample_file, pairedness):
        logging.info("[#] Processing Sample file %r ..." % sample_file)
        df = pd.read_csv(sample_file, sep=",", names=['sample_id', 'population_id'])
        self.data.attrs['sample_ids'] = df['sample_id'].to_list() # Order as they appear in file
        self.data.attrs['population_ids'] = df['population_id'].to_list() 
        self.data.attrs['pop_ids'] = sorted(set(df['population_id'].to_list())) # Sorted pop_ids
        self.data.attrs['pop_ids_by_sample_id'] = {sample_id: pop_id for sample_id, pop_id in zip(self.data.attrs['sample_ids'], self.data.attrs['population_ids'])}
        self.data.attrs['sample_sets'] = [_ for _ in itertools.combinations(self.data.attrs['sample_ids'], pairedness)] # DO NOT SORT, otherwise there could be erroneously polarised sample sets (pop2, pop1)
        self.data.attrs['idx_cartesian_sample_sets'] = [idx for idx, sample_set in enumerate(self.data.attrs['sample_sets']) if (len(set([self.data.attrs['pop_ids_by_sample_id'][sample_id] for sample_id in sample_set])) == len(self.data.attrs['pop_ids']))]
        # TEST FOR CORRECT SAMPLE_SET POLARISATION
        # _pops = []
        # for sample_set in self.data.attrs['sample_sets']:
        #     _pops.append(tuple(
        #         [self.data.attrs['pop_ids_by_sample_id'][sample_set[0]],
        #         self.data.attrs['pop_ids_by_sample_id'][sample_set[1]]]
        #         ))
        # print(collections.Counter(_pops))
        logging.info("[+] Found %s samples from %s populations. Generated %s sets of %s samples" % (
            len(self.data.attrs['sample_ids']), 
            len(self.data.attrs['population_ids']),
            len(self.data.attrs['sample_sets']),
            pairedness
            ))
        self.data.attrs['sample_f'] = sample_file
        self.data.attrs['pairedness'] = pairedness

    def _parse_vcf_file(self, vcf_file):
        sample_ids = self.data.attrs['sample_ids']
        sequence_ids = self.data.attrs['sequence_ids']
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            for seq_id in tqdm(sequence_ids, total=len(sequence_ids), desc="[%] Parsing input files ", ncols=100):
                zarr_key = ''
                for key, data in allel.read_vcf(vcf_file, region=seq_id, samples=sample_ids, fields=['samples', 'calldata/GT', 'variants/POS']).items():
                    if key == 'samples':
                        self.data.attrs['sample_ids_vcf'] = list(data)
                        #print("vcf samples", list(data))
                        self.data.attrs['sample_ids_to_vcf_idx'] = {sample_id: idx for idx, sample_id in enumerate(data)}
                    elif key == 'calldata/GT':
                        zarr_key = "%s/gt" % seq_id
                        #print(key, data.shape)
                        self.data.create_dataset(zarr_key, data=data)
                        
                    elif key == 'variants/POS':
                        zarr_key = "%s/pos" % seq_id
                        self.data.create_dataset(zarr_key, data=(np.array(data) - 1)) # port to BED (0-based) coordinates
                    else:
                        logging.error("[X] Unknown key %r" % key)
                        exit()
                unique_pos, counts_pos = np.unique(self.data["%s/pos" % seq_id], return_counts=True)
                duplicates = unique_pos[counts_pos > 1]
                if duplicates.any():
                    sys.exit("\n[X] Non-unique positions along sequence %r\n\tPositions: %s" % (seq_id, ", ".join([str(x) for x in list(duplicates)])))
            self.data.attrs['vcf_f'] = vcf_file

    def make_blocks(self, parameterObj):
        # if self.has_blocks() and not parameterObj.overwrite:
        #     sys.exit("[X] %s contains %s blocks (-l %s -m %s -r %s -u %s -i %s" % (
        #         store.   
        #         ))

        self.data.attrs['block_length'] = parameterObj.block_length
        self.data.attrs['block_gap_run'] = parameterObj.block_gap_run
        self.data.attrs['block_span'] = parameterObj.block_span
        self.data.attrs['mutypes_count'] = 4 # should be calculated from possible folded genotypes and pairedness

        sample_ids = self.data.attrs['sample_ids']
        sequence_ids = self.data.attrs['sequence_ids']
        sample_sets = self.data.attrs['sample_sets']
        #print(self.attrs())
        logging.info("[#] Processing BED file %r ..." % self.data.attrs['bed_f'])
        df = pd.read_csv(self.data.attrs['bed_f'], sep="\t", usecols=[0, 1, 2, 4], names=['sequence_id', 'start', 'end', 'samples'], 
            dtype={'sequence_id': str, 'start': np.int, 'end': np.int, 'samples': str})
        # remove sequence_ids that are not in sequence_names_array, sort, reset index
        intervals_df = df[df['sequence_id'].isin(sequence_ids)].sort_values(['sequence_id', 'start'], ascending=[True, True]).reset_index(drop=True)
        # get length column
        intervals_df['length'] = intervals_df['end'] - intervals_df['start'] 
        # get coverage matrix and drop columns of samples that are not in sample_ids_array
        intervals_df = pd.concat([intervals_df, intervals_df.samples.str.get_dummies(sep=',').filter(sample_ids)], axis=1).drop(columns=['samples'])
        # remove those intervals including less than two sample_ids (query sample_id columns : intervals_df[intervals_df.columns.intersection(sample_ids)])
        intervals_df = intervals_df.loc[(intervals_df[intervals_df.columns.intersection(sample_ids)].sum(axis=1) > 1)]
        interval_spans, block_spans = [], []
        with tqdm(total=(len(sequence_ids) * len(sample_sets)), desc="[%] Calculating bSFSs ", ncols=100, unit_scale=True) as pbar:
            for seq_id in sequence_ids:        
                interval_span, block_span = [], []
                _intervals_df = intervals_df[intervals_df['sequence_id'] == seq_id]
                _pos = self.data["%s/pos" % seq_id] # zarr.core.array
                genotypeArray = allel.GenotypeArray(self.data["%s/gt" % seq_id])
                for sample_set_idx, sample_set in enumerate(sample_sets):
                    sample_set_intervals_df = _intervals_df[_intervals_df[sample_set].all(axis='columns')]
                    interval_span.append(sample_set_intervals_df['length'].sum())
                    block_sites = cut_blocks(
                        sample_set_intervals_df.start, 
                        sample_set_intervals_df.end, 
                        parameterObj.block_length, 
                        parameterObj.block_span, 
                        parameterObj.block_gap_run
                        )
                    block_span.append(np.sum(((block_sites[:,-1] - block_sites[:,0]) + 1))) # block_space == span !
                    self.data.create_dataset("%s/%s/blocks/starts" % (seq_id, sample_set_idx), data=block_sites[:,0])
                    self.data.create_dataset("%s/%s/blocks/ends" % (seq_id, sample_set_idx), data=(block_sites[:,-1] + 1))  
                    idx_block_sites_in_pos = np.isin(block_sites, _pos, assume_unique=True) # will crash if non-unique pos
                    # idx_block_sites_in_pos = np.isin(block_sites, _pos) # will work even if non-unique pos
                    idx_pos_in_block_sites = np.isin(_pos, block_sites, assume_unique=True) # will crash if non-unique pos
                    # idx_pos_in_block_sites = np.isin(_pos, block_sites) # will work even if non-unique pos
                    # make room for reading in array
                    sample_set_idxs = [self.data.attrs['sample_ids_to_vcf_idx'][sample_id] for sample_id in sample_set] # list of indices
                    sample_set_genotypes = genotypeArray.subset(idx_pos_in_block_sites, sample_set_idxs)
                    block_sites[idx_block_sites_in_pos] = genotype_to_mutype_array(sample_set_genotypes)
                    block_sites[~idx_block_sites_in_pos] = 2 # monomorphic = 2 (0 = multiallelic, 1 = missing)
                    multiallelic, missing, monomorphic, variation = block_sites_to_variation_arrays(block_sites, self.data.attrs['mutypes_count'])

                    self.data.create_dataset("%s/%s/blocks/variation" % (seq_id, sample_set_idx), data=variation)
                    self.data.create_dataset("%s/%s/blocks/multiallelic" % (seq_id, sample_set_idx), data=multiallelic.flatten())
                    self.data.create_dataset("%s/%s/blocks/missing" % (seq_id, sample_set_idx), data=missing.flatten())
                    pbar.update(1)
                self.data.create_dataset("%s/interval_span" % seq_id, data=interval_span)
                self.data.create_dataset("%s/block_span" % seq_id, data=block_span)
                interval_spans.append(interval_span)
                block_spans.append(block_span)
        #print("interval_span", interval_span)
        interval_bases_df = pd.DataFrame(interval_span)
        #print("interval_bases_df", interval_bases_df)
        #print("block_span", block_span)
        blocked_bases_df = pd.DataFrame(block_span)
        #print("blocked_bases_df", blocked_bases_df)
        blocked_fraction_df = blocked_bases_df / interval_bases_df
        blocked_bases_df['sum'] = blocked_bases_df.sum(axis=1)
        blocked_fraction_df['mean'] = blocked_fraction_df.mean(axis=1)
        self.data.attrs['blocks'] = True

    def make_windows(self, parameterObj):
        sample_sets_idxs = self.data.attrs['idx_cartesian_sample_sets']
        with tqdm(total=(len(self.data.attrs['sequence_ids']) * len(sample_sets_idxs)), desc="[%] Making windows ", ncols=100, unit_scale=True) as pbar: 
            for seq_id in self.data.attrs['sequence_ids']: 
                mutypes, sample_set_covs, starts, ends = [], [], [], []
                for sample_set_idx in sample_sets_idxs:
                    multiallelic = np.array(self.data["%s/%s/blocks/multiallelic" % (seq_id, sample_set_idx)])
                    missing = np.array(self.data["%s/%s/blocks/missing" % (seq_id, sample_set_idx)])
                    valid = np.less_equal(missing, parameterObj.max_missing) & np.less_equal(multiallelic, parameterObj.max_multiallelic)
                    mutype = np.array(self.data["%s/%s/blocks/variation" % (seq_id, sample_set_idx)])
                    mutypes.append(mutype[valid])
                    start = np.array(self.data["%s/%s/blocks/starts" % (seq_id, sample_set_idx)])
                    starts.append(start[valid])
                    end = np.array(self.data["%s/%s/blocks/ends" % (seq_id, sample_set_idx)])
                    ends.append(end[valid])
                    sample_set_covs.append(np.full_like(end[valid], sample_set_idx)) 
                    pbar.update()
                mutype_array = np.concatenate(mutypes, axis=0)
                start_array = np.concatenate(starts, axis=0)
                end_array = np.concatenate(ends, axis=0)
                sample_set_array = np.concatenate(sample_set_covs, axis=0)
                window_mutypes, window_min, window_max, window_mean, window_median = cut_windows(mutype_array, sample_sets_idxs, start_array, end_array, sample_set_array, num_blocks=parameterObj.window_size, num_steps=parameterObj.window_step)
                window_id = np.array(["_".join([seq_id, _start, _end]) for (_start, _end) in zip(window_min.astype(str), window_max.astype(str))])
                self.data.create_dataset("%s/windows/variation" % seq_id, data=window_mutypes)
                self.data.create_dataset("%s/windows/window_id" % seq_id, data=window_id)
                self.data.create_dataset("%s/windows/midpoint_mean" % seq_id, data=window_mean)
                self.data.create_dataset("%s/windows/midpoint_median" % seq_id, data=window_median)
        self.data.attrs['windows'] = True

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
        sample_sets_idxs = self.data.attrs['idx_cartesian_sample_sets']
        starts, ends, mutypes, sample_set = [], [], [], []
        with tqdm(total=(len(self.data.attrs['sequence_ids']) * len(sample_sets_idxs)), desc="[%] Writing bSFSs ", ncols=100, unit_scale=True) as pbar: 
            for seq_id in self.data.attrs['sequence_ids']: 
                for sample_set_idx in sample_sets_idxs:
                    multiallelic = np.array(self.data["%s/%s/blocks/multiallelic" % (seq_id, sample_set_idx)])
                    missing = np.array(self.data["%s/%s/blocks/missing" % (seq_id, sample_set_idx)])
                    valid = np.less_equal(missing, parameterObj.max_missing) & np.less_equal(multiallelic, parameterObj.max_multiallelic)
                    mutype = np.array(self.data["%s/%s/blocks/variation" % (seq_id, sample_set_idx)])
                    mutypes.append(mutype[valid])
                    start = np.array(self.data["%s/%s/blocks/starts" % (seq_id, sample_set_idx)])
                    starts.append(start[valid])
                    end = np.array(self.data["%s/%s/blocks/ends" % (seq_id, sample_set_idx)])
                    sample_set.append(np.full_like(end[valid], sample_set_idx)) 
                    ends.append(end[valid])
                    pbar.update()
        mutypes_array = np.concatenate(mutypes, axis=0)
        #print("mutypes_array", mutypes_array)
        # popgen
        pi_1, pi_2, d_xy, f_st, fgv = calculate_popgen_from_array(mutypes_array, (self.data.attrs['block_length'] * mutypes_array.shape[0]))
        print("[+] Pi_%s = %s; Pi_%s = %s; D_xy = %s; F_st = %s; FGVs = %s / %s blocks (%s)" % (self.data.attrs['pop_ids'][0], pi_1, self.data.attrs['pop_ids'][1], pi_2, d_xy, f_st, fgv, mutypes_array.shape[0], format_percentage(fgv / mutypes_array.shape[0]))) 
        
        # mutuple barchart
        mutypes, counts = np.unique(mutypes_array, return_counts=True, axis=0)
        mutype_counter = collections.Counter({tuple(i):j for i,j in zip(mutypes, counts)})
        plot_mutuple_barchart('%s.mutuple_barchart.png' % self.prefix, mutype_counter)

        # mutuple tally
        bsfs = np.concatenate([counts[:, np.newaxis], mutypes], axis =-1)
        header = ['count'] + [x+1 for x in range(self.data.attrs['mutypes_count'])]
        pd.DataFrame(data=bsfs, columns=header, dtype='int64').to_hdf("%s.blocks.h5" % self.prefix, 'bsfs', format='table')

        # block coordinates
        #header = ['block_id', 'start', 'end', 'sample_set', 'multiallelic', 'missing']
        #pd.DataFrame(data=bsfs, columns=header, dtype='int64').to_hdf("%s.blocks.h5" % self.prefix, 'bsfs', format='table')

    def dump_windows(self, parameterObj):
        window_info_rows = []
        window_mutuple_tally = []
        for sequence_id in tqdm(self.data.attrs['sequence_ids'], total=len(self.data.attrs['sequence_ids']), desc="[%] Generating output ", ncols=100):
            variations = self.data["%s/windows/variation" % sequence_id]
            window_ids = self.data["%s/windows/window_id" % sequence_id]
            midpoint_means = self.data["%s/windows/midpoint_mean" % sequence_id]
            midpoint_medians = self.data["%s/windows/midpoint_median" % sequence_id]
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
        window_bsfs_df.to_hdf("%s.window_bsfs.h5" % self.prefix, 'bsfs', format='table')
        
        window_info_cols = ['window_id', 'sequence_id', 'midpoint_mean', 'midpoint_median', 'pi_%s' % self.data.attrs['pop_ids'][0], 'pi_%s' % self.data.attrs['pop_ids'][1], 'd_xy', 'f_st', 'fgv']
        window_info_df = pd.DataFrame(window_info_rows, columns=window_info_cols)
        window_info_df.to_hdf("%s.windows_info.h5" % self.prefix, 'info', format='table')
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

    def _block_stats(self):
        pass
        # logging.info(f"[+] {block_spaces[-1]:,} b in blocks ({(block_spaces[-1] / interval_spaces[-1]):.1%} of total bases in intervals)")