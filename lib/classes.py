import collections
import itertools
from tqdm import tqdm
from lib.functions import get_pop_ids_by_sample_ids_from_csv, plot_mutuple_barchart, format_bases, format_fraction, create_hdf5_store, tabulate_df, poolcontext, format_percentage, plot_distance_scatter, plot_fst_genome_scan, plot_pi_genome_scan, plot_pi_scatter, plot_sample_barchart
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
import sys
from timeit import default_timer as timer
'''
[To do]

# unify output as hd5
- Make module to create datastore
- Add blocks
- Add variants
- Add filters 
- Add windows




# filtering
- variants:
    plot distribution of full_mutypes

- filtering DSL require -f 'missing <= 1'
    >>> a = 'missing <= 1'
    >>> a
    'missing <= 1'
    >>> a.split("&&")
    ['missing <= 1']
    >>> a.split("=>")
    ['missing <= 1']
    operators = ["==", '<=', '>=', '!=']
    comparators = ['eq', 'le', 'ge', 'ne']
    for operator, comparator in zip(operators, comparators):

    >>> mutype, value = a.split("<=")
    ['missing ', ' 1']
    >>> a.split("<=")
    label = 'le'
    operator = "<="
    value = 
'''
COLOURS = ['orange', 'dodgerblue']
FULL_MUTYPE_ORDER = ['hetA', 'fixed', 'hetB', 'hetAB', 'missing', 'multiallelic']
MUTYPE_ORDER = ['hetA', 'fixed', 'hetB', 'hetAB']
MUTYPE_OTHER = ['missing', 'multiallelic']
GT_ORDER = ['TOTAL', 'MISS', 'HOM', 'HET']

'''
# setup 
    - reads vcf into store
    - reads bed_f into store
    - generates pair_ids

# blocks
    - reads bed_f and makes blocks for each pair

# variants
    - queries variants for all bed regions for each pair

'''
store_meta_by_key = {
    'sequence_ids': '/meta/sequence_ids',
    'sequence_lengths': '/meta/sequence_lengths',
    'sample_ids': '/meta/sample_ids',
    'population_ids': '/meta/population_ids'
}

def data_key(key):
    if key == 'gt':
        pass
         
store_data_by_key = {
    'pos': '/data/variants/POS',
    'gt': '/data/calldata/GT'
}

np.set_printoptions(threshold=20)


import matplotlib.pyplot as plt
import numpy as np

def plot_histogram(x, out_f):
    fig, ax = plt.subplots(figsize=(14, 5))
    hist, bins = np.histogram(x, bins=50)
    width = 0.7 * (bins[1] - bins[0])
    center = (bins[:-1] + bins[1:]) / 2
    ax.bar(center, hist, align='center', width=width)
    fig.savefig('%s.png' % out_f, format="png")

def is_hom(gts):
    # print("gts", gts.shape, type(gts), gts[:2])                
    allele1 = gts[..., 0, np.newaxis]
    # print("allele1", allele1.shape, type(allele1), allele1[:2])
    other_alleles = gts[..., 1:]
    # print("other_alleles", other_alleles.shape, type(other_alleles), other_alleles[:2])
    tmp = (allele1 >= 0) & (allele1 == other_alleles)
    out = np.all(tmp, axis=-1)
    # print("out", out.shape, type(out), out[:2])
    return out

def is_hom2(gts):
    print("gts", gts.shape, type(gts), gts[:1])                
    allele1 = gts[..., 0, np.newaxis]
    print("allele1", allele1.shape, type(allele1), allele1[:1])
    other_alleles = gts[..., 1:]
    print("other_alleles", other_alleles.shape, type(other_alleles), other_alleles[:1])
    tmp = (allele1 >= 0) & (allele1 == other_alleles)
    print("tmp", tmp.shape, type(tmp))
    out = np.all(tmp, axis=-1)
    print("out", out.shape, type(out), out[:1])
    return out

def consecutive(data, stepsize=1):
    return np.split(data, np.where(np.diff(data) != stepsize)[0]+1)

def create_ranges(aranges):
    # does the transformation form 0-based (BED) to 1-based (VCF) coordinate system 
    # https://stackoverflow.com/a/47126435
    l = aranges[:,1] - aranges[:,0]
    clens = l.cumsum()
    ids = np.ones(clens[-1], dtype=int)
    ids[0] = aranges[0,0]
    ids[clens[:-1]] = aranges[1:, 0] - aranges[:-1, 1]+1
    return ids.cumsum()

def bool_to_unsigned_int(boolean):
    return np.sum(2**np.arange(boolean.shape[1])*boolean, axis=-1)

def make_blocks(sample_set_intervals_df, block_length, block_span, block_gap_run):
    sites = create_ranges(np.array((sample_set_intervals_df.start, sample_set_intervals_df.end)).T)
    return block(sites, block_length, block_span, block_gap_run)

def block(sites, block_length, block_span, block_gap_run):
    block_sites = np.concatenate([x[:block_length * (x.shape[0] // block_length)].reshape(-1, block_length) for x in np.split(sites, np.where(np.diff(sites) > block_gap_run)[0]+1)])
    return block_sites[(((block_sites[:,-1] - block_sites[:,0]) + 1) <= block_span)] # return blocks where ((block_end - block_start) + 1) <= block_span

class Store(object):
    def __init__(self, parameterObj):
        self.path = None
        self.data = None
        if parameterObj.module == 'setup':
            self._from_input(parameterObj)
        else:
            self._from_zarr(parameterObj)

    def _from_input(self, parameterObj):
        self.path = self._get_path(parameterObj.outprefix)
        self.data = zarr.open(self.path, mode='w')
        self._parse_genome_file(parameterObj.genome_file)
        self._parse_sample_file(
            parameterObj.sample_file, 
            parameterObj.pairedness
            )
        self._parse_vcf_file(parameterObj.vcf_file)
        self._parse_bed_file(
            parameterObj.bed_file, 
            parameterObj.block_length, 
            parameterObj.block_span, 
            parameterObj.block_gap_run
        )

    def _from_zarr(self, parameterObj):
        self.path = parameterObj.zstore
        self.data = zarr.open(self.path, mode='r')

    def tree(self):
        return self.data.tree()

    def attrs(self):
        return "\n".join(["\t".join([k, str(len(v)), str(type(v)), str(v)]) for k, v in self.data.attrs.asdict().items()])

    def _parse_genome_file(self, genome_file):
        logging.info("[#] Processing Genome file %r..." % genome_file)
        df = pd.read_csv(genome_file, sep="\t", names=['sequence_id', 'sequence_length'], header=None)
        self.data.attrs['sequence_ids'] = df['sequence_id'].to_list() #df['sequence_id'].to_numpy(dtype=str)
        self.data.attrs['sequence_length'] = df['sequence_length'].to_list()
        logging.info("[+] Found %s sequences of a total length of %s b..." % (len(self.data.attrs['sequence_ids']), sum(self.data.attrs['sequence_length'])))
        for sequence_id in self.data.attrs['sequence_ids']:
            self.data.create_group('%s/' % sequence_id)
        
    def _parse_sample_file(self, sample_file, pairedness):
        logging.info("[#] Processing Sample file %r ..." % sample_file)
        df = pd.read_csv(sample_file, sep=",", names=['sample_id', 'population_id'])
        self.data.attrs['sample_ids'] = df['sample_id'].to_list() #.to_numpy(dtype=str)
        self.data.attrs['population_ids'] = df['population_id'].to_list() #.to_numpy(dtype=str)  
        # self.data.attrs['sample_sets'] = ["_".join(_) for _ in itertools.combinations(sorted(self.data.attrs['sample_ids']), pairedness)]
        self.data.attrs['sample_sets'] = [_ for _ in itertools.combinations(sorted(self.data.attrs['sample_ids']), pairedness)]
        logging.info("[+] Found %s samples from %s populations. Generated %s sets of %s samples" % (
            len(self.data.attrs['sample_ids']), 
            len(self.data.attrs['population_ids']),
            len(self.data.attrs['sample_sets']),
            pairedness
            ))

    def get_gts(self, sequence_ids, start, end, sample_ids):
        pass

    def _parse_bed_file(self, bed_file, block_length, block_span, block_gap_run):
        '''
        https://stupidpythonideas.blogspot.com/2014/01/grouping-into-runs-of-adjacent-values.html
        https://stackoverflow.com/questions/2154249/identify-groups-of-continuous-numbers-in-a-list
        https://stackoverflow.com/questions/4494404/find-large-number-of-consecutive-values-fulfilling-condition-in-a-numpy-array
        https://stackoverflow.com/questions/32357823/splitting-a-list-of-integers-by-breaks-in-the-data-when-the-break-size-is-variab
        '''
        sample_ids = self.data.attrs['sample_ids']
        sequence_ids = self.data.attrs['sequence_ids']
        sample_sets = self.data.attrs['sample_sets']

        logging.info("[#] Processing BED file %r ..." % bed_file)
        df = pd.read_csv(bed_file, sep="\t", usecols=[0, 1, 2, 4], names=['sequence_id', 'start', 'end', 'samples'], 
            dtype={'sequence_id': str, 'start': np.int, 'end': np.int, 'samples': str})
        # remove sequence_ids that are not in sequence_names_array, sort, reset index
        intervals_df = df[df['sequence_id'].isin(sequence_ids)].sort_values(['sequence_id', 'start'], ascending=[True, True]).reset_index(drop=True)
        # get length column
        intervals_df['length'] = intervals_df['end'] - intervals_df['start'] 
        # get distance column
        # intervals_df['distance'] = np.where((intervals_df['sequence_id'] == intervals_df['sequence_id'].shift(-1)), intervals_df['start'].shift(-1) - intervals_df['end'], block_gap_run + 1).astype(np.int64)
        # get coverage matrix and drop columns of samples that are not in sample_ids_array
        intervals_df = pd.concat([intervals_df, intervals_df.samples.str.get_dummies(sep=',').filter(sample_ids)], axis=1).drop(columns=['samples'])
        #print(intervals_df)
        # remove those intervals including less than two sample_ids (query sample_id columns : intervals_df[intervals_df.columns.intersection(sample_ids)])
        intervals_df = intervals_df.loc[(intervals_df[intervals_df.columns.intersection(sample_ids)].sum(axis=1) > 1)]
        with tqdm(total=(len(sequence_ids) * len(sample_sets)), desc="[%] ", ncols=100, unit_scale=True) as pbar:
            for seq_id in sequence_ids:
                _intervals_df = intervals_df[intervals_df['sequence_id'] == seq_id]
                #print(_intervals_df)
                
                interval_spaces = []
                block_spaces = []
                
                _pos = self.data["%s/pos" % seq_id] # zarr.core.array
                genotypeArray = allel.GenotypeArray(self.data["%s/gt" % seq_id])
                #print(genotypeArray)
                #print(len(sample_sets))
                #for sample_set in tqdm(sample_sets, total=len(sample_sets), desc="[%] ", ncols=100):
                for sample_set in sample_sets:
                    #print(sample_set)
                    # blocks
                    sample_set_intervals_df = _intervals_df[_intervals_df[sample_set].all(axis='columns')]
                    interval_spaces.append(sample_set_intervals_df['length'].sum())
                    block_sites = make_blocks(sample_set_intervals_df, block_length, block_span, block_gap_run)
                    block_spaces.append(np.sum(((block_sites[:,-1] - block_sites[:,0]) + 1))) # block_space == span !                
                    # logging.info(f"[+] {block_spaces[-1]:,} b in blocks ({(block_spaces[-1] / interval_spaces[-1]):.1%} of total bases in intervals)")
    
                    # variants
                    sample_selection = [self.data.attrs['sample_ids_to_vcf_idx'][sample_id] for sample_id in sample_set]
                    print(sample_set, sample_selection)
                    sample_set_genotype = genotypeArray.subset(None, sample_selection)
                    #print(sample_set_genotype)
                    #print(sample_set_genotype.shape)
                    pbar.update(1)
                self.data.create_dataset("%s/raw_sites" % seq_id, data=interval_spaces)
                self.data.create_dataset("%s/block_sites" % seq_id, data=block_spaces)
            
        #### 

        #    hom_gt = is_hom(_gt)
        #    #print(hom_gt[:10], hom_gt.shape)
            
        #    #print(miss_gt[:10], miss_gt.shape)
        #    for pair_tuple in itertools.combinations(sorted(sample_ids), 2):
        #        a_sample_id = pair_tuple[0]
        #        a_vcf_idx = self.data.attrs['sample_ids_to_vcf_idx'][pair_tuple[0]]
        #        b_sample_id = pair_tuple[1]
        #        b_vcf_idx = self.data.attrs['sample_ids_to_vcf_idx'][pair_tuple[1]]
        #        pair_id = "%s_%s" % (a_sample_id, b_sample_id)
        #        pair_intervals_df = _intervals_df[_intervals_df[[a_sample_id, b_sample_id]].all(axis='columns')]
        #        block_sites = make_blocks(pair_intervals_df, max_gap=1, block_length=5)
        #        #print(block_sites, block_sites.shape)
        #        block_starts = block_sites[:, 0]
        #        self.data.create_dataset("%s/blocks/%s/block_starts/" % (seq_id, pair_id), data=block_starts)
        #        #print('block_starts', block_starts, block_starts.shape)
        #        block_ends = block_sites[:, -1]
        #        self.data.create_dataset("%s/blocks/%s/block_ends/" % (seq_id, pair_id), data=block_ends)
        #        #print('block_ends', block_ends, block_ends.shape)
        #        #print("----")
#        #        is_MISS = np.logical_or(miss_gt[:, a_vcf_idx], miss_gt[:, b_vcf_idx])
        #        #print('is_MISS', is_MISS, is_MISS.shape)
        #        is_SAME = np.all(_gt[:, a_vcf_idx] == _gt[:, b_vcf_idx], axis=-1)
        #        #print('is_SAME', is_SAME, is_SAME.shape)
        #        is_a_HOM = hom_gt[:, a_vcf_idx]
        #        #print('is_a_HOM', is_a_HOM, is_a_HOM.shape)
        #        is_b_HOM = hom_gt[:, b_vcf_idx]
        #        #print('is_b_HOM', is_b_HOM, is_b_HOM.shape)
        #        boolean = np.concatenate([np.transpose([is_a_HOM]), np.transpose([is_b_HOM]), np.transpose([is_SAME])], axis=1)
        #        #print('boolean', boolean, boolean.shape)
        #        unsigned_int = bool_to_unsigned_int(boolean)
        #        #print('unsigned_int', unsigned_int, unsigned_int.shape)
        #        #test = np.where(np.in1d(block_sites.ravel(),_pos), unsigned_int, block_sites.ravel())
        #        #print('np.in1d(block_sites, _pos)', np.in1d(block_sites, _pos), np.in1d(block_sites, _pos).shape)
        #        #print('np.in1d(_pos, block_sites)', np.in1d(_pos, block_sites), np.in1d(_pos, block_sites).shape)
        #        # np.in1d(block_sites, _pos) [False False False ... False False False] (29056,)
        #        # np.in1d(_pos, block_sites) [False False False ... False False False] (2648,)
        #        block_sites.ravel()[~np.in1d(block_sites, _pos)] = 7
        #        block_sites.ravel()[np.in1d(block_sites, _pos)] = unsigned_int[np.in1d(_pos, block_sites)]
        #        uniques, unique_counts = np.unique(block_sites, return_counts=True, axis = 0)
        #        print(pair_id, "Counts")
        #        print(uniques, uniques.shape)
        #        print(unique_counts, unique_counts.shape)
        #        def count_unique_by_row(a):
        #            weight = 1j*np.linspace(0, a.shape[1], a.shape[0], endpoint=False)
        #            b = a + weight[:, np.newaxis]
        #            u, ind, cnt = np.unique(b, return_index=True, return_counts=True)
        #            b = np.zeros_like(a)
        #            np.put(b, ind, cnt)
        #            return b
#        #        '''
        #        >>> x
        #        array([[0, 0, 0, 1, 1],
        #               [0, 1, 0, 0, 1],
        #               [0, 1, 1, 0, 0],
        #               [1, 0, 0, 0, 1],
        #               [0, 1, 1, 1, 1]])
        #        >>> c = count_unique_by_row(x)
        #        array([[3, 0, 0, 2, 0],
        #               [3, 2, 0, 0, 0],
        #               [3, 2, 0, 0, 0],
        #               [2, 3, 0, 0, 0],
        #               [1, 4, 0, 0, 0]])
        #        np.concatenate([x[np.nonzero(count_unique_by_row(x))], x)
        #        c.non_zero_count()
#
        #        '''
                #print('block_sites', block_sites, block_sites.shape)
                # print('_pos', _pos, _pos.shape)
                
                # mask_pos = np.in1d(block_sites, _pos)
                # print('mask_pos', mask_pos, mask_pos.shape)
                # print('block_sites.ravel()[mask_pos]', block_sites.ravel()[mask_pos], block_sites.ravel()[mask_pos].shape)
                # block_sites[mask_pos] = unsigned_int[mask_pos]
                
                
                
                #print("----")
                


                # import numpy as np
                # gt = np.random.randint(-1, 2, (2648, 7, 2))
                # pos = np.sort(np.random.choice(10000, 2648, replace=False))
                # interval_range = np.sort(np.random.choice(10000, 6000, replace=False))
                # def make_blocks(data, max_gap=5, block_length=10):
                #    return np.concatenate([x[:block_length * (x.shape[0] // block_length)].reshape(-1, block_length) for x in np.split(data, np.where(np.diff(data) > max_gap)[0]+1)])
                # blocks = make_blocks(interval_range, max_gap=5, block_length=10)
                # block_span = blocks[:, -1] - blocks[:, 0]
                # blocks.ravel()[np.in1d(blocks,pos)] = 0
                # access sample genotypes : gt[np.in1d(pos,blocks),2]
                
                # identical genotypes : np.all(gt[:,1] == gt[:,5], axis=-1)[:10]
                # set blocks to blocks.ravel()[np.in1d(blocks,pos)] = ... 
                
                # HOM_MATRIX = is_hom(gt)
                # MISSING_MATRIX = np.any(gt < 0, axis=-1)
                # a, b = pair_idx
                # is_MISS = np.logical_or(MISSING_MATRIX[:,a], MISSING_MATRIX[:,b])
                # is_SAME = np.all(gt[:,a] == gt[:,b], axis=-1)
                # is_a_HOM = HOM_MATRIX[a]
                # is_b_HOM = HOM_MATRIX[b]




                #self.data.create_dataset(store_meta_by_key['coverages'], data=coverages)
        #    _coverages_df = coverages_df[coverages_df['sequence_id'] == seq_id]
         #   print(_coverages_df)
            #coverages_df.loc[mask].reset_index(drop=True).to_numpy()
            #
            #print(coverages)
            #intervals = intervals.loc[mask].reset_index(drop=True).to_numpy()
            #print(intervals)
            #self.data.create_dataset(store_meta_by_key['intervals'], data=intervals)
        #store.create_dataset('data/coverages', data=coverages)

        #intervals, coverages = parse_intervals(parameterObj.bed_file, sequence_names, sample_ids)
        ##store.create_dataset('data/intervals', data=intervals)
        #store.create_dataset('data/coverages', data=coverages)

    def _parse_vcf_file(self, vcf_file):
        sample_ids = self.data.attrs['sample_ids']
        sequence_ids = self.data.attrs['sequence_ids']
        for seq_id in tqdm(sequence_ids, total=len(sequence_ids), desc="[%] ", ncols=100):
            zarr_key = ''
            for key, data in allel.read_vcf(vcf_file, region=seq_id, samples=sample_ids, fields=['samples', 'calldata/GT', 'variants/POS']).items():
                if key == 'samples':
                    self.data.attrs['sample_ids_vcf'] = list(data)
                    self.data.attrs['sample_ids_to_vcf_idx'] = {sample_id: idx for idx, sample_id in enumerate(data)}
                elif key == 'calldata/GT':
                    zarr_key = "%s/gt" % seq_id
                    #print(key, data.shape)
                    self.data.create_dataset(zarr_key, data=data)
                elif key == 'variants/POS':
                    zarr_key = "%s/pos" % seq_id
                    self.data.create_dataset(zarr_key, data=data)
                else:
                    logging.error("[X] Unknown key %r" % key)
                    exit()

    def _get_path(self, outprefix):
        path = "%s.z" % outprefix
        if os.path.isdir(path):
            logging.info("[!] ZARR store %r exists. Deleting ..." % path)
            shutil.rmtree(path)
        logging.info("[+] Generating ZARR store %r" % path)
        return path

def parse_genome_file(genome_file):
    logging.info("[#] Processing Genome file %r..." % genome_file)
    df = pd.read_csv(genome_file, sep="\t", names=['sequence_id', 'sequence_length'], header=None)
    logging.info("[+] Found %s sequences" % len(df.index))
    return df['sequence_id'].to_numpy(dtype=str), df['sequence_length'].to_numpy()

def parse_sample_file(sample_file):
    logging.info("[#] Processing Sample file %r ..." % sample_file)
    df = pd.read_csv(sample_file, sep=",", names=['sample_id', 'population_id'])
    logging.info("[+] Found %s samples from %s populations" % (len(df.index), len(df.population_id.unique())))
    return df['sample_id'].to_numpy(dtype=str), df['population_id'].to_numpy(dtype=str)  

def parse_intervals(bed_file, sequence_names_array, sample_ids_array):
    logging.info("[#] Processing BED file %r ..." % bed_file)
    df = pd.read_csv(bed_file, sep="\t", usecols=[0, 1, 2, 4], names=['sequence_id', 'start', 'end', 'samples'], 
        dtype={'sequence_id': str, 'start': np.int, 'end': np.int, 'samples': str})
    # remove sequence_ids that are not in sequence_names_array, sort, reset index
    df = df[df['sequence_id'].isin(sequence_names_array)].sort_values(['sequence_id', 'start'], ascending=[True, True]).reset_index(drop=True)
    # get coverage matrix and drop columns of samples that are not in sample_ids_array
    coverages = df.samples.str.get_dummies(sep=',').filter(sample_ids_array).reset_index(drop=True)
    print(coverages)
    # drop samples column (no longer needed)
    intervals = df.drop(columns=['samples'])
    # get length column
    intervals['length'] =  intervals['end'] - intervals['start'] 
    # remove all intervals/coverages for which the filtered coverages are all zero
    mask = (coverages!=0).any(1)
    coverages = coverages.loc[mask].reset_index(drop=True)
    print(coverages)
    intervals = intervals.loc[mask].reset_index(drop=True)
    print(intervals)
    return intervals.to_numpy(), coverages.to_numpy()

def load_StoreObj(parameterObj):
    #storeObj = StoreObj()

    return StoreObj()

class StoreObj(object):
    def __init__(self):
        self.zarr_dir = None
        self.dataset = None
        self.bed_f = None
        self.sample_ids_by_pop_id = {} # 'pop1', 'pop2', ...
        self.pop_ids_order = []
        #self.outprefix = parameterObj.outprefix
        #if parameterObj.zarr_dir:
        #    self._load_store(parameterObj)
        #if parameterObj.vcf_file:
        #    self._initiate_store(parameterObj)
        #if parameterObj.bed_file:
        #    self._parse_intervals(parameterObj)

    def parse_sample_file(self, parameterObj):
        # parse sample CSV
        samples_df = pd.read_csv(\
            parameterObj.sample_file, \
            sep=",", \
            usecols=[0, 1], \
            names=['sample_id', 'population_id'], \
            header=None, \
            dtype={ \
                'sample_id': 'category', \
                'population_id': 'category' \
                } \
            )
        # Error if not 2 populations
        if not len(samples_df.groupby('population_id').count().index) == 2:
             exit('[X] Invalid number of populations: %s (must be 2)' % len(samples_df.groupby('population_id').count().index))
        sorted_samples_df = samples_df.sort_values(['population_id', 'sample_id'])
        print(sorted_samples_df)
        population_idx = -1
        sample_ids_by_population_id = collections.defaultdict(list)
        for sample_idx, (sample_id, population_id) in enumerate(sorted_samples_df.values.tolist()):
            if not population_id in self.population_idx_by_population_id:
                population_idx += 1
                self.add_population(population_idx, population_id)
            self.add_sample_to_populationObj(population_id, sample_id)
            sample_ids_by_population_id[population_idx].append(sample_id)
        for pair_idx, pair_id in enumerate([(x) for x in itertools.product(*sorted(sample_ids_by_population_id.values()))]):
            self.add_pair(pair_idx, pair_id)

    def parse_genome_file(self, parameterObj):
        df = pd.read_csv(parameterObj.genome_file, sep="\t", usecols=[0, 1], names=['sequence_id', 'length'], header=None)
        for sequence_idx, (sequence_id, length_str) in enumerate(df.values.tolist()):
            try:
                length = int(length_str)
            except TypeError:
                exit("[X] Line %s: Second column of --genome_file %s must be integers, not '%s'" % (sequence_idx, parameterObj.genome_file, length_str))
            self.add_sequence(sequence_idx, sequence_id, length)


    def tree(self):
        return self.dataset.tree()

    def _initiate_store(self, parameterObj):
        self.zarr_dir = str(pathlib.Path(parameterObj.outprefix).with_suffix('.zarr'))
        self.dataset = self._get_dataset(overwrite=True)
        self._parse_data(parameterObj)
        self.sample_ids_by_pop_id = self._get_sample_ids_by_pop_id()
        self.pop_ids_order = sorted(set([pop_id for pop_id in self.yield_pop_ids()]))

    def _parse_intervals(self, parameterObj):
        logging.info("[#] Processing BED file ...")
        print(parameterObj.bed_file)
        bed_df = pd.read_csv( 
            parameterObj.bed_file, 
            sep="\t", 
            usecols=[0, 1, 2, 4], 
            names=['sequence_id', 'start', 'end', 'samples'], 
            skiprows=1, 
            header=None, 
            dtype={ 
                'sequence_id': 'category', 
                'start': np.int, 
                'end': np.int, 
                'samples': 'category'})
        # MIGHT NOT BE NECESSARY
        # filter rows based on sequence_ids, sort by sequence_id, start
        # bed_df = bed_df[bed_df['sequence_id'].isin(entityCollection.sequence_idx_by_sequence_id)].sort_values(['sequence_id', 'start'], ascending=[True, True])
        # compute length 
        bed_df['length'] =  bed_df['end'] - bed_df['start'] 
        print(bed_df)
        bed_length_total = int(bed_df['length'].sum())
        print(bed_length_total)
        # compute cov-matrix (sample_ids are column names)
        cov_df = bed_df.samples.str.get_dummies(sep=',')
        # for each combination of sample_ids ...
        print(self.sample_ids_by_pop_id)
        for a, b in itertools.combinations(self.sample_ids_by_pop_id['all'], 2):
            pair_id = frozenset(a, b)
            print(pair_id)
            boolean_mask = cov_df[[a, b]].all(axis='columns')
            interval_space = interval_df[boolean_mask]
            block_arrays = np.split(interval_space, block_length)

        # https://stupidpythonideas.blogspot.com/2014/01/grouping-into-runs-of-adjacent-values.html
        # https://stackoverflow.com/questions/2154249/identify-groups-of-continuous-numbers-in-a-list
        # https://stackoverflow.com/questions/4494404/find-large-number-of-consecutive-values-fulfilling-condition-in-a-numpy-array
        a = np.array((bed_df.start, bed_df.end)).T
        r = create_ranges(a)
        def consecutive(data, stepsize=1):
            return np.split(data, np.where(np.diff(data) != stepsize)[0]+1)
        a = np.array([0, 47, 48, 49, 50, 97, 98, 99])
        consecutive(a)

    def _load_store(self, parameterObj):
        self.zarr_dir = parameterObj.zarr_dir
        self._add_sample_information(parameterObj)
        self.dataset = self._get_dataset()
        self.sample_ids_by_pop_id = self._get_sample_ids_by_pop_id()
        self.pop_ids_order = sorted(set([pop_id for pop_id in self.yield_pop_ids()]))

    def yield_pop_ids(self):
        self._yield('/pop_ids')

    def yield_sample_ids(self):
        self._yield('/pop_ids')

    def yield_seq_ids(self):
        self._yield('/seq_ids')

    def yield_seq_lengths(self):
        self._yield('/seq_lengths')

    def _yield(self, key):
        for value in list(self.dataset[key]):
            yield value

    def _get_sample_ids_by_pop_id(self):
        sample_ids_by_pop_id = collections.defaultdict(list)
        sample_ids_by_pop_id['all'] = list(range(len(self.sample_ids)))
        for sample_id, pop_id in zip(range(len(self.sample_ids)), self.pop_ids):
            sample_ids_by_pop_id[pop_id].append(sample_id)
        return sample_ids_by_pop_id

    def _get_dataset(self, overwrite=False):
        if overwrite:
            if os.path.isdir(self.zarr_dir):
                logging.info("[!] ZARR store %r exists. Deleting ..." % self.zarr_dir)
                shutil.rmtree(self.zarr_dir)
            logging.info("[+] Generating ZARR store %r" % self.zarr_dir)
            return zarr.open(self.zarr_dir, mode='w')
        return zarr.open(self.zarr_dir, mode='a')

    def _parse_data(self, parameterObj):
        self.sample_ids = self._add_sample_ids(allel.read_vcf_headers(parameterObj.vcf_file).samples)
        seq_ids = self._add_sequences_from_header(allel.read_vcf_headers(parameterObj.vcf_file).headers)
        self._add_variants(parameterObj.vcf_file, seq_ids)
        self._add_sample_information(parameterObj, self.sample_ids)

    def _add_sample_ids(self, sample_ids):
        self.dataset.create_dataset('sample_ids', data=np.array(sample_ids))
        return sample_ids

    def _add_sequences_from_header(self, header_lines):
        pattern = re.compile(r'##contig=<ID=(\S+),length=(\d+)>')
        seq_ids, lengths = [], []
        for header_line in header_lines:
            if header_line.startswith("##contig"):
                seq_id, length = re.match(pattern, header_line).groups()
                seq_ids.append(seq_id)
                lengths.append(int(length))
        self.dataset.create_dataset('seq_ids', data=np.array(seq_ids))
        self.dataset.create_dataset('seq_lengths', data=np.array(lengths))
        return seq_ids

    def _add_variants(self, vcf_file, seq_ids):
        for seq_id in tqdm(seq_ids, total=len(seq_ids), desc="[%] ", ncols=100):
            allel.vcf_to_zarr( 
                vcf_file, 
                self.zarr_dir, 
                group="%s" % seq_id,
                region=seq_id,
                fields=[
                    'variants/POS',
                    'variants/REF',
                    'variants/ALT',
                    'variants/QUAL',
                    'variants/is_snp',
                    'variants/numalt',
                    'variants/DP',
                    'calldata/GT',
                    'calldata/DP'
                ], 
                overwrite=True)

    def _add_sample_information(self, parameterObj, sample_ids):
        pop_id_by_sample_id = get_pop_ids_by_sample_ids_from_csv(parameterObj)
        self.dataset.create_dataset('pop_ids', data=np.array([pop_id_by_sample_id[sample_id] for sample_id in sample_ids]))

def mutype_counter_to_mutuple(mutype_counter, full=True):
    if full:
        return tuple([mutype_counter[mutype] for mutype in FULL_MUTYPE_ORDER])
    return tuple([mutype_counter[mutype] for mutype in MUTYPE_ORDER])

def mutype_counter_to_dict(mutype_counter, full=True):
    order = MUTYPE_ORDER
    if full:
        order = FULL_MUTYPE_ORDER
    return dict({mutype: mutype_counter[mutype] for mutype in order})

def gt_counter_to_dict(gt_counter): 
    return dict({gt: gt_counter[gt] for gt in GT_ORDER})

class EntityCollection(object):
    def __init__(self):
        # populations
        self.population_idx_by_population_id = {}
        self.populationObjs = []
        # pairs
        self.pair_idx_by_pair_id = {}
        self.pairObjs = []
        # sequences
        self.sequence_idx_by_sequence_id = {}
        self.sequenceObjs = []
        # blocks
        self.block_idx_by_block_id = {}
        self.blockObjs = []
        # samples
        self.population_id_by_sample_id = {}
        # windows
        self.windowObjs = []
        
    def __str__(self):
        output = []
        output.append("# %s" % str(self.population_idx_by_population_id))
        for populationObj in self.populationObjs:
            output.append(str(populationObj))
        output.append("# %s" % str(self.pair_idx_by_pair_id))
        for pairObj in self.pairObjs:
            output.append(str(pairObj))
        output.append("# %s" % str(self.sequence_idx_by_sequence_id))
        for sequenceObj in self.sequenceObjs:
            output.append(str(sequenceObj))
        output.append("# %s" % str(self.block_idx_by_block_id))
        for blockObj in self.blockObjs:
            output.append(str(blockObj))
        return "%s\n" % "\n".join(output)

    def count(self, entity_type):
        if entity_type == 'populations':
            return len(self.populationObjs)
        elif entity_type == 'samples':
            return sum([len(populationObj.sample_ids) for populationObj in self.populationObjs])
        elif entity_type == 'pairs':
            return len(self.pairObjs)
        elif entity_type == 'sequences':
            return len(self.sequenceObjs)
        elif entity_type == 'bases':
            return sum([sequenceObj.length for sequenceObj in self.sequenceObjs])
        elif entity_type == 'blocks':
            return len(self.blockObjs)
        elif entity_type == 'windows':
            return len(self.windowObjs)
        else:
            return 0
    
    def sample_string_to_sample_ids(self, sample_string):
        return (sample_id for sample_id in sample_string.split(","))

    def sample_string_to_pair_idxs(self, sample_string):
        # works only for two pops ...
        pair_idxs = frozenset(filter(lambda x: x >= 0, [self.pair_idx_by_pair_id.get(frozenset(x), -1) for x in itertools.combinations(sample_string.split(","), 2)])) 
        if pair_idxs:
            return pair_idxs
        return np.nan

    def count_pair_idxs(self, pair_idxs):
        if not pair_idxs is np.nan:
            return len(list(pair_idxs))
        return 0

    def pair_idxs_to_sample_ids(self, pair_idxs):
        return set(itertools.chain.from_iterable([self.pairObjs[pair_idx].id for pair_idx in pair_idxs]))

    def add_population(self, population_idx, population_id):
        self.populationObjs.append(PopulationObj(population_idx, population_id))
        self.population_idx_by_population_id[population_id] = population_idx

    def add_pair(self, pair_idx, pair_id):
        self.pairObjs.append(PairObj(pair_idx, pair_id))
        self.pair_idx_by_pair_id[frozenset(pair_id)] = pair_idx        

    def add_sequence(self, sequence_idx, sequence_id, length):
        self.sequenceObjs.append(SequenceObj(sequence_idx, sequence_id, length))
        self.sequence_idx_by_sequence_id[sequence_id] = sequence_idx

    def add_blockObjs(self, blockObjs):
        blockObjs.sort(key=lambda x: (x.sequence_id, x.start)) 
        for idx, blockObj in enumerate(blockObjs):
            if blockObj.span < 64 or blockObj.span > 80:
                print(blockObj.id, blockObj.span)
            self.blockObjs.append(blockObj)
            self.block_idx_by_block_id[blockObj.id] = len(self.blockObjs) - 1
    
    def add_sample_to_populationObj(self, population_id, sample_id):
        population_idx = self.population_idx_by_population_id[population_id]
        self.populationObjs[population_idx].sample_ids.append(sample_id)
        self.population_id_by_sample_id[sample_id] = population_id

    def add_windowObjs(self, windowObjs):
        self.windowObjs = windowObjs

    def load_blocks(self, parameterObj, purpose='variants'):
        block_cols = ['block_id', 'length', 'span', 'sample_ids', 'pair_idxs', 'distance'] # maybe other colums need not be saved ...
        blocks_hdf5_store = pd.HDFStore(parameterObj.blocks_file)
        block_bed_df = pd.read_hdf(blocks_hdf5_store, key='bed')   
        block_df = pd.read_hdf(blocks_hdf5_store, 'block').reindex(columns=block_cols)
        blocks_hdf5_store.close()
        block_lengths = block_df.length.unique().tolist()
        if not len(block_lengths) == 1:
            exit('[X] Non uniform block length found : %s' % block_lengths)
        parameterObj.block_length = block_lengths[0]
        bed_tuples_by_block_id = collections.defaultdict(list)
        for block_id, sequence_id, start, end in block_bed_df.values.tolist():
            bed_tuples_by_block_id[block_id].append((sequence_id, start, end))
        blockObjs = []
        # parallelisation possible
        for block_id, length, span, sample_ids, _pair_idxs, distance in tqdm(block_df[block_cols].values.tolist(), total=len(block_df.index), desc="[%] ", ncols=100):
            if block_id in bed_tuples_by_block_id: # only Blocks in BED file are instanciated!
                pair_idxs = [int(pair_idx) for pair_idx in _pair_idxs.split(",")]
                blockObj = BlockObj(block_id, parameterObj.block_length)
                for bed_tuple in bed_tuples_by_block_id[block_id]:
                    bedObj = BedObj(bed_tuple[0], bed_tuple[1], bed_tuple[2], pair_idxs, bed_tuple[2] - bed_tuple[1])
                    blockObj.add_bedObj(bedObj, parameterObj, self)
                blockObj.mutype_counter_by_pair_idx = {pair_idx: collections.Counter() for pair_idx in pair_idxs}
                blockObjs.append(blockObj)
            else:
                print("[-] Block %s is missing from %s. Block will be ignored." % (block_id, parameterObj.block_bed))
        self.add_blockObjs(blockObjs)

        if purpose == 'windows':
            block_region_ids = (block_df["distance"].fillna(parameterObj.max_block_distance + 1).shift() > float(parameterObj.max_block_distance)).cumsum() # generate indices for splitting
            block_id_batches = []
            for idx, block_region_df in block_df.groupby(block_region_ids):
                #print("####", block_region_df)
                if len(block_region_df) > parameterObj.window_size: # remove regions below window_size
                    #block_region_df = block_region_df.drop(columns=['length', 'span', 'sample_ids', 'pair_idxs', 'distance'], axis=1) # remove distance/sample_idsx columns
                    block_id_batches.append(block_region_df['block_id'].tolist())
            if len(block_id_batches) == 0:
                exit("[X] Insufficient consecutive blocks with parameters '--window_size' (%s) and '--max_block_distance' (%s)" % (parameterObj.window_size, parameterObj.max_block_distance) )
            return block_id_batches

    def get_mutype_counters_by_block_id(self, parameterObj):
        variants_hdf5_store = pd.HDFStore(parameterObj.variants_file)
        variant_cols = ['block_id', 'pair_idx', 'hetA', 'fixed', 'hetB', 'hetAB']
        variants_block_df = pd.read_hdf(variants_hdf5_store, key='blocks').reindex(columns=variant_cols)
        variants_hdf5_store.close()
        mutype_counters_by_block_id = collections.defaultdict(list)
        for block_id, pair_idx, hetA, fixed, hetB, hetAB in tqdm(variants_block_df.values.tolist(), total=len(variants_block_df.index), desc="[%] ", ncols=100):
            mutype_counters_by_block_id[block_id].append(collections.Counter({'hetA': hetA, 'fixed': fixed, 'hetB': hetB, 'hetAB': hetAB}))
        return mutype_counters_by_block_id

    def transform_coordinates(self, parameterObj, coordinateTransformObj):
        params = [(blockObj, coordinateTransformObj) for blockObj in self.blockObjs]
        if parameterObj.threads < 2:
            with tqdm(total=len(params), desc="[%] ", ncols=100, unit_scale=True) as pbar:
                for param in params:
                    #print(blockObj, blockObj.void) 
                    self.transform_coordinates_blockObj(param)
                    #print(blockObj, blockObj.void)
                    pbar.update()
        else:
            with poolcontext(processes=parameterObj.threads) as pool:
                with tqdm(total=len(params), desc="[%] ", ncols=100, unit_scale=True) as pbar:
                    for blockObj in pool.imap_unordered(self.transform_coordinates_blockObj, params):
                        pbar.update()

    def transform_coordinates_blockObj(self, params):
        blockObj, coordinateTransformObj = params
        # print(">", blockObj.sequence_id, blockObj.start, blockObj.end)
        new_sequence_id, new_start, new_end = coordinateTransformObj.transform_coordinates(blockObj.sequence_id, blockObj.start, blockObj.end)
        # print("<", new_sequence_id, new_start, new_end)
        new_bed_tuples = []
        if new_sequence_id is None:
            blockObj.void = True
        else:
            blockObj.sequence_id = new_sequence_id
            blockObj.start = new_start
            blockObj.end = new_end
            for bed_tuple in blockObj.bed_tuples:
                new_bed_sequence_id, new_bed_start, new_bed_end = coordinateTransformObj.transform_coordinates(bed_tuple[0], bed_tuple[1], bed_tuple[2])
                new_bed_tuples.append((new_bed_sequence_id, new_bed_start, new_bed_end))
            new_bed_tuples.sort(key=lambda i: (i[0], i[1]))
            blockObj.bed_tuples = new_bed_tuples

    def parse_sample_file(self, parameterObj):
        if not parameterObj.sample_file:
            exit('[X] File cannot be read: %s' % parameterObj.args['sample_file'])
        # parse sample CSV
        samples_df = pd.read_csv(\
            parameterObj.sample_file, \
            sep=",", \
            usecols=[0, 1], \
            names=['sample_id', 'population_id'], \
            header=None, \
            dtype={ \
                'sample_id': 'category', \
                'population_id': 'category' \
                } \
            )
        if not len(samples_df.groupby('population_id').count().index) == 2:
             exit('[X] Invalid number of populations: %s (must be 2)' % len(samples_df.groupby('population_id').count().index))
        sorted_samples_df = samples_df.sort_values(['population_id', 'sample_id'])
        population_idx = -1
        sample_ids_by_population_id = collections.defaultdict(list)
        for sample_idx, (sample_id, population_id) in enumerate(sorted_samples_df.values.tolist()):
            if not population_id in self.population_idx_by_population_id:
                population_idx += 1
                self.add_population(population_idx, population_id)
            self.add_sample_to_populationObj(population_id, sample_id)
            sample_ids_by_population_id[population_idx].append(sample_id)
        for pair_idx, pair_id in enumerate([(x) for x in itertools.product(*sorted(sample_ids_by_population_id.values()))]):
            self.add_pair(pair_idx, pair_id)

    def parse_genome_file(self, parameterObj):
        df = pd.read_csv(parameterObj.genome_file, sep="\t", usecols=[0, 1], names=['sequence_id', 'length'], header=None)
        for sequence_idx, (sequence_id, length_str) in enumerate(df.values.tolist()):
            try:
                length = int(length_str)
            except TypeError:
                exit("[X] Line %s: Second column of --genome_file %s must be integers, not '%s'" % (sequence_idx, parameterObj.genome_file, length_str))
            self.add_sequence(sequence_idx, sequence_id, length)

    #def calculate_variation(self, mutype_counter, sites):
    #    # if missing, assume invariant
    #    pi_1 = float("%.8f" % ((mutype_counter['hetA'] + mutype_counter['hetAB']) / sites))
    #    pi_2 = float("%.8f" % ((mutype_counter['hetB'] + mutype_counter['hetAB']) / sites))
    #    d_xy = float("%.8f" % ((((mutype_counter['hetA'] + mutype_counter['hetB'] + mutype_counter['hetAB']) / 2.0) + mutype_counter['fixed']) / sites))
    #    mean_pi = (pi_1 + pi_2) / 2.0
    #    total_pi = (d_xy + mean_pi) / 2.0 # special case of pairwise Fst
    #    f_st = np.nan
    #    if (total_pi):
    #        f_st = float("%.8f" % ((total_pi - mean_pi) / total_pi)) # special case of pairwise Fst
    #    return pi_1, pi_2, d_xy, f_st

    def calculate_variation_from_df(self, df, sites):
        # print(df)
        pi_1 = float("%.8f" % ((df.hetA.sum() + df.hetAB.sum()) / sites))
        pi_2 = float("%.8f" % ((df.hetB.sum() + df.hetAB.sum()) / sites))
        d_xy = float("%.8f" % ((((df.hetA.sum() + df.hetB.sum() + df.hetAB.sum()) / 2.0) + df.fixed.sum()) / sites))
        mean_pi = (pi_1 + pi_2) / 2.0
        total_pi = (d_xy + mean_pi) / 2.0 # special case of pairwise Fst
        f_st = np.nan
        if (total_pi):
            f_st = float("%.8f" % ((total_pi - mean_pi) / total_pi)) # special case of pairwise Fst
        return (pi_1, pi_2, d_xy, f_st)

    def generate_modify_output(self, parameterObj):
        print("[#] Generating output...")
        self.generate_block_output(parameterObj, mode='modify')

    def generate_window_output(self, parameterObj, entityCollection, mutype_counters_by_block_id):
        mutuples = []
        #variants_df.set_index(keys=['block_id'], drop=False, inplace=True)
        mutuple_counters_by_window_id = collections.defaultdict(collections.Counter)
        window_vals = []
        for windowObj in tqdm(self.windowObjs, total=len(self.windowObjs), desc="[%] ", ncols=100, unit_scale=True):
            window_mutuples = []
            for block_id in windowObj.block_ids:
                mutype_counters = mutype_counters_by_block_id[block_id]
                for mutype_counter in mutype_counters:
                    mutype_tuple = mutype_counter_to_mutuple(mutype_counter, full=False)
                    # window_tally
                    mutuple_counters_by_window_id[windowObj.id][mutype_tuple] += 1
                    # global tally
                    mutuples.append(mutype_tuple)
                    # window list
                    window_mutuples.append(mutype_tuple)

            window_df = pd.DataFrame(window_mutuples, columns=MUTYPE_ORDER)
            #print(window_df)
            pi_1, pi_2, dxy, fst = self.calculate_variation_from_df(window_df, sites=(len(window_df.index) * parameterObj.block_length))
            mean_sample_fraction = (sum(windowObj.sample_counts) / len(self.population_id_by_sample_id)) / parameterObj.window_size
            window_vals.append([windowObj.id, windowObj.sequence_id, windowObj.start, windowObj.end, windowObj.span, windowObj.centre, pi_1, pi_2, dxy, fst, mean_sample_fraction])
        
        print("[#] Creating dataframe of window mutuple tallies...")
        window_mutuple_tally_cols = ['window_id', 'count'] + MUTYPE_ORDER 
        window_mutuple_tally_vals = []
        #print(mutuple_counters_by_window_id)
        for window_id, mutype_counters in mutuple_counters_by_window_id.items():
            #print(window_id, mutype_counters)
            for mutype, count in mutype_counters.most_common():
                window_mutuple_tally_vals.append([window_id, count] + list(mutype))
        window_mutuple_tally_df = pd.DataFrame(window_mutuple_tally_vals, columns=window_mutuple_tally_cols)

        print("[#] Creating dataframe of window metrics...")
        window_cols = [
            'window_id',
            'sequence_id',
            'start',
            'end', 
            'span',
            'centre',
            'pi_%s' % (self.populationObjs[0].id),
            'pi_%s' % (self.populationObjs[1].id),
            'dxy',
            'fst',
            'sample_cov'
            ]
        window_df = pd.DataFrame(window_vals, columns=window_cols)
        plot_fst_genome_scan(window_df, '%s.fst_genome_scan.png' % parameterObj.dataset, self.sequenceObjs)
        plot_pi_genome_scan(window_df, '%s.pi_genome_scan.png' % parameterObj.dataset, self.sequenceObjs)
        plot_pi_scatter(window_df, '%s.pi_scatter.png' % parameterObj.dataset)
        # storage
        window_hdf5_store = create_hdf5_store(
            out_f='%s.windows.h5' % (parameterObj.prefix), 
            path=parameterObj.path
            )
        window_df.to_hdf(window_hdf5_store, 'window_metrics', append=True)
        window_mutuple_tally_df.to_hdf(window_hdf5_store, 'mutypes', append=True)
        # tabulate_df(window_df, columns=window_cols, title="Windows")
        window_hdf5_store.close()

    def generate_variant_output(self, parameterObj):
        # collate data
        block_vals = []
        counts_by_site_by_sample_id_by_population_id = collections.defaultdict(dict)
        mutuples = []
        print("[#] Analysing variants in blocks...")
        for blockObj in tqdm(self.blockObjs, total=len(self.blockObjs), desc="[%] ", ncols=100, unit_scale=True):
            for sample_id in blockObj.sample_ids:
                population_id = self.population_id_by_sample_id[sample_id]
                gt_counter = blockObj.gt_counter_by_sample_id.get(sample_id, collections.Counter())
                if not sample_id in counts_by_site_by_sample_id_by_population_id[population_id]:
                    counts_by_site_by_sample_id_by_population_id[population_id][sample_id] = collections.Counter()
                counts_by_site_by_sample_id_by_population_id[population_id][sample_id] += gt_counter
            for pair_idx in blockObj.pair_idxs:
                mutype_counter = blockObj.mutype_counter_by_pair_idx.get(pair_idx, collections.Counter())
                mutuples.append(mutype_counter_to_mutuple(mutype_counter, full=False))
                block_vals.append([blockObj.id, pair_idx] + list(mutype_counter_to_mutuple(mutype_counter, full=True)))

        global_mutype_counter = collections.Counter(mutuples)
        print("[+] Monomorphic blocks: %s (%s of blocks)" % (global_mutype_counter[(0, 0, 0, 0)], format_percentage(global_mutype_counter[(0, 0, 0, 0)] / len(mutuples) )))
        print("[#] Creating dataframe of global mutuple tallies...")
        global_mutuple_tally_cols = ['count'] + MUTYPE_ORDER 
        global_mutuple_tally_vals = []
        for mutype, count in global_mutype_counter.most_common():
            global_mutuple_tally_vals.append([count] + list(mutype))
        global_mutuple_tally_df = pd.DataFrame(global_mutuple_tally_vals, columns=global_mutuple_tally_cols)
        # Mutype plot
        plot_mutuple_barchart(
            '%s.mutuple_barchart.png' % parameterObj.dataset,
            global_mutype_counter
            )
        print("[#] Creating dataframe of blocks...")
        block_idx = ['block_id', 'pair_idx']
        block_df = pd.DataFrame(block_vals, columns=block_idx + FULL_MUTYPE_ORDER)
        # add four gamete violations (FGV)
        block_df['FGV'] = np.where(( (block_df['fixed'] > 1) & (block_df['hetAB'] > 1) ), True, False) 

        print("[#] Creating dataframe for samples...")
        sample_df = pd.concat(
            {population_id: pd.DataFrame(dict_of_dicts).T for population_id, dict_of_dicts in counts_by_site_by_sample_id_by_population_id.items()}, 
            axis=0, 
            names=['population_id', 'sample_id'], 
            sort=False).fillna(0).astype(int)
        sample_df['sites'] = sample_df['blocks'] * parameterObj.block_length
        sample_df['heterozygosity'] = sample_df['HET'] / sample_df['sites']
        sample_df['missingness'] = sample_df['MISS'] / sample_df['sites']
        sample_df.drop(columns=['HET', 'HOM', 'MISS'], axis=1, inplace=True)
        tabulate_df(sample_df, columns = (sample_df.index.names + list(sample_df.columns)), title="Sample metrics: %s" % parameterObj.args.get('--prefix', 'gimble'))        

        print("[#] Creating dataframe for populations...")
        population_df = sample_df.groupby('population_id').mean()
        tabulate_df(population_df, columns=['population_id'] + list(population_df.columns), title="Population metrics: %s" % parameterObj.args.get('--prefix', 'gimble'))        

        print("[#] Creating dataframe for dataset...")
        #print(block_df)
        mutype_df = block_df.drop(columns=MUTYPE_OTHER + ['FGV'], axis=1)[block_idx + MUTYPE_ORDER].set_index(['block_id', 'pair_idx'])
        #print(mutype_df)
        pi_1, pi_2, dxy, fst = self.calculate_variation_from_df(mutype_df, sites=(parameterObj.block_length * len(mutype_df.index)))
        #print(parameterObj.block_length*len(mutype_df))
        #print((len(self.blockObjs) * parameterObj.block_length))
        sites = (len(self.blockObjs) * parameterObj.block_length)
        dataset_cols = [
            'blocks',
            'sites',
            'FGVs',
            'pi_%s' % (self.populationObjs[0].id),
            'pi_%s' % (self.populationObjs[1].id),
            'dxy',
            'fst'
            ]
        dataset_vals = [
            len(self.blockObjs),
            sites,
            (len(block_df[block_df['FGV'] == True]) / sites),
            pi_1,
            pi_2,
            dxy,
            fst
        ]
        dataset_df = pd.DataFrame.from_dict(dict(zip(dataset_cols, dataset_vals)), orient='index').T
        tabulate_df(dataset_df, columns=dataset_cols, title="Dataset metrics: %s" % parameterObj.args.get('--prefix', 'gimble'))        
        
        # storage
        variant_hdf5_store = create_hdf5_store(
            out_f='%s.variants.h5' % (parameterObj.prefix), 
            path=parameterObj.path, 
            )
        population_df.to_hdf(variant_hdf5_store, 'populations', append=True)
        sample_df.to_hdf(variant_hdf5_store, 'samples', append=True)
        block_df.to_hdf(variant_hdf5_store, 'blocks', append=True)
        global_mutuple_tally_df.to_hdf(variant_hdf5_store, 'mutypes', append=True)
        dataset_df.to_hdf(variant_hdf5_store, 'dataset', append=True)
        variant_hdf5_store.close()

        
        #for mutype in FULL_MUTYPE_ORDER:
            #print(mutype, np.unique(block_df[mutype], return_index=False, return_inverse=False, return_counts=True, axis=0))
        #non_missing_mutypes = block_df.loc[block_df['missing'] <= 4,:][MUTYPE_ORDER]
        #values, counts = np.unique(non_missing_mutypes[MUTYPE_ORDER].values, return_index=False, return_inverse=False, return_counts=True, axis=0)
        
    def generate_block_output(self, parameterObj, mode='blocks'):
        block_bed_cols = ['block_id', 'sequence_id', 'bed_start', 'bed_end']
        block_bed_vals = []
        block_cols = ['block_id', 'sequence_id', 'block_start', 'block_end', 'length', 'span', 'sample_ids', 'pair_idxs', 'count_samples', 'count_pairs']
        block_vals = []
        bases_blocked_by_sequence_id = collections.Counter()
        bases_blocked_by_sample_id = collections.Counter()
        bases_blocked_by_pair_id = collections.Counter()
        bases_blocked_by_sample_count = collections.Counter()
        void_count = 0
        for blockObj in tqdm(self.blockObjs, total=len(self.blockObjs), desc="[%] ", ncols=100, unit_scale=True):
            if blockObj.void:
                void_count += 1
            else:
                block_vals.append([
                    blockObj.id, 
                    blockObj.sequence_id,
                    blockObj.start,
                    blockObj.end,
                    blockObj.length, 
                    blockObj.span, 
                    ",".join([str(x) for x in sorted(blockObj.sample_ids)]), 
                    ",".join([str(x) for x in sorted(blockObj.pair_idxs)]),
                    len(blockObj.sample_ids), 
                    len(blockObj.pair_idxs)
                    ])
                for bed_tuple in blockObj.bed_tuples:
                    block_bed_vals.append([blockObj.id, bed_tuple[0],bed_tuple[1],bed_tuple[2]])
                bases_blocked_by_sequence_id[blockObj.sequence_id] += parameterObj.block_length
                # collect base counts
                for sample_count in range(len(blockObj.sample_ids), 0, -1):
                    bases_blocked_by_sample_count[sample_count] += parameterObj.block_length
                for sample_id in blockObj.sample_ids:
                    bases_blocked_by_sample_id[sample_id] += parameterObj.block_length
                for pair_idx in blockObj.pair_idxs:
                    bases_blocked_by_pair_id[self.pairObjs[pair_idx].id] += parameterObj.block_length
        block_bed_df = pd.DataFrame(block_bed_vals, columns=block_bed_cols)
        #tabulate_df(block_bed_df, columns=block_bed_cols, title="Block BED intervals")
        block_df = pd.DataFrame(block_vals, columns=block_cols)
        #tabulate_df(block_df, columns=block_cols, title="Blocks")
        out_f = '%s.blocks.h5' % parameterObj.prefix
        if mode == 'modify':
            print("[#] %s of %s blocks (%s) are being retained ..." % (
                (len(self.blockObjs) - void_count), 
                (len(self.blockObjs)), 
                (format_percentage((len(self.blockObjs) - void_count) / len(self.blockObjs)))))
            block_bed_df.sort_values(['sequence_id', 'bed_start'], ascending=[False, True], inplace=True)
            block_df.sort_values(['sequence_id', 'block_start'], ascending=[False, True], inplace=True)
            out_f = '%s.blocks.h5' % parameterObj.dataset
        block_hdf5_store = create_hdf5_store(
            out_f=out_f, 
            path=parameterObj.path
            )
        block_bed_df.to_hdf(block_hdf5_store, 'bed', append=True)
        block_df['distance'] = pd.to_numeric(np.where((block_df['sequence_id'] == block_df['sequence_id'].shift(-1)), block_df['block_start'].shift(-1) - block_df['block_end'], np.nan)) # compute distance to next interval)
        plot_distance_scatter(
            '%s.distance.png' % parameterObj.dataset,
            "Distance between blocks (in b)", 
            block_df['distance'].dropna().tolist()
            )
        barchart_y_vals, barchart_x_vals, barchart_labels, barchart_colours, barchart_populations = [], [], [], [], []
        for idx, (sample_id, population_id) in enumerate(self.population_id_by_sample_id.items()):
            barchart_populations.append(self.populationObjs[self.population_idx_by_population_id[population_id]].id)
            barchart_colours.append(self.populationObjs[self.population_idx_by_population_id[population_id]].colour)
            barchart_y_vals.append(bases_blocked_by_sample_id[sample_id])
            barchart_x_vals.append(idx)
            barchart_labels.append(sample_id)
        plot_sample_barchart(
            '%s.blocks_per_sample.png' % parameterObj.dataset,
            "Bases in blocks by sample", 
            barchart_y_vals, barchart_x_vals, barchart_labels, barchart_colours, barchart_populations
            )
#        #plot_shared_blocks(
        #   '%s.shared_blocks_between_samples.png' % parameterObj.dataset,
        #   "Distance between blocks (in b)", 
        #   block_df['distance'].dropna().tolist()
         #  )
        block_df.to_hdf(block_hdf5_store, 'block', append=True)
        block_hdf5_store.close()
        # call plot from here
        #print(bases_blocked_by_sequence_id)
        #print(bases_blocked_by_sample_id)
        #print(bases_blocked_by_pair_id)
        #print(bases_blocked_by_sample_count)
        #print(bases_blocked_by_pair_count)

class BedObj(object):
    __slots__ = ['sequence_id', 'start', 'end', 'pair_idxs', 'length']
    
    def __init__(self, sequence_id, start, end, pair_idxs, length):
        self.sequence_id = sequence_id
        self.start = int(start)
        self.end = int(end)
        self.length = int(length)
        self.pair_idxs = set(pair_idxs)

    def __str__(self):
        return "\t".join([self.sequence_id, str(self.start), str(self.end), str(self.length), str(self.pair_idxs)]) 

class WindowObj(object):
    #__slots__ = ["sequence_id", 
    #             "start", 
    #             "end", 
    #             "id", 
    #             "span", 
    #             "centre", 
    #             "block_ids", 
    #             "sample_counts"
    #             ]

    def __init__(self, blockObjs):
        self.sequence_id = blockObjs[0].sequence_id
        self.start = blockObjs[0].start
        self.end = blockObjs[-1].end
        self.id = "%s_%s_%s" % (self.sequence_id, self.start, self.end)
        self.span = self.end - self.start
        self.centre = self.start + (self.span / 2)
        self.block_ids = [blockObj.id for blockObj in blockObjs]
        #self.variant_weights = [len(blockObj.pair_idxs) for blockObj in blockObjs]
        self.sample_counts = [len(blockObj.sample_ids) for blockObj in blockObjs]

class CoordinateTransformObj(object):
    def __init__(self):
        self.tuples_by_sequence_id = collections.defaultdict(list)

    def add_tuple(self, sequence_id, sequence_start, sequence_end, sequence_orientation, chrom_id, chrom_start, chrom_end):
        self.tuples_by_sequence_id[sequence_id].append((sequence_id, sequence_start, sequence_end, sequence_orientation, chrom_id, chrom_start, chrom_end))

    def transform_coordinates(self, old_id, old_start, old_end):
        new_id, new_start, new_end = None, None, None
        if old_id in self.tuples_by_sequence_id:
            for interval in self.tuples_by_sequence_id[old_id]:
                sequence_id, sequence_start, sequence_end, sequence_orientation, chrom_id, chrom_start, chrom_end = interval
                if old_start >= sequence_start and old_end <= sequence_end:
                    new_id = chrom_id 
                    if sequence_orientation == "+":
                        new_start = chrom_start + (old_start - sequence_start)
                        new_end = chrom_start + (old_end - sequence_start)
                        return (new_id, new_start, new_end)
                    else:
                        new_start = chrom_start + (sequence_end - old_end)
                        new_end = chrom_start + (sequence_end - old_start)
                        return (new_id, new_start, new_end)
        return (new_id, new_start, new_end)

class BlockObj(object):

    __slots__ = [
        "sequence_id", 
        "id", 
        "pair_idxs", 
        "sample_ids", 
        "start", 
        "end", 
        "length", 
        "span", 
        "score", 
        "needed", 
        "bed_tuples", 
        "mutype_counter_by_pair_idx", 
        "gt_counter_by_sample_id", 
        "void" 
        ]

    def __init__(self, block_id, block_length):
        self.sequence_id = block_id.split(".")[0]
        self.id = block_id
        self.pair_idxs = None
        self.sample_ids = None
        self.start = None
        self.end = None
        self.length = 0
        self.span = 0 
        self.score = 0.0
        self.needed = block_length
        self.bed_tuples = [] # list of tuples (sequence_id, start, end) of consecutive regions

        self.mutype_counter_by_pair_idx = {} # dict of counters 
        self.gt_counter_by_sample_id = {} # dict of counters
        self.void = False

    def __str__(self):
        return "[B] ID=%s %s %s %s LEN=%s SPAN=%s SCORE=%s [P]=%s\n%s\n%s" % (
            self.id,
            self.sequence_id, 
            self.start, 
            self.end, 
            format_bases(self.length), 
            format_bases(self.span), 
            format_fraction(self.score, 1), 
            self.pair_idxs, 
            str(self.mutype_counter_by_pair_idx),
            str(self.gt_counter_by_sample_id)
            )

    def __nonzero__(self):
        if self.length:
            return True
        return False

    def add_bedObj(self, bedObj, parameterObj, entityCollection):
        '''
        Function for adding a bedObj to the blockObj
        [parameters]
            - bedObj to be added
            - parameterObj
        [returns]
            a) None (if bedObj has been consumed)
            b) bedObj
                b1) original bedObj (if span-violation)
                b2) remainder bedObj (if bedObj.length > blockObj.needed)
        [comments]
            - span-violation:
                if blockObj.span > parameterObj.max_block_length:
                    - original bedObj is returned, blockObj.score is set to 0.0
            - blockObj.needed: allows distinction between 
                a) finished block: blockObj.needed = 0
                b) virgin block: blockObj.needed = parameterObj.block_length
                c) started block: 0 < blockObj.needed < parameterObj.block_length
            - blockObj.score: 
                a) if blockObj.needed == parameterObj.block_length (virgin block):
                    blockObj.score = (bedObj.pair_idxs / parameterObj.pairs_count) * (min(bedObj.length, required_length) / parameterObj.block_length)
                b) if blockObj.needed < parameterObj.block_length (non-virgin block):
                    blockObj.score = (len(blockObj.pair_idxs.intersection(bedObj.pair_idxs)) / parameterObj.pairs_count) * (min(bedObj.length, required_length) / parameterObj.block_length)
                c) if span-violation:
                    blockObj.score = 0.0
        '''
        interval_length = min(self.needed, bedObj.length)
        block_end = bedObj.start + interval_length
        try:
            _span = block_end - self.start # TypeError: int() argument must be a string, a bytes-like object or a number, not 'NoneType' 
        except TypeError:
            self.start = bedObj.start
            _span = block_end - self.start
        if parameterObj.max_block_length and _span > parameterObj.max_block_length:
            self.score = 0.0
            return bedObj
        try:
            self.pair_idxs = self.pair_idxs.intersection(bedObj.pair_idxs) # AttributeError: 'NoneType' object has no attribute 'intersection'
        except AttributeError:
            self.pair_idxs = bedObj.pair_idxs
        self.end = block_end 
        self.span = _span
        self.length += interval_length
        self.needed -= interval_length
        self.score = (len(self.pair_idxs) / len(entityCollection.pairObjs)) * (self.length / parameterObj.block_length)
        self.sample_ids = entityCollection.pair_idxs_to_sample_ids(self.pair_idxs)   
        try:
            last_tuple = self.bed_tuples.pop()
            if (last_tuple[2] - bedObj.start) == 0: # no gap
                self.bed_tuples.append((bedObj.sequence_id, last_tuple[1], block_end)) 
            else: # gap
                self.bed_tuples.append(last_tuple)
                self.bed_tuples.append((bedObj.sequence_id, bedObj.start, block_end)) 
        except IndexError:
            self.bed_tuples.append((bedObj.sequence_id, bedObj.start, block_end))
        self.sequence_id = bedObj.sequence_id
        if interval_length < bedObj.length:
            return BedObj(bedObj.sequence_id, (bedObj.start + interval_length), bedObj.end, bedObj.pair_idxs, (bedObj.length - interval_length))    
        else:
            return None

class PopulationObj(object):
    def __init__(self, idx, population_id):
        self.idx = idx
        self.id = population_id
        self.sample_ids = []
        self.colour = COLOURS[idx]
        
    def __str__(self):
        return "[Population] %s %s %s" % (self.idx, self.id, ", ".join(self.sample_ids))

class PairObj(object):
    def __init__(self, idx, pair_id):
        self.idx = idx
        self.id = pair_id # tuple, use this for sample_ids since order is preserved!

    def __str__(self):
        return "[Pair] %s %s" % (self.idx, self.id)

class SequenceObj(object):
    def __init__(self, sequence_idx, sequence_id, length):
        self.idx = sequence_idx
        self.id = sequence_id
        self.length = length

    def __str__(self):
        return "[Sequence] %s %s %s" % (self.idx, self.id, format_bases(self.length))