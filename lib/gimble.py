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
import configparser
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
from matplotlib.path import Path
import matplotlib.patches as patches
import numpy as np
import cerberus
import lib.simulate
import hashlib 
from timeit import default_timer as timer
import fractions

# np.set_printoptions(threshold=sys.maxsize)


'''
[Rules for better living]

- gimbleStore.data.attrs (meta): ZARR JSON encoder does not like numpy/pandas dtypes, have to be converted to python dtypes

- Coordinate systems:
    {GIMBLE}                    : 0-based       | {GIMBLE}  
    -----------------------------------------------------------
    {GENOME_FILE} : LENGTH      : 0-based   =>  | [0, l)
    {BED_FILE}    : START:END   : 0-based   =>  | [START, END)
    {VCF_FILE}    : POS         : 1-based   =>  | [POS-1)
'''

'''
[To Do]
- generalise prefligth-checks for relevant modules

- inference
    - calculate different Pi's

-info:
    - window
    
- QC plots
    - variants 
        - plot barcharts of HOMREF/HOMALT/HET/MISS/MULTI as proportion of total records
    - intervals
        - plot length
        - plot distance
    - mutuples
        - mutuple pcp

- data dumping 
- metrics dumping

'''

PURPLE = '#4F3D63'

DFRM = '[──'
DFRT = '├──'
DPPM = '    '
DFRL = '└──'
DPPL = '│   '
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

SIGNS = {'T': '├──', 'F': '└──', 'S': '    ', 'P': '│   ', 'B': '───'}

class ReportObj(object):
    '''Info class for making reports
    '''
    def __init__(self, width=80):
        '''
        Parameters 
        ----------
        width : int, width of report in characters, default=80
        
        Returns
        -------
        out : instance
        '''
        self.width = width
        self.out = []
    
    def add_line(self, prefix="[+]", branch="", left='', center='', right='', fill=' '):
        '''
        Writes line to Info Object
        Lines is composed of {prefix} {branch} {left} {fill} {center} {fill} {right}
        SIGNS = {'T': '├──', 'F': '└──', 'S': '    ', 'P': '│   ', 'B': '───'}
        '''
        start = "%s %s" % (prefix, "".join([SIGNS[b] for b in branch])) if branch else prefix
        w = self.width - len(left) - len(right) - len(start) 
        self.out.append("%s %s %s %s" % (start, left, (" %s " % center).center(w, fill) if center else fill * w, right))

    def __add__(self, other):
        if isinstance(other, ReportObj):
            if self.width == other.width:
                result = ReportObj(width=self.width)
                result.out = self.out + other.out
                return result
            else:
                raise ValueError("addition only supported for ReportObj of equal width")    
        else:
            raise ValueError("addition only supported for ReportObj's")

    def __radd__(self, other):
        return self.__add__(other)


    def __repr__(self):
        return "\n".join(self.out)

def recursive_get_size(path: str) -> int:
    """Gets size in bytes of the given path, recursing into directories."""
    if os.path.isfile(path):
        return os.path.getsize(path)
    if not os.path.isdir(path):
        return 0
    return sum(
        recursive_get_size(os.path.join(path, name))
        for name in os.listdir(path))

def parse_csv(csv_f='', dtype=[], usecols=[], sep=',', header=None):
    '''dtypes := "object", "int64", "float64", "bool", "datetime64", "timedelta", "category"'''
    df = pd.read_csv(csv_f, sep=sep, usecols=usecols, names=list(dtype.keys()), header=header, dtype=dtype)
    if df.isnull().values.any():
        sys.exit("[X] Bad file format %r." % csv_f)
    return df

def format_bases(bases: int) -> str:
    return "%s b" % format(bases, ',d')

def format_percentage(fraction: float, precision=2) -> str:
    return "{:.{}%}".format(fraction, precision)

def format_proportion(fraction, precision=2):
    if fraction == '-':
        return '-'
    return "{:.{}f}".format(fraction, precision)

def format_count(count):
    return "%s" % str(format(count, ',d'))

def format_bytes(size, precision=1):
    power = 2**10
    n = 0
    power_labels = {0 : ' B', 1: 'KB', 2: 'MB', 3: 'GB', 4: 'TB'}
    while size > power:
        size /= power
        n += 1
    return "{:.{}f} {}".format(size, precision, power_labels[n])

def get_n50_from_lengths(lengths):
    length_sorted = sorted(lengths, reverse=True)
    cum_sum = np.cumsum(length_sorted)
    half = int(sum(lengths) / 2)
    cum_sum_2 =min(cum_sum[cum_sum >= half])
    n50_idx = np.where(cum_sum == cum_sum_2)
    return length_sorted[int(n50_idx[0][0])]
    
def check_unique_pos(pos_array):
    unique_pos, counts_pos = np.unique(pos_array, return_counts=True)
    duplicates = unique_pos[counts_pos > 1]
    if duplicates.any():
        print("\n[-] %s VCF records with non-unique positions found. Rescuing records by shifting position... (abort if this is not desired)" % (len(duplicates)))
        pos_array = fix_pos_array(pos_array)
    return pos_array

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

def szudzik_pairing(folded_minor_allele_counts):
    # adapted from: https://drhagen.com/blog/superior-pairing-function/
    return np.where(
            (folded_minor_allele_counts[:,0] >= folded_minor_allele_counts[:,1]),
            np.square(folded_minor_allele_counts[:,0]) + folded_minor_allele_counts[:,0] + folded_minor_allele_counts[:,1],
            folded_minor_allele_counts[:,0] + np.square(folded_minor_allele_counts[:,1])
            )
    #if isinstance(folded_minor_allele_counts, np.ndarray):
    #    # assumes folded_minor_allele_counts is array with shape (n,2)
    #    return np.where(
    #        (folded_minor_allele_counts[:,0] >= folded_minor_allele_counts[:,1]),
    #        np.square(folded_minor_allele_counts[:,0]) + folded_minor_allele_counts[:,0] + folded_minor_allele_counts[:,1],
    #        folded_minor_allele_counts[:,0] + np.square(folded_minor_allele_counts[:,1])
    #        )
    #elif isinstance(folded_minor_allele_counts, tuple):
    #    # assumes folded_minor_allele_counts is tuple of size 2
    #    a, b = folded_minor_allele_counts
    #    if a >= b:
    #        return (a**2) + a + b 
    #    return a + (b**2)
    #else:
    #    pass

def get_coverage_counts(coverages, idxs, num_blocks):
    # not used at the moment
    num_sample_sets = idxs[-1] + 1 # this is correct, don't worry about it ...
    temp = coverages + (num_sample_sets * np.arange(coverages.shape[0]))[:, None]
    blocks_per_sample_set = np.bincount(temp.ravel(), minlength=(num_sample_sets * coverages.shape[0])).reshape(-1, num_sample_sets)
    # remove columns/sample_sets that only contain zeroes and normalise
    return blocks_per_sample_set[:, ~(blocks_per_sample_set == 0).all(0)] / (num_blocks / len(idxs)) / num_blocks 

def block_sites_to_variation_arrays(block_sites, cols=np.array([1,2,3]), max_type_count=7):
    temp_sites = block_sites + (max_type_count * np.arange(block_sites.shape[0], dtype=np.int64).reshape(block_sites.shape[0], 1))
    # return multiallelic, missing, monomorphic, variation
    return np.hsplit(np.bincount(temp_sites.ravel(), minlength=(block_sites.shape[0] * max_type_count)).reshape(-1, max_type_count), cols)

# def calculate_popgen_from_array(mutype_array, sites):
#     # print('# Mutypes: 0=MULTI, 1=MISS, 2=MONO, 3=HetB, 4=HetA, 5=HetAB, 6=Fixed')
#     pi_1 = float("%.8f" % np.divide(np.sum(mutype_array[:, 1]) + np.sum(mutype_array[:, 2]), sites)) # average heterozygosity
#     pi_2 = float("%.8f" % np.divide(np.sum(mutype_array[:, 0]) + np.sum(mutype_array[:, 2]), sites)) # average heterozygosity
#     d_xy = float("%.8f" % np.divide(np.divide(np.sum(mutype_array[:, 0]) + np.sum(mutype_array[:, 1]) + np.sum(mutype_array[:, 2]), 2.0) + np.sum(mutype_array[:, 3]), sites))
#     mean_pi = (pi_1 + pi_2) / 2.0
#     total_pi = (d_xy + mean_pi) / 2.0 # special case of pairwise Fst
#     f_st = np.nan
#     if (total_pi):
#         f_st = float("%.8f" % ((total_pi - mean_pi) / total_pi)) # special case of pairwise Fst
#     fgv = len(mutype_array[(mutype_array[:, 2] > 0) & (mutype_array[:, 3] > 0)])
#     return (pi_1, pi_2, d_xy, f_st, fgv)

def genotype_to_mutype_array(sa_genotype_array, idx_block_sites_in_pos, block_sites, debug=True):
    '''
    - possible errors:
        - if whole sequence has only missing GTs in genotypes for a sample set, then np_allele_count_array will be empty ... (should never happen)
    '''
    np_genotype_array = np.array(sa_genotype_array)
    #print('np_genotype_array', np_genotype_array.shape, np_genotype_array)
    np_allele_count_array = np.ma.masked_equal(sa_genotype_array.count_alleles(), 0, copy=True) 
    #print('np_allele_count_array', np_allele_count_array.shape, np_allele_count_array)
    #print('np.any(np_allele_count_array)', np.any(np_allele_count_array))
    allele_map = np.ones((np_allele_count_array.shape), dtype=np.int64) * np.arange(np_allele_count_array.shape[-1], dtype=np.int64)
    #print('allele_map', allele_map.shape, allele_map)
    #print('np.any(allele_map)', np.any(allele_map))
    if np.any(np_allele_count_array) and np.any(allele_map):
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
    #print(folded_minor_allele_counts)
    folded_minor_allele_counts[np.any(sa_genotype_array.is_missing(), axis=1)] = np.ones(2) * -1        # -1, -1 for missing => -1
    folded_minor_allele_counts[(np_allele_count_array.count(axis=1) > 2)] = np.ones(2) * (-1, -2)       # -1, -2 for multiallelic => -2

    block_sites[idx_block_sites_in_pos] = szudzik_pairing(folded_minor_allele_counts) + 2               # add 2 so that not negative for bincount
    block_sites[~idx_block_sites_in_pos] = 2                                                            # monomorphic = 2 (0 = multiallelic, 1 = missing)
    if debug == True:
        print("# block_sites as mutype array", block_sites)
    if debug == True:
        block_sites_pos = block_sites.flatten()
        pos_df = pd.DataFrame(block_sites_pos[idx_block_sites_in_pos.flatten()], dtype='int64', columns=['pos'])
        genotypes_df = pd.DataFrame(np_genotype_array.reshape(np_genotype_array.shape[0], 4), dtype='i4', columns=['a1', 'a2', 'b1', 'b2'])        
        block_sites_df = pos_df.join(genotypes_df)
        folded_minor_allele_count_df = pd.DataFrame(folded_minor_allele_counts, dtype='int8', columns=['fmAC_a', 'fmAC_b'])
        block_sites_df = block_sites_df.join(folded_minor_allele_count_df)
        variants = pd.DataFrame(block_sites[idx_block_sites_in_pos], dtype='int', columns=['SVar'])
        block_sites_df = block_sites_df.join(variants)
        print('# Mutypes: 0=MULTI, 1=MISS, 2=MONO, 3=HetB, 4=HetA, 5=HetAB, 6=Fixed')
        print(block_sites_df)
    return block_sites

def _harmonic(a, b):
    if b-a == 1:
        return fractions.Fraction(1,a)
    m = (a+b)//2
    return _harmonic(a,m) + _harmonic(m,b)

def harmonic(n):
    '''https://fredrik-j.blogspot.com/2009/02/how-not-to-compute-harmonic-numbers.html'''
    return _harmonic(1,n+1)

def create_ranges(aranges):
    # does the transformation form 0-based (BED) to 1-based (VCF) coordinate system 
    # https://stackoverflow.com/a/47126435
    l = (aranges[:, 1] - aranges[:, 0])
    clens = l.cumsum()
    if np.any(clens):
        ids = np.ones(clens[-1], dtype=np.int64)
        ids[0] = aranges[0, 0]
        ids[clens[:-1]] = aranges[1:, 0] - aranges[:-1, 1] + 1
        return ids.cumsum()
    return None

def cut_windows(mutype_array, idxs, start_array, end_array, num_blocks=10, num_steps=3):
    coordinate_sorted_idx = np.argsort(end_array)
    mutype_array_sorted = mutype_array.take(coordinate_sorted_idx, axis=0)
    window_idxs = np.arange(mutype_array_sorted.shape[0] - num_blocks + 1)[::num_steps, None] + np.arange(num_blocks)
    window_mutypes = mutype_array_sorted.take(window_idxs, axis=0)
    block_starts = start_array.take(coordinate_sorted_idx, axis=0).take(window_idxs, axis=0)
    window_starts = np.min(start_array.take(coordinate_sorted_idx, axis=0).take(window_idxs, axis=0), axis=1).T
    block_ends = end_array.take(coordinate_sorted_idx, axis=0).take(window_idxs, axis=0)#[:,-1]
    window_ends = np.max(end_array.take(coordinate_sorted_idx, axis=0).take(window_idxs, axis=0), axis=1).T
    window_midpoints = (block_starts / 2) + (block_ends / 2)
    window_pos_mean = np.mean(window_midpoints, axis=1).T
    window_pos_median = np.median(window_midpoints, axis=1).T
    return window_mutypes, window_starts, window_ends, window_pos_mean, window_pos_median

def cut_blocks(interval_starts, interval_ends, block_length, block_span, block_gap_run):
    sites = create_ranges(np.array((interval_starts, interval_ends), dtype=np.int64).T)
    if sites is None: 
        return None
    block_sites = np.concatenate([
        x[:block_length * (x.shape[0] // block_length)].reshape(-1, block_length) 
            for x in np.split(sites, np.where(np.diff(sites) > block_gap_run)[0] + 1)
        ]) 
    block_span_valid_mask = (((block_sites[:, -1] - block_sites[:, 0]) + 1) <= block_span)
    return block_sites[block_span_valid_mask]


def bsfs_to_2d(bsfs):
    """Returns 2D array with counts, mutuples. First column is counts.

    Parameters 
    ----------
    bsfs : ndarray, int, ndim (4)  
    
    Returns
    -------
    out : ndarray, int, ndim (2)

    """
    non_zero_idxs = np.nonzero(bsfs)
    return np.concatenate([bsfs[non_zero_idxs].reshape(non_zero_idxs[0].shape[0], 1), np.array(non_zero_idxs).T], axis=1)

def bsfs_to_counter(bsfs):
    """Returns (dict of) collections.Counter based on bsfs_array of 4 or 5 dimensions.

    Parameters 
    ----------
    bsfs : ndarray, int, ndim (4) or (5)  
    
    Returns
    -------
    out : (dict of) collections.Counter

    """
    non_zero_idxs = np.nonzero(bsfs)
    if bsfs.ndim == 5:
        out = collections.defaultdict(collections.Counter())
        for idx, count in zip(tuple(np.array(non_zero_idxs).T), bsfs[non_zero_idxs]):
            window_idx, mutuple = idx[0], idx[1:]
            out[window_idx][mutuple] = count
    elif bsfs.ndim == 4:
        out = collections.defaultdict(lambda: collections.Counter())
        for mutuple, count in zip(tuple(np.array(non_zero_idxs).T), bsfs[non_zero_idxs]):
            out[window_idx][mutuple] = count
    else:
        raise ValueError('bsfs must have 4 or 5 dimensions')
    return out

def ordered_intersect(a=[], b=[], order='a'):
    A, B = a, b
    if order == 'b' :
        B, A = a, b
    return [_a for _a in A if _a in set(B)]

class CustomNormalizer(cerberus.Validator):
    def __init__(self, *args, **kwargs):
        super(CustomNormalizer, self).__init__(*args, **kwargs)
        self.valid_pop_ids = kwargs['valid_pop_ids']
        self.valid_sync_pops = kwargs['valid_sync_pops']
        
    def _normalize_coerce_float_or_list(self, value):
        try:
            return float(value)
        except:
            values = value.strip('()[]').split(",")
            if len(values) == 3:
                try:
                    return [float(v) for v in values]
                except ValueError:
                    return None
            elif len(values) == 2:
                try:
                    values = [float(v) for v in values]
                    mid = np.mean(values)
                    return [mid,]+values
                except ValueError:
                    return None
                values = []
            elif len(values) == 5:
                valid_scales = ['lin', 'log']
                if not values[4].strip(' ') in valid_scales:
                    return None
                try:
                    return [float(v) for v in values[:3]]+[int(values[3]),values[4].strip(' ')]
                except ValueError:
                    return None
            else:
                return None

    def _normalize_coerce_float_or_empty(self, value):
        try:
            return float(value)
        except:
            if value.strip(' ') =='':
                return ''
            else:
                return None

    def _normalize_coerce_int_or_empty(self, value):
        try:
            return int(value)
        except:
            if value.strip(' ')=='':
                return ''
            else:
                return None

    def _validate_notNoneInt(self, notNoneNumeric, field, value):
        """
        {'type':'boolean'}
        """
        if value == None and notNoneNumeric:
            self._error(field, "Must be an int value or empty")
    
    def _validate_notNoneFloat(self, notNoneNumeric, field, value):
        """
        {'type':'boolean'}
        """
        if value == None and notNoneNumeric:
            self._error(field, "Must be a float or empty")

    def _validate_notNone(self, notNone, field, value):
        """
        {'type':'boolean'}
        """
        if not value and notNone:
            self._error(field, "Must be a float, or a list of length 3 or 5.")

    def _validate_isPop(self, isPop, field, value):
        """
        {'type':'boolean'}
        """
        if value.strip(" ") not in self.valid_pop_ids:
            self._error(field, "Must be either A, B or A_B")

    def _validate_isPopSync(self, isPopSync, field, value):
        """
        {'type':'boolean'}
        """
        if value.strip(" ") !='':
            if value.strip(" ") not in self.valid_sync_pops:
                self._error(field, "Must be either A,A_B, A,B or A_B,B")

    def _validate_isPath(self, isPath, field, value):
        """
        {'type':'boolean'}
        """
        if value.strip(" ") != '' and not os.path.isfile(value):
            self._error(field, 'Must be a valid path to the recombination map.')

class ParameterObj(object):
    '''Superclass ParameterObj'''
    def __init__(self, params):
        self._PATH = params['path']
        self._VERSION = params['version']
        self._MODULE = params['module']
        self._CWD = params['cwd']
        self.config = None

    def __repr__(self):
        return("[+] VER := %s\n[+] CWD := %s\n[+] CMD := %s\n" % (
            self._VERSION, 
            self._CWD, 
            self._get_cmd()
            ))

    def _cast_to_repeated_list(self, x, repeat=1):
        if isinstance(x, list):
            return x
        elif isinstance(x, str):
            return [x]*repeat
        try:
            return list(x)*repeat
        except TypeError:
            return [x]*repeat

    def _dict_product(self):
        if len(self.config["parameters"])>0:
            return [dict(zip(self.config["parameters"], x)) for x in itertools.product(*self.config["parameters"].values())]

    def _dict_zip(self, pdict):
        return [dict(zip(pdict, x)) for x in zip(*pdict.values())]

    def _expand_params(self):
        if len(self.config['parameters'])>0:
            for key, value in self.config['parameters'].items():
                if isinstance(value, float) or isinstance(value, int):
                    self.config['parameters'][key] = [value,]
                elif key=='recombination':
                    pass
                else:
                    if len(value) == 5:
                        midv, minv, maxv, n, scale = value
                        if scale.startswith('lin'):
                            sim_range = self._expand_params_lin(midv, minv, maxv, n)
                        elif scale.startswith('log'):
                            sim_range = self._expand_params_log(midv, minv, maxv, n)
                        else:
                            raise ValueError
                        self.config['parameters'][key] = sim_range
                    elif len(value) <4:
                        self.config['parameters'][key] = np.unique(value)
                    else:
                        sys.exit("[X] Uncaught error in config file configuration.")

    def _expand_params_lin(self, pcentre, pmin, pmax, psamples):
        starts = [pmin, pcentre]
        ends = [pcentre, pmax]
        nums = [round(psamples/2) + 1, round(psamples/2)] if psamples % 2 == 0 else [round(psamples/2) + 1, round(psamples/2) + 1]
        return np.unique(np.concatenate([np.linspace(start, stop, num=num, endpoint=True, dtype=np.float64) for start, stop, num in zip(starts, ends, nums)]))    

    def _expand_params_log(self, pcentre, pmin, pmax, psamples):
        plogcentre, plogmin, plogmax = np.log10([pcentre, pmin, pmax])
        temp =  np.logspace(plogmin, plogmax, num=psamples, endpoint=True, dtype=np.float64)
        left = (temp<pcentre).sum()
        right = psamples-left
        left = np.logspace(plogmin, plogcentre, num=left, endpoint=False, dtype=np.float64)
        right = np.logspace(plogcentre, plogmax, num=right, endpoint=True, dtype=np.float64)
        return np.unique(np.concatenate([left, right]))

    def _get_blocks_length(self, zstore=None):
        blocks_length_zarr = None
        blocks_length_ini = self._get_int(self.config['mu']['blocklength'], ret_none=True)
        if zstore:
            z = zarr.open(self.zstore, mode='r')
            blocks_length_zarr =  z['seqs'].attrs.get('blocks_length')
        if blocks_length_zarr and blocks_length_ini:
            if blocks_length_zarr != blocks_length_ini:
                print("[-] block length in ini file and used in GStore differ. Using blocks_length_ini.")
                return blocks_length_ini
        if blocks_length_ini:
            return blocks_length_ini
        elif blocks_length_zarr:
            return blocks_length_zarr
        else:
            sys.exit("[X] Blocklength needs to be specified in ini file.")

    def _get_int(self, string, ret_none=False):
        try:
            return int(string)
        except TypeError:
            if not ret_none:
                sys.exit("[X] %r can't be converted to interger." % string)
            return None

    def _get_float(self, string, ret_none=False):
        try:
            return float(string)
        except TypeError:
            if not ret_none:
                sys.exit("[X] %r can't be converted to float." % string)
            return None

    def _get_cmd(self):
        return "%s/gimble %s %s" % (
            self._PATH, 
            self._MODULE, 
            "".join(["--%s " % " ".join((k, str(v))) for k,v in self.__dict__.items() if not k.startswith("_")]))

    def _get_path(self, infile, path=False):
        if infile is None:
            return None
        _path = pathlib.Path(infile).resolve()
        if not _path.exists():
            sys.exit("[X] File not found: %r" % str(infile))
        if path:
            return _path
        return str(_path)

    def _get_prefix(self, prefix):
        if prefix is None:
            return None
        path = pathlib.Path(prefix).resolve()
        new_path = pathlib.Path(prefix+'.z').resolve()
        if new_path.exists():
            sys.exit("[X] zstore already exists. Specify using -z or provide new prefix.")
        parent_path = pathlib.Path(prefix).resolve().parent
        if not parent_path.exists():
            sys.exit("[X] File not found: %r" % str(infile))
        return str(path)

    def _get_pops_to_sync(self):
        reference, to_be_synced = None, None
        syncing =  self.config['populations']['sync_pop_sizes']
        if syncing:
            if len(syncing)>0:
                syncing = syncing.split(',')
                reference = syncing[0]
                to_be_synced = syncing[1:]
                if any(isinstance(self.config['parameters'][f'Ne_{pop}'], list) or isinstance(self.config['parameters'][f'Ne_{pop}'], float) for pop in to_be_synced):
                    print(f"[-] Ne_{', Ne_'.join(to_be_synced)} is specified in config file but synced with Ne_{reference}.")
        return (reference, to_be_synced)

    def _get_threads(self, num):
        try:
            num=int(num)
        except TypeError:
            return (1,1)
        if num == 1:
            return (1,1)
        else: 
            split = [num%i for i in (2,3)]
            threads = min(split)
            gridThreads = num//threads
        return (threads, gridThreads)

    def _remove_pop_from_dict(self, toRemove):
        if toRemove:
            for pop in toRemove:
                del self.config['parameters'][f'Ne_{pop}']

    def _sync_pop_sizes(self, reference, toBeSynced):
        if toBeSynced and reference:
            for pop in toBeSynced:
                for paramCombo in self.parameter_combinations:
                    paramCombo[f'Ne_{pop}'] = paramCombo[f'Ne_{reference}']

    def _sync_pop_sizes_optimise(self, reference, toBeSynced):
        if toBeSynced and reference:
            for pop in toBeSynced:
                self.config['parameters'][f'Ne_{pop}']=self.config['parameters'][f'Ne_{reference}']

    def _verify_parent(self, infile):
        if infile is None:
            return None
        path = pathlib.Path(infile).resolve()
        parent = path.parent
        if not parent.exists():
            sys.exit("[X] Parent directory not found: %r" % str(infile))
        return str(path)        

    def _parse_config(self, config_file):
        '''validates types in INI config, returns dict with keys, values as sections/params/etc
        - does not deal with missing/incompatible values (should be dealt with in ParameterObj subclasses)
            https://docs.python-cerberus.org/en/stable/usage.html
        '''
        raw_config = configparser.ConfigParser(inline_comment_prefixes='#', allow_no_value=True)
        raw_config.optionxform=str # otherwise keys are lowercase
        raw_config.read(config_file)
        config = {s: dict(raw_config.items(s)) for s in raw_config.sections()}
        
        possible_values_config = configparser.ConfigParser(allow_no_value=True, strict=False, comment_prefixes=None)
        possible_values_config.optionxform=str # otherwise keys are lowercase
        possible_values_config.read(config_file)
        possible_values_dict = {s: dict(possible_values_config.items(s)) for s in possible_values_config.sections()}
        valid_pop_ids = [population.strip(" ") for population in possible_values_dict["populations"]["# possible values reference_pop"].split("|")]
        sample_pop_ids = [population for population in valid_pop_ids if "_" not in population]
        schema = {
            'gimble': {
                'type': 'dict',
                'schema': {
                    'version': {'required': True, 'empty':False, 'type': 'string'},
                    'precision': {'required': True, 'empty':False, 'type': 'integer', 'coerce': int},
                    'random_seed': {'required': True, 'empty':False, 'type': 'integer', 'coerce': int}}},
            'populations': {
                'required': True, 'type': 'dict', 
                'schema': {
                    'A': {'required':True, 'empty':True, 'type':'string'},
                    'B': {'required':True, 'empty':True, 'type':'string'},
                    'reference_pop': {'required':True, 'empty':False, 'isPop':True},
                    'sync_pop_sizes': {'required':False, 'empty':True, 'isPopSync':True},
                    }},
            'k_max': {
                'type': 'dict', 
                'valuesrules': {'required': True, 'empty':False, 'type': 'integer', 'min': 1, 'coerce':int}},
            'mu': {
                'type':'dict', 
                'schema':{
                    'mu': {'required': False, 'empty':True, 'type': 'float', 'coerce':float},
                    'blocklength': {'required': False, 'empty':True, 'type': 'integer', 'coerce':int}
                    }},
            'parameters': {
                'type': 'dict', 'required':True, 'empty':False, 
                'valuesrules': {'coerce':'float_or_list', 'notNone':True}
            },
            'simulations':{
                'required':False, 'empty':True}
            }
        if self._MODULE=='simulate':
            schema['simulations'] = {
                'type': 'dict',
                'required':True,
                'schema': {
                    'ploidy': {'empty':False, 'required':True, 'type':'integer', 'min':1, 'coerce':int},
                    'blocks': {'empty':False, 'type': 'integer', 'min':1, 'coerce':int},
                    'replicates': {'empty': False, 'type': 'integer', 'min':1, 'coerce':int},
                    'sample_size_A': {'empty':False, 'type': 'integer', 'min':1, 'coerce':int},
                    'sample_size_B': {'empty':False, 'type': 'integer', 'min':1, 'coerce':int},
                    'recombination_rate': {'empty': True, 'notNoneFloat':True, 'coerce':'float_or_empty', 'min':0.0},
                    'recombination_map': {'empty': True, 'type': 'string', 'isPath':True, 'dependencies':['number_bins', 'cutoff', 'scale']},
                    'number_bins': {'empty': True, 'notNoneInt': True, 'coerce':'int_or_empty', 'min':1},
                    'cutoff': {'empty': True, 'notNoneFloat': True, 'coerce':'float_or_empty', 'min':0},
                    'scale': {'empty':True, 'type':'string', 'allowed':['lin', 'log']}
                }}
            schema['mu'] = {
                'type':'dict', 
                'schema':{
                    'mu': {'required': True, 'empty':False, 'type': 'float', 'coerce':float},
                    'blocklength': {'required': False, 'empty':True, 'min':1, 'type': 'integer', 'coerce':int}
                    }}

        sync_pops = config["populations"]["sync_pop_sizes"].strip(" ")
        valid_sync_pops = [population.strip(" ") for population in possible_values_dict["populations"]["# possible values sync_pop_sizes"].split("|")]
        if sync_pops in valid_sync_pops:
            equal_pops = sync_pops.strip(' ').split(',')
            if len(equal_pops) == 2:
                first, last = equal_pops
                config['parameters'][f'Ne_{last}'] = config['parameters'][f'Ne_{first}']
            elif len(equal_pops) == 3:
                config['parameters']['Ne_A_B'] = config['parameters']['Ne_A']
                config['parameters']['Ne_B'] = config['parameters']['Ne_A']
            else:
                sys.exit('[X] Provided sync_pop_sizes are invalid.')
        validator = CustomNormalizer(schema, valid_pop_ids=valid_pop_ids, valid_sync_pops=valid_sync_pops)
        validator.validate(config)
        output = ["[X] INI Config file format error(s) ..."]
        if not validator.validate(config):
            print("[X] validator.errors", validator.errors)
            sys.exit()
        config = validator.normalized(config)
        config['populations']['sample_pop_ids'] = sample_pop_ids
        print("[+] Config file validated.")
        return config

    def _process_config(self):
        if self._MODULE in ['makegrid', 'inference', 'simulate']:
            self.config['mu']['blockslength'] = self._get_blocks_length(self.zstore)
            self.config['parameters']['mu'] = self.config['mu']['mu']
            self._expand_params()
            self.reference, self.toBeSynced = self._get_pops_to_sync()
            self._remove_pop_from_dict(self.toBeSynced)
            self.parameter_combinations = self._dict_product()
            self._sync_pop_sizes(self.reference, self.toBeSynced)
        elif self._MODULE=='optimise':
            self.config['mu']['blockslength'] = self._get_blocks_length(self.zstore)
            self.config['parameters']['mu'] = self.config['mu']['mu']
            #parameters either float or [mid, min, max]
            self.reference, self.toBeSynced = self._get_pops_to_sync()
            self._sync_pop_sizes_optimise(self.reference, self.toBeSynced)
            self.parameter_combinations = self._return_boundaries()
            #ready to be scaled
        else:
            sys.exit("[X] gimble.py_processing_config: Not implemented yet.")

    def _return_boundaries(self, length_boundary_set=3):
        parameter_combinations = {k:self._cast_to_repeated_list(v, length_boundary_set)[:length_boundary_set] for k,v in self.config['parameters'].items()}    
        return self._dict_zip(parameter_combinations)

class Store(object):
    def __init__(self, prefix=None, path=None, create=False, overwrite=False):
        self.prefix = prefix if not prefix is None else str(pathlib.Path(path).resolve().stem)
        self.path = path if not path is None else "%s.z" % prefix
        self.data = self._init_data(create, overwrite)
        if create:
            self._init_meta(overwrite=overwrite)

    def tree(self):
        print(self.data.tree())
    
    def log_stage(self, parameterObj):
        '''only remembers last command per module'''
        self.data.attrs[parameterObj._MODULE] = parameterObj._get_cmd()
    
    def get_stage(self, stage):
        return self.data.attrs[stage]

    def has_stage(self, stage):
        return stage in self.data.attrs

    def setup_sim(self, parameterObj):
        print("[#] Preparing store...")
        self._init_meta(overwrite=True)

    def setup_debug(self, parameterObj):
        print("[#] Preparing store...")
        self._init_meta(overwrite=True)
        self._set_bsfs(parameterObj)
        # inference can be done on ETP-counts

    def setup_seq(self, parameterObj):
        print("[#] Preparing store...")
        self._init_meta(overwrite=True)
        print("[#] Processing GENOME_FILE %r..." % parameterObj.genome_f)
        self._set_sequences(parameterObj)
        print("[#] Processing SAMPLE_FILE %r..." % parameterObj.sample_f)
        self._set_samples(parameterObj)
        print("[#] Processing VCF_FILE %r..." % parameterObj.vcf_f)
        self._set_variants(parameterObj)
        print("[#] Processing BED_FILE %r..." % parameterObj.bed_f)
        self._set_intervals(parameterObj)
        self.log_stage(parameterObj)

    def blocks(self, parameterObj):
        print("[#] Preflight...")
        self._preflight_blocks(parameterObj)
        print("[#] Making blocks...")
        #self._make_blocks_threaded(parameterObj)
        self._make_blocks(parameterObj)
        #self.plot_bsfs_pcp(sample_set='X')
        #self._calculate_block_pop_gen_metrics()
        #self._plot_blocks(parameterObj)
        self.log_stage(parameterObj)

    def windows(self, parameterObj):
        print("[#] Preflight...")
        self._preflight_windows(parameterObj)
        print("[#] Making windows...")
        self._make_windows(parameterObj, sample_sets='X')
        self.log_stage(parameterObj)

    def simulate(self, parameterObj):
        print("[#] Preflight...")
        self._preflight_simulate(parameterObj)
        print("[+] Checks passed.")
        lib.simulate.run_sim(parameterObj, self)
        self.log_stage(parameterObj)

    def query(self, parameterObj):
        print("[#] Preflight...")
        self._preflight_query(parameterObj)
        print("[#] Query...")
        if parameterObj.blocks:
            self._write_block_bed(parameterObj, cartesian_only=True)
        if parameterObj.windows:
            self._write_window_bed(parameterObj, cartesian_only=True)

    def _validate_seq_names(self, sequences=None):
        """Returns valid seq_names in sequences or raises ValueError."""
        meta = self.data['seqs'].attrs
        if not sequences:
            return meta['seq_names']
        if set(sequences).issubset(set(meta['seq_names'])):
            return sequences
        else:
            raise ValueError("%s not a subset of %s" % (sequences, meta['seq_names']))

    def _get_invert_population_flag(self, population_by_letter=None):
        """Returns True if populations need inverting, and False if not or population_by_letter is None. 
        Raises ValueError if population_by_letter of data and config differ"""
        meta = self.data['seqs'].attrs
        if population_by_letter:
            if not population_by_letter['A'] in meta['population_by_letter'].values() or not population_by_letter['B'] in meta['population_by_letter'].values():
                raise ValueError("population names in config (%r) and gimble-store (%r) must match" % (string(set(population_by_letter.values())), string(set(meta['population_by_letter'].values()))))
            if not population_by_letter['A'] == meta['population_by_letter']['A']:
                return True
        return False

    def _get_sample_set_idxs(self, query=None):
        """Returns list of sample_set_idxs.

        Parameters 
        ----------
        query : string or None
                'X' - inter-population sample_sets
                'A' - intra-population sample_sets of population A
                'B' - intra-population sample_sets of population B 
                None - all sample_sets
        
        Returns
        -------
        out : list of strings
            sample_set_idxs that can be used to access data in gimble store
        """
        meta = self.data['seqs'].attrs
        if not query:
            return [str(idx) for idx in range(len(meta['sample_sets']))]
        elif query == 'X':
            return [str(idx) for (idx, is_cartesian) in enumerate(meta['sample_sets_inter']) if is_cartesian]
        elif query == 'A':
            return [str(idx) for (idx, is_intra_A) in enumerate(meta['sample_sets_intra_A']) if is_intra_A]
        elif query == 'B':
            return [str(idx) for (idx, is_intra_B) in enumerate(meta['sample_sets_intra_B']) if is_intra_B]
        else:
            raise ValueError("'query' must be 'X', 'A', 'B', or None")

    def _get_window_bsfs(self, sample_sets='X', sequences=None, population_by_letter=None, kmax_by_mutype=None):
        """Return bsfs_array of 5 dimensions (fifth dimension is the window-idx across ALL sequences in sequences).
        [ToDo] Ideally this should work with regions, as in CHR:START-STOP.
        [ToDo] put in sample set context (error if not samples set).
    
        Parameters 
        ----------
        sequences : list of strings or None
            Only make bSFS based on variation on these sequences seq_names.

        population_by_letter : dict (string -> string) or None
            Mapping of population IDs to population letter in model (from INI file).

        kmax : dict (string -> int) or None
            Mapping of kmax values to mutypes.
        
        Returns
        -------
        bsfs : ndarray, int, ndim (1 + mutypes). First dimension is window idx. 
        """
        sequences = self._validate_seq_names(sequences)
        invert_population_flag = self._get_invert_population_flag(population_by_letter)
        max_k = np.array(list(kmax_by_mutype.values())) + 1 if kmax_by_mutype else None 
        variations = []
        for seq_name in tqdm(sequences, total=len(sequences), desc="[%] Querying data ", ncols=100):
            variation = np.array(self.data["seqs/%s/windows/variation" % seq_name], dtype=np.int64)
            variations.append(variation)
        bsfs = np.concatenate(variations, axis=0)
        if invert_population_flag:
            bsfs[0], bsfs[1] = bsfs[1], bsfs[0]
        bsfs = bsfs.reshape((bsfs.shape[0] * bsfs.shape[1], bsfs.shape[2]))
        index = np.repeat(np.arange(variation.shape[0]), variation.shape[1]).reshape(variation.shape[0] * variation.shape[1], 1)
        mutuples, counts = np.unique(
            np.concatenate([index, np.clip(bsfs, 0, max_k)], axis=-1).reshape(-1, bsfs.shape[-1] + 1),
            return_counts=True, axis=0)
        # define out based on max values for each column
        out = np.zeros(tuple(np.max(mutuples, axis=0) + 1), np.int64)
        # assign values
        out[tuple(mutuples.T)] = counts
        return out

    def set_bsfs(self, data_type='etp', bsfs=None):
        '''Saves array in GStore and saves parameters (as strings) in meta. Recovery of parameters is possible via key_to_params().''' 
        key = hashlib.md5(str({k: v for k, v in locals().items() if not k == 'self'}).encode()).hexdigest()
        bsfs_key = 'bsfs/%s' % key
        self.data.create_dataset(bsfs_key, data=bsfs, overwrite=True)

    def get_bsfs(self, data_type=None, sequences=None, sample_sets=None, population_by_letter=None, kmax_by_mutype=None):
        """main method for accessing data (needs error-logic)"""
        key = hashlib.md5(str({k: v for k, v in locals().items() if not k == 'self'}).encode()).hexdigest()
        bsfs_key = 'bsfs/%s' % key
        if bsfs_key in self.data:
            return np.array(self.data[bsfs_key], dtype=np.int64)
        if data_type == 'blocks':
            if not self.has_stage('blocks'):
                sys.exit("[X] No blocks in GStore.")
            bsfs = self._get_block_bsfs(sample_sets=sample_sets, population_by_letter=population_by_letter, kmax_by_mutype=kmax_by_mutype)
        elif data_type == 'windows':
            if not self.has_stage('windows'):
                sys.exit("[X] No blocks in GStore.")
            bsfs = self._get_window_bsfs(population_by_letter=population_by_letter, kmax_by_mutype=kmax_by_mutype)
        else:
            raise ValueError("data_type must be 'blocks' or 'windows'")
        self.data.create_dataset(bsfs_key, data=bsfs, overwrite=True)
        return bsfs

    def _get_block_bsfs(self, sequences=None, sample_sets=None, population_by_letter=None, kmax_by_mutype=None):
        """Returns bsfs_array of 4 dimensions.

        Parameters 
        ----------
        sequences : list of strings or None
            If supplied, bSFSs are based only on variation on those sequences
        
        sample_sets : string or None
                None - all sample_sets (default)
                'X' - inter-population sample_sets
                'A' - intra-population sample_sets of population A
                'B' - intra-population sample_sets of population B 
            If supplied, bSFSs are based only on variation of those sample_sets
        
        population_by_letter : dict (string -> string) or None
            Mapping of population IDs to population letter in model (from INI file).

        kmax : dict (string -> int) or None
            Mapping of kmax values to mutypes.

        Returns
        -------
        out : ndarray, int, ndim (mutypes)

        """
        sample_set_idxs = self._get_sample_set_idxs(query=sample_sets)
        sequences = self._validate_seq_names(sequences)
        invert_population_flag = self._get_invert_population_flag(population_by_letter)
        max_k = np.array(list(kmax_by_mutype.values())) + 1 if kmax_by_mutype else None 
        variations = []
        for seq_name in sequences: 
            for sample_set_idx in sample_set_idxs:
                variation_key = 'seqs/%s/blocks/%s/variation' % (seq_name, sample_set_idx)
                variations.append(np.array(self.data[variation_key], dtype=np.int64))
        variation = np.concatenate(variations, axis=0) # concatenate is faster than offset-indexes
        if invert_population_flag:
            variation[0], variation[1] = variation[1], variation[0]        
        # count mutuples (clipping at k_max, if supplied)
        mutuples, counts = np.unique(np.clip(variation, 0, max_k), return_counts=True, axis=0)
        # define out based on max values for each column
        out = np.zeros(tuple(np.max(mutuples, axis=0) + 1), np.int64)
        # assign values
        out[tuple(mutuples.T)] = counts
        return out

    def _get_setup_report(self, width, blocks=False):
        meta = self.data['seqs'].attrs
        reportObj = ReportObj(width=width)
        reportObj.add_line(prefix="[+]", left='[', center='Setup', right=']', fill='=')
        reportObj.add_line(prefix="[+]", left='seqs')
        right = "%s in %s sequence(s) (n50 = %s)" % (
            format_bases(sum(meta['seq_lengths'])), 
            format_count(len(meta['seq_lengths'])), 
            format_bases(meta['seq_n50']))
        reportObj.add_line(prefix="[+]", branch="T", fill=".", left='genome', right=right)
        right = "%s samples in %s populations" % (format_count(len(meta['samples'])), format_count(len(meta['population_ids'])))
        reportObj.add_line(prefix="[+]", branch="T", left='populations', fill=".", right=right)
        sample_counts_by_population = collections.Counter(meta['populations'])
        for idx, (letter, population_id) in enumerate(meta['population_by_letter'].items()):
            left = "%s = %s" % (letter, population_id)
            right = "%s" % format_count(sample_counts_by_population[population_id])
            branch = "P%s" % ("F" if idx == len(meta['population_by_letter']) -1 else "T")
            reportObj.add_line(prefix="[+]", branch=branch, left=left, right=right)
        reportObj.add_line(prefix="[+]", branch="T", left="sample sets", fill=".", right=format_count(len(meta['sample_sets'])))
        reportObj.add_line(prefix="[+]", branch="PT", left="INTER-population sample-sets (X)", right=format_count(len(self._get_sample_set_idxs(query='X'))))
        reportObj.add_line(prefix="[+]", branch="PT", left="INTRA-population sample-sets (A)", right=format_count(len(self._get_sample_set_idxs(query='A'))))
        reportObj.add_line(prefix="[+]", branch="PF", left="INTRA-population sample-sets (B)", right=format_count(len(self._get_sample_set_idxs(query='B'))))
        reportObj.add_line(prefix="[+]", branch="T", left="variants", fill=".", right=("%s (%s per 1 kb)" % (
               format_count(meta['variants_counts']),
               format_proportion(1000 * meta['variants_counts'] / sum(meta['seq_lengths'])))))
        reportObj.add_line(prefix="[+]", branch="PP", right="".join([c.rjust(8) for c in ["HOMREF", "HOMALT", "HET", "MISS"]]))
        for idx, sample in enumerate(meta['samples']):
            variant_idx = meta['variants_idx_by_sample'][sample]
            branch = "PF" if idx == len(meta['variants_idx_by_sample']) - 1 else "PT"
            left = sample
            right = "%s %s %s %s" % (
            (format_percentage(meta['variants_counts_hom_ref'][variant_idx] / meta['variants_counts']) if meta['variants_counts'] else format_percentage(0)).rjust(7),
            (format_percentage(meta['variants_counts_hom_alt'][variant_idx] / meta['variants_counts']) if meta['variants_counts'] else format_percentage(0)).rjust(7),
            (format_percentage(meta['variants_counts_het'][variant_idx] / meta['variants_counts']) if meta['variants_counts'] else format_percentage(0)).rjust(7),
            (format_percentage(meta['variants_counts_missing'][variant_idx] / meta['variants_counts']) if meta['variants_counts'] else format_percentage(0)).rjust(7))
            reportObj.add_line(prefix="[+]", branch=branch, left=left, right=right)
        reportObj.add_line(prefix="[+]", branch="F", left='intervals', fill=".", 
            right = "%s intervals across %s (%s of genome)" % (
                format_count(meta['intervals_count']),
                format_bases(meta['intervals_span']),
                format_percentage(meta['intervals_span'] / sum(meta['seq_lengths']))))
        return reportObj
    
    def _get_storage_report(self, width):
        reportObj = ReportObj(width=width)
        reportObj.add_line(prefix="[+]", left='[', center='Storage', right=']', fill='=')
        total_size = recursive_get_size(self.path)
        total_right = "%s | %s" % (format_bytes(total_size), format_percentage(total_size/total_size))
        reportObj.add_line(prefix="[+]", left=pathlib.Path(self.path).name, right=total_right, fill='.')
        for idx, group in enumerate(self.data):
            size = recursive_get_size(pathlib.Path(self.path) / pathlib.Path(group))
            percentage = format_percentage(size/total_size)
            branch = branch="T" if idx < len(self.data) - 1 else 'F'
            reportObj.add_line(prefix="[+]", branch=branch, left=group, right='%s | %s' % (format_bytes(size), percentage.rjust(7)), fill='.')
        return reportObj

    def _get_blocks_report(self, width):
        meta = self.data['seqs'].attrs
        reportObj = ReportObj(width=width)
        reportObj.add_line(prefix="[+]", left='[', center='Blocks', right=']', fill='=')
        blocks_count_total = sum(meta['blocks_by_sample_set_idx'].values())
        blocks_raw_count_total = sum(meta['blocks_raw_per_sample_set_idx'].values())
        block_validity = blocks_count_total / blocks_raw_count_total if blocks_count_total else 0
        reportObj.add_line(prefix="[+]", left='blocks')
        reportObj.add_line(prefix="[+]", branch='T', fill=".", 
            left="'-l %s -m %s -u %s -i %s'" % (meta['blocks_length'],meta['blocks_span'],meta['blocks_max_missing'],meta['blocks_max_multiallelic']),
            right=' %s blocks (%s discarded)' % (format_count(blocks_count_total), format_percentage(1 - block_validity)))
        column_just = 14
        reportObj.add_line(prefix="[+]", branch="P", right="".join([c.rjust(column_just) for c in ["X", "A", "B"]]))
        bsfs_X = bsfs_to_2d(self.get_bsfs(data_type='blocks', sample_sets='X'))
        X_hetB, X_hetA, X_hetAB, X_fixed = (np.sum(bsfs_X[:,0] * bsfs_X[:,1]), np.sum(bsfs_X[:,0] * bsfs_X[:,2]), 
                                            np.sum(bsfs_X[:,0] * bsfs_X[:,3]), np.sum(bsfs_X[:,0] * bsfs_X[:,4]))
        bsfs_A = bsfs_to_2d(self.get_bsfs(data_type='blocks', sample_sets='A'))
        A_hetB, A_hetA, A_hetAB, A_fixed = (np.sum(bsfs_A[:,0] * bsfs_A[:,1]), np.sum(bsfs_A[:,0] * bsfs_A[:,2]), 
                                            np.sum(bsfs_A[:,0] * bsfs_A[:,3]), np.sum(bsfs_A[:,0] * bsfs_A[:,4]))
        total_SA = np.sum(bsfs_A[:,0, None] * bsfs_A[:,1:])
        bsfs_B = bsfs_to_2d(self.get_bsfs(data_type='blocks', sample_sets='B'))
        B_hetB, B_hetA, B_hetAB, B_fixed = (np.sum(bsfs_B[:,0] * bsfs_B[:,1]), np.sum(bsfs_B[:,0] * bsfs_B[:,2]), 
                                            np.sum(bsfs_B[:,0] * bsfs_B[:,3]), np.sum(bsfs_B[:,0] * bsfs_B[:,4]))
        total_SB = np.sum(bsfs_A[:,0, None] * bsfs_A[:,1:])
        reportObj.add_line(prefix="[+]", branch='T', left='interval coverage', right="".join(
            [format_percentage(c).rjust(column_just) for c in [
                (meta['blocks_length'] * np.sum(bsfs_X[:,0]) / meta['intervals_span'] / len(self._get_sample_set_idxs("X"))),
                (meta['blocks_length'] * np.sum(bsfs_A[:,0]) / meta['intervals_span'] / len(self._get_sample_set_idxs("A"))),
                (meta['blocks_length'] * np.sum(bsfs_B[:,0]) / meta['intervals_span'] / len(self._get_sample_set_idxs("B")))
                ]]))
        reportObj.add_line(prefix="[+]", branch='T', left='total blocks', right="".join([format_count(c).rjust(column_just) for c in [np.sum(bsfs_X[:,0]), np.sum(bsfs_A[:,0]), np.sum(bsfs_B[:,0])]]))
        reportObj.add_line(prefix="[+]", branch='T', left='invariant blocks', right="".join(
            [format_percentage(c).rjust(column_just) for c in [
                bsfs_X[0,0] / np.sum(bsfs_X[:,0]), 
                bsfs_A[0,0] / np.sum(bsfs_A[:,0]), 
                bsfs_B[0,0] / np.sum(bsfs_B[:,0])
                ]]))
        reportObj.add_line(prefix="[+]", branch='T', left='four-gamete-violation blocks', right="".join(
            [format_percentage(c, precision=2).rjust(column_just) for c in [
                np.sum(bsfs_X[(bsfs_X[:,3]>0) & (bsfs_X[:,4]>0)][:,0]) / np.sum(bsfs_X[:,0]), 
                np.sum(bsfs_A[(bsfs_A[:,3]>0) & (bsfs_A[:,4]>0)][:,0]) / np.sum(bsfs_A[:,0]), 
                np.sum(bsfs_B[(bsfs_B[:,3]>0) & (bsfs_B[:,4]>0)][:,0]) / np.sum(bsfs_B[:,0])]]))
        heterozygosity_XA = (X_hetA + X_hetAB) / (meta['blocks_length'] * np.sum(bsfs_X[:,0]))
        heterozygosity_XB = (X_hetB + X_hetAB) / (meta['blocks_length'] * np.sum(bsfs_X[:,0]))
        dxy_X = ((X_hetA + X_hetB + X_hetAB) / 2.0 + X_fixed) / (meta['blocks_length'] * np.sum(bsfs_X[:,0]))
        mean_pi = (heterozygosity_XA + heterozygosity_XB) / 2.0
        total_pi = (dxy_X + mean_pi) / 2.0 
        fst_X = (dxy_X - mean_pi) / (dxy_X + mean_pi) if (total_pi) else np.nan
        pi_A = float(fractions.Fraction(1, 2) * (A_hetA + A_hetB) + fractions.Fraction(2, 3) * (A_hetAB + A_fixed)) / (meta['blocks_length'] * np.sum(bsfs_A[:,0]))
        watterson_theta_A = total_SA / float(harmonic(3)) / (meta['blocks_length'] * np.sum(bsfs_A[:,0]))
        heterozygosity_A = (A_hetA + A_hetAB) / (meta['blocks_length'] * np.sum(bsfs_A[:,0]))
        pi_B = float(fractions.Fraction(1, 2) * (B_hetA + B_hetB) + fractions.Fraction(2, 3) * (B_hetAB + B_fixed)) / (meta['blocks_length'] * np.sum(bsfs_B[:,0]))
        watterson_theta_B = total_SB / float(harmonic(3)) / (meta['blocks_length'] * np.sum(bsfs_B[:,0]))
        heterozygosity_B = (B_hetA + B_hetAB) / (meta['blocks_length'] * np.sum(bsfs_B[:,0]))
        reportObj.add_line(prefix="[+]", branch='T', left='heterozygosity (A)', right="".join([format_proportion(c, precision=5).rjust(column_just) for c in [heterozygosity_XA,heterozygosity_A, '-']]))
        reportObj.add_line(prefix="[+]", branch='T', left='heterozygosity (B)', right="".join([format_proportion(c, precision=5).rjust(column_just) for c in [heterozygosity_XB,'-', heterozygosity_B]]))
        reportObj.add_line(prefix="[+]", branch='T', left='D_xy', right="".join([format_proportion(c, precision=5).rjust(column_just) for c in [dxy_X,'-', '-']]))
        reportObj.add_line(prefix="[+]", branch='T', left='F_st', right="".join([format_proportion(c, precision=5).rjust(column_just) for c in [fst_X,'-', '-']]))
        reportObj.add_line(prefix="[+]", branch='T', left='Pi', right="".join([format_proportion(c, precision=5).rjust(column_just) for c in ['-',pi_A, pi_B]]))
        reportObj.add_line(prefix="[+]", branch='F', left='Watterson theta', right="".join([format_proportion(c, precision=5).rjust(column_just) for c in ['-',watterson_theta_A, watterson_theta_B]]))
        return reportObj

    def _get_windows_report(self, width):
        meta = self.data['seqs'].attrs
        reportObj = ReportObj(width=width)
        reportObj.add_line(prefix="[+]", left='[', center='Windows', right=']', fill='=')
        reportObj.add_line(prefix="[+]", left='windows')
        reportObj.add_line(prefix="[+]", branch='F', fill=".", left="'-w %s -s %s'" % (meta['window_size'], meta['window_step']), right=' %s windows' % (
            format_count(meta['window_count'])))
        return reportObj

    def info(self, query=None):
        width = 100
        report = self._get_storage_report(width)
        if self.has_stage('setup'):
            report += self._get_setup_report(width)    
        if self.has_stage('blocks'):
            report += self._get_blocks_report(width)
        if self.has_stage('windows'):
            report += self._get_windows_report(width)
        print(report)

    def _count_groups(self, name):
        return len(list(self.data[name]))

    def _init_data(self, create, overwrite):
        if create:
            if os.path.isdir(self.path):
                print("[-] Directory %r already exists." % self.path)
                if not overwrite:
                    print("[X] Please specify '-f' to overwrite.")
                    sys.exit(1)
                print("[+] Deleting existing directory %r" % self.path)
                shutil.rmtree(self.path)
            print("[+] Creating GStore in %r" % self.path)
            return zarr.open(str(self.path), mode='w')
        #print("[+] Loading GStore from %r" % self.path)
        return zarr.open(str(self.path), mode='r+')
    
    def _init_meta(self, overwrite=False):
        attrs_by_group = {
            'seqs' : {
                'vcf_f': None, 
                'sample_f': None, 
                'genome_f': None, 
                'bed_f': None,
                'seq_names': [], 
                'seq_lengths': [], 
                'seq_n50': 0,
                'samples': [], 
                'populations': [], 
                'population_ids': [], 
                'spacing' : 16,
                'sample_sets': [],
                'sample_sets_intra_A': [],
                'sample_sets_intra_B': [],
                'sample_sets_inter': [],
                'population_by_sample': {},
                'population_by_letter': {},
                'variants_counts': [], 
                'variants_idx_by_sample': {}, 
                'variants_counts_hom_ref': [],
                'variants_counts_hom_alt': [],
                'variants_counts_het': [],
                'variants_counts_missing': [],
                'intervals_count': 0, 
                'intervals_span': 0, 
                'intervals_span_sample': [],
                'intervals_idx_by_sample': {},
                'mutypes_count': 4,
                'blocks_length': 0, 
                'blocks_span': 0, 
                'blocks_max_missing': 0, 
                'blocks_max_multiallelic': 0, 
                'blocks_by_sample_set_idx': {},
                'blocks_raw_per_sample_set_idx': {},
                'blocks_maxk_by_mutype': {},
                'window_size': 0, 
                'window_step': 0, 
                'window_count': 0, 
            },
            'inference': {},
            'grids': {},
            'sims': {},
            'bsfs': {}
        }
        for group, attrs in attrs_by_group.items():
            self.data.require_group(group, overwrite=overwrite)
            self.data[group].attrs.put(attrs_by_group[group])

    def _is_zarr_group(self, name, subgroup=None):
        '''needed?'''
        if not subgroup:
            return name in list(self.data.group_keys())
        else:
            return name in list(self.data[subgroup].group_keys())

    def _return_group_last_integer(self, name):
        '''needed?'''
        try:
            all_groups = [int([namestring for namestring in groupnames.split('_')][-1]) for groupnames in list(self.data[name])]
        except KeyError:
            return 0
        if len(all_groups):
            return max(all_groups)+1
        else:
            return 0

    def _set_sequences(self, parameterObj):
        meta = self.data['seqs'].attrs
        sequences_df = parse_csv(
            csv_f=parameterObj.genome_f, 
            sep="\t", 
            usecols=[0,1], 
            dtype={'sequence_id': 'category', 'sequence_length': 'int64'}, 
            header=None)
        meta['seq_names'] = sequences_df['sequence_id'].to_list()
        meta['seq_lengths'] = sequences_df['sequence_length'].to_list()
        meta['seq_n50'] = get_n50_from_lengths(meta['seq_lengths'])
        for sequence_id in meta['seq_names']:
            self.data.create_group('seqs/%s/' % sequence_id)
        meta['genome_f'] = parameterObj.genome_f

    def _set_samples(self, parameterObj):
        meta = self.data['seqs'].attrs
        samples_df = parse_csv(
            csv_f=parameterObj.sample_f, 
            sep=",", 
            usecols=[0,1], 
            dtype={'samples': 'object', 'populations': 'category'}, 
            header=None)
        meta['samples'] = samples_df['samples'].to_list()
        meta['populations'] = samples_df['populations'].to_list()
        meta['population_ids'] = sorted(set(samples_df['populations'].to_list()))
        meta['population_by_letter'] = {letter: population_id for population_id, letter in zip(meta['population_ids'], string.ascii_uppercase)}
        meta['population_by_sample'] = {sample: population for sample, population in zip(meta['samples'], meta['populations'])}
        meta['sample_sets'] = [
            tuple(sorted(x, key=(meta['population_by_sample'].get if meta['population_by_sample'][x[0]] != meta['population_by_sample'][x[1]] else None))) 
                for x in itertools.combinations(meta['population_by_sample'].keys(), 2)]
        longest_sample_string = max([len(", ".join(sample_set)) for sample_set in meta['sample_sets']]) + 2
        meta['spacing'] = longest_sample_string if longest_sample_string > meta['spacing'] else meta['spacing']
        meta['sample_sets_inter'] = [
            False if len(set([meta['population_by_sample'][sample] for sample in sample_set])) == 1 else True 
                for sample_set in meta['sample_sets']]
        meta['sample_sets_intra_A'] = [
            all([meta['population_by_sample'][sample] == meta['population_ids'][0] for sample in sample_set]) for sample_set in meta['sample_sets']]
        meta['sample_sets_intra_B'] = [
            all([meta['population_by_sample'][sample] == meta['population_ids'][1] for sample in sample_set]) for sample_set in meta['sample_sets']]
        meta['sample_f'] = parameterObj.sample_f

    def _set_variants(self, parameterObj):
        meta = self.data['seqs'].attrs
        sequences = meta['seq_names']
        samples = meta['samples']
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            gt_key, pos_key, sample_key = 'calldata/GT', 'variants/POS', 'samples'
            samples_gt_order = allel.read_vcf(parameterObj.vcf_f, fields=[sample_key])[sample_key]
            query_samples = ordered_intersect(a=samples_gt_order, b=samples, order='a')
            # Check if all samples were found
            if set(query_samples) != set(samples):
                sys.exit("[X] The following samples in SAMPLE_FILE were not found in VCF_FILE: %s" % (
                    ", ".join(list(set(samples).difference(set(query_samples))))
                    ))
            # Set up counts arrays
            count_shape = (len(meta['seq_names']), len(query_samples))
            count_records = np.zeros(count_shape[0], dtype=np.int64)
            count_called = np.zeros(count_shape, dtype=np.int64)
            count_hom_ref = np.zeros(count_shape, dtype=np.int64)
            count_hom_alt = np.zeros(count_shape, dtype=np.int64)
            count_het = np.zeros(count_shape, dtype=np.int64)
            count_missing = np.zeros(count_shape, dtype=np.int64)
            for idx, sequence in tqdm(enumerate(sequences), total=len(sequences), desc="[%] Reading variants...", ncols=100):
                vcf_data = allel.read_vcf(parameterObj.vcf_f, region=sequence, samples=query_samples, fields=[gt_key, pos_key])
                if vcf_data:
                    # genotypes
                    gt_matrix = vcf_data[gt_key]
                    # counts
                    sa_genotype_matrix = allel.GenotypeArray(gt_matrix)
                    count_records[idx] = gt_matrix.shape[0]
                    count_called[idx,:] = sa_genotype_matrix.count_called(axis=0)
                    count_hom_ref[idx,:] = sa_genotype_matrix.count_hom_ref(axis=0)
                    count_hom_alt[idx,:] = sa_genotype_matrix.count_hom_alt(axis=0)
                    count_het[idx,:] = sa_genotype_matrix.count_het(axis=0)
                    count_missing[idx,:] = sa_genotype_matrix.count_missing(axis=0)
                    # positions
                    pos_array = check_unique_pos(vcf_data[pos_key] - 1) # port to BED (0-based) coordinates
                    self.data.create_dataset("seqs/%s/variants/pos" % sequence, data=pos_array, dtype=np.int64)
                    self.data.create_dataset("seqs/%s/variants/matrix" % sequence, data=gt_matrix, dtype=np.int64)
            meta['variants_idx_by_sample'] = {query_sample: idx for idx, query_sample in enumerate(query_samples)}
        meta['vcf_f'] = parameterObj.vcf_f
        meta['variants_counts'] = int(np.sum(count_records)) # ZARR JSON encoder does not like numpy dtypes
        meta['variants_counts_called'] = [int(x) for x in np.sum(count_called, axis=0)] # ZARR JSON encoder does not like numpy dtypes
        meta['variants_counts_hom_ref'] = [int(x) for x in np.sum(count_hom_ref, axis=0)] # ZARR JSON encoder does not like numpy dtypes
        meta['variants_counts_hom_alt'] = [int(x) for x in np.sum(count_hom_alt, axis=0)] # ZARR JSON encoder does not like numpy dtypes
        meta['variants_counts_het'] = [int(x) for x in np.sum(count_het, axis=0)] # ZARR JSON encoder does not like numpy dtypes
        meta['variants_counts_missing'] = [int(x) for x in np.sum(count_missing, axis=0)] # ZARR JSON encoder does not like numpy dtypes
        # QC plots 

    def _set_intervals(self, parameterObj):
        meta = self.data['seqs'].attrs
        query_sequences = set(meta['seq_names'])
        df = parse_csv(
            csv_f=parameterObj.bed_f, 
            sep="\t", 
            usecols=[0, 1, 2, 4], 
            dtype={'sequence': 'category', 'start': 'int64', 'end': 'int64', 'samples': 'category'},
            header=None)
        intervals_df = df[df['sequence'].isin(set(meta['seq_names']))].sort_values(['sequence', 'start'], ascending=[True, True]).reset_index(drop=True)
        intervals_df = pd.concat([intervals_df, intervals_df.samples.str.get_dummies(sep=',').filter(meta['samples'])], axis=1).drop(columns=['samples'])
        intervals_df_samples = [sample for sample in intervals_df.columns[3:]]
        query_samples = ordered_intersect(a=intervals_df_samples, b=meta['samples'], order='a')
        intervals_df['length'] = (intervals_df['end'] - intervals_df['start'])
        # Check if all samples were found
        if set(query_samples) != set(meta['samples']):
                sys.exit("[X] The following samples in SAMPLE_FILE were not found in BED_FILE: %s" % (
                    ", ".join(list(set(meta['samples_sorted']).difference(set(query_samples))))))
        # Set up counts arrays
        count_shape = (len(meta['seq_names']), len(query_samples))
        count_bases_samples = np.zeros(count_shape, dtype=np.int64)
        for idx, (sequence, _df) in tqdm(enumerate(intervals_df.groupby(['sequence'], observed=True)), total=len(query_sequences), desc="[%] Reading intervals...", ncols=100):
            interval_matrix = _df[query_samples].to_numpy()
            length_matrix = np.repeat(_df['length'].to_numpy(), interval_matrix.shape[1]).reshape(interval_matrix.shape)
            length_matrix[interval_matrix == 0] = 0 # sets length to 0 if interval not present in interval_matrix 
            count_bases_samples[idx,:] = np.sum(length_matrix, axis=0)
            self.data.create_dataset("seqs/%s/intervals/matrix" % sequence, data=interval_matrix)
            self.data.create_dataset("seqs/%s/intervals/starts" % sequence, data=_df['start'].to_numpy())
            self.data.create_dataset("seqs/%s/intervals/ends" % sequence, data=_df['end'].to_numpy())
        meta['intervals_span_sample'] = [int(x) for x in np.sum(count_bases_samples, axis=0)] # JSON encoder does not like numpy dtypes   
        meta['intervals_count'] = len(intervals_df.index)
        meta['intervals_span'] = int(intervals_df['length'].sum())
        meta['intervals_idx_by_sample'] = {sample: idx for idx, sample in enumerate(query_samples)}
        meta['bed_f'] = parameterObj.bed_f
        # QC plots
        #intervals_df['distance'] = np.where((intervals_df['sequence'] == intervals_df['sequence'].shift(-1)), (intervals_df['start'].shift(-1) - intervals_df['end']) + 1, np.nan)
        #distance_counter = intervals_df['distance'].dropna(how="any", inplace=False).value_counts()
        #length_counter = intervals_df['length'].value_counts()
        #distance_f = "%s.intervals.distance.png" % parameterObj.outprefix
        #plot_loglog(distance_counter, 'Distance to downstream BED interval', distance_f)
        #length_f = "%s.intervals.length.png" % parameterObj.outprefix
        #plot_loglog(length_counter, 'Length of BED interval', length_f)
        #count_sequences = intervals_df['sequence'].nunique()
        #count_intervals = len(intervals_df.index)
        #count_samples = len(query_samples)

    def dump_bsfs(self, parameterObj):
        meta = self.data['seqs'].attrs
        header = ['count'] + [x+1 for x in range(meta['mutypes_count'])]
        bsfs_2d = bsfs_to_2d(self.get_bsfs(data_type='blocks', sample_sets='X'))
        prefix = "%s.l_%s.m_%s.i_%s.u_%s" % (self.prefix, meta['blocks_length'], meta['blocks_span'], meta['blocks_max_missing'], meta['blocks_max_multiallelic'])
        pd.DataFrame(data=bsfs_2d, columns=header, dtype='int64').to_csv("%s.inter.blocks.tsv" % prefix, index=False, sep='\t')
        #bsfs = self.get_bsfs(self, data_type='blocks', sample_sets='A')
        #pd.DataFrame(data=bsfs, columns=header, dtype='int64').to_hdf("%s.intra_A.blocks.h5" % prefix, 'tally', format='table')
        #pd.DataFrame(data=bsfs_2d, columns=header, dtype='int64').to_csv("%s.intra_A.blocks.tsv" % prefix, index=False, sep='\t')
        #bsfs = self.get_bsfs(self, data_type='blocks', sample_sets='B')
        #pd.DataFrame(data=bsfs, columns=header, dtype='int64').to_hdf("%s.intra_B.blocks.h5" % prefix, 'tally', format='table')
        #pd.DataFrame(data=bsfs, columns=header, dtype='int64').to_csv("%s.intra_B.blocks.tsv" % prefix, index=False, sep='\t')

    def _plot_blocks(self, parameterObj):
        # mutuple barchart
        # meta = self.data['seqs'].attrs
        mutypes_inter_key = 'seqs/bsfs/inter/mutypes' 
        counts_inter_key = 'seqs/bsfs/inter/counts' 
        mutypes_inter = self.data[mutypes_inter_key]
        counts_inter = self.data[counts_inter_key]
        self.plot_bsfs_pcp('%s.bsfs_pcp.png' % self.prefix, mutypes_inter, counts_inter)
        
    def _write_block_bed(self, parameterObj, cartesian_only=True):
        '''new gimblestore'''
        meta = self.data['seqs'].attrs
        sample_set_idxs = [idx for (idx, is_cartesian) in enumerate(meta['sample_sets_inter']) if is_cartesian] if cartesian_only else range(len(meta['sample_sets']))
        blocks_count_total = sum([meta['blocks_by_sample_set_idx'][str(idx)] for idx in sample_set_idxs])
        starts = np.zeros(blocks_count_total, dtype=np.int64)
        ends = np.zeros(blocks_count_total, dtype=np.int64)
        # dynamically set string dtype for sequence names
        MAX_SEQNAME_LENGTH = max([len(seq_name) for seq_name in meta['seq_names']])
        sequences = np.zeros(blocks_count_total, dtype='<U%s' % MAX_SEQNAME_LENGTH) 
        sample_sets = np.zeros(blocks_count_total, dtype=np.int64) 
        # if extended_bed
        variation = np.zeros((blocks_count_total, meta['mutypes_count']), dtype=np.int64)
        missing = np.zeros(blocks_count_total, dtype=np.int64) 
        multiallelic = np.zeros(blocks_count_total, dtype=np.int64) 
        with tqdm(total=(len(meta['seq_names']) * len(sample_set_idxs)), desc="[%] Preparing data...", ncols=100, unit_scale=True) as pbar: 
            offset = 0
            for seq_name in meta['seq_names']: 
                for sample_set_idx in sample_set_idxs:
                    start_key = 'seqs/%s/blocks/%s/starts' % (seq_name, sample_set_idx)
                    end_key = 'seqs/%s/blocks/%s/ends' % (seq_name, sample_set_idx)
                    start_array = np.array(self.data[start_key])
                    block_count = start_array.shape[0]
                    starts[offset:offset+block_count] = start_array
                    ends[offset:offset+block_count] = np.array(self.data[end_key])
                    sequences[offset:offset+block_count] = np.full_like(block_count, seq_name, dtype='<U%s' % MAX_SEQNAME_LENGTH)
                    sample_sets[offset:offset+block_count] = np.full_like(block_count, sample_set_idx)
                    if parameterObj.extended_bed:
                        variation_key = 'seqs/%s/blocks/%s/variation' % (seq_name, sample_set_idx)
                        missing_key = 'seqs/%s/blocks/%s/missing' % (seq_name, sample_set_idx)
                        multiallelic_key = 'seqs/%s/blocks/%s/missing' % (seq_name, sample_set_idx)
                        variation[offset:offset+block_count] = np.array(self.data[variation_key])
                        missing[offset:offset+block_count] = np.array(self.data[missing_key]).flatten()
                        multiallelic[offset:offset+block_count] = np.array(self.data[multiallelic_key]).flatten()
                    offset += block_count
                    pbar.update()
        columns = ['sequence', 'start', 'end', 'sample_set']
        if not parameterObj.extended_bed:
            int_bed = np.vstack([starts, ends, sample_sets]).T    
        else:
            int_bed = np.vstack([starts, ends, sample_sets, missing, multiallelic, variation.T]).T
            mutypes_count = ["m_%s" % str(x+1) for x in range(meta['mutypes_count'])]
            columns += ['missing', 'multiallelic'] + mutypes_count    
        # header
        header = ["# %s" % parameterObj._VERSION]
        header += ["# %s = %s" % (sample_set_idx, ", ".join(meta['sample_sets'][sample_set_idx])) for sample_set_idx in sample_set_idxs] 
        header += ["# %s" % "\t".join(columns)]  
        out_f = '%s.blocks.bed' % self.prefix
        with open(out_f, 'w') as out_fh:
            out_fh.write("\n".join(header) + "\n")
        # bed
        bed_df = pd.DataFrame(data=int_bed, columns=columns[1:])
        bed_df['sequence'] = sequences
        bed_df.sort_values(['sequence', 'start'], ascending=[True, True]).to_csv(out_f, mode='a', sep='\t', index=False, header=False, columns=columns)

    def _preflight_query(self, parameterObj):
        if parameterObj.blocks:
            if not self.has_stage('blocks'):
                sys.exit("[X] GStore %r has no blocks. Please run 'gimble blocks'." % self.path)
        if parameterObj.windows:
            if not self.has_stage('windows'):
                sys.exit("[X] GStore %r has no windows. Please run 'gimble windows'." % self.path)

    def _preflight_windows(self, parameterObj):
        if not self.has_stage('blocks'):
            sys.exit("[X] GStore %r has no blocks. Please run 'gimble blocks'." % self.path)
        if self.has_stage('windows'):
            if not parameterObj.overwrite:
                sys.exit("[X] GStore %r already contains windows.\n[X] These windows => %r\n[X] Please specify '--force' to overwrite." % (self.path, self.get_stage('windows')))
            print('[-] GStore %r already contains windows. But these will be overwritten...' % (self.path))
    
    def _preflight_blocks(self, parameterObj):
        if not self.has_stage('setup'):
            sys.exit("[X] GStore %r has no data. Please run 'gimble setup'." % self.path)
        if self.has_stage('blocks'):
            if not parameterObj.overwrite:
                sys.exit("[X] GStore %r already contains blocks.\n[X] These blocks => %r\n[X] Please specify '--force' to overwrite." % (self.path, self.get_stage('blocks')))
            print('[-] GStore %r already contains blocks. But these will be overwritten...' % (self.path))
            # wipe bsfs, since new blocks....
            self.data.create_group("bsfs/", overwrite=True)

    def _preflight_simulate(self, parameterObj):
        if 'sims' not in self.data.group_keys():
            self._init_meta(overwrite=False, module='sims')

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

    def _make_windows(self, parameterObj, sample_sets='X'):
        meta = self.data['seqs'].attrs
        meta['window_size'] = parameterObj.window_size
        meta['window_step'] = parameterObj.window_step
        meta['window_count'] = 0
        sample_set_idxs = self._get_sample_set_idxs(sample_sets)
        with tqdm(meta['seq_names'], total=(len(meta['seq_names']) * len(sample_set_idxs)), desc="[%] Constructing windows ", ncols=100, unit_scale=True) as pbar: 
            for seq_name in meta['seq_names']:
                variation, starts, ends = [], [], []
                for sample_set_idx in sample_set_idxs:
                    variation_key = 'seqs/%s/blocks/%s/variation' % (seq_name, sample_set_idx)
                    variation.append(np.array(self.data[variation_key]))
                    start_key = 'seqs/%s/blocks/%s/starts' % (seq_name, sample_set_idx)
                    starts.append(np.array(self.data[start_key]))
                    end_key = 'seqs/%s/blocks/%s/ends' % (seq_name, sample_set_idx)
                    ends.append(np.array(self.data[end_key]))
                    pbar.update()
                variation_array = np.concatenate(variation, axis=0)
                start_array = np.concatenate(starts, axis=0)
                end_array = np.concatenate(ends, axis=0)
                # window_variation : shape = (windows, blocklength, 4)
                window_variation, window_starts, window_ends, window_pos_mean, window_pos_median = cut_windows(variation_array, sample_set_idxs, start_array, end_array, num_blocks=parameterObj.window_size, num_steps=parameterObj.window_step)
                #b, counts = np.unique(variation, return_counts=True, axis=0)
                self.data.create_dataset("seqs/%s/windows/variation" % seq_name, data=window_variation, overwrite=True)
                self.data.create_dataset("seqs/%s/windows/starts" % seq_name, data=window_starts, overwrite=True)
                self.data.create_dataset("seqs/%s/windows/ends" % seq_name, data=window_ends, overwrite=True)
                self.data.create_dataset("seqs/%s/windows/pos_mean" % seq_name, data=window_pos_mean, overwrite=True)
                self.data.create_dataset("seqs/%s/windows/pos_median" % seq_name, data=window_pos_median, overwrite=True)
                meta['window_count'] += window_variation.shape[0]
        #window_info_rows = []
        #window_mutuple_tally = []
        ## window bsfs
        #for seq_name in tqdm(meta['seq_names'], total=len(meta['seq_names']), desc="[%] Generating output ", ncols=100):
        #    variations = self.data["seqs/%s/windows/variation" % seq_name]
        #    print('variations', variations.shape, variations.info, variations.hexdigest())
        #    for window_id, variation, midpoint_mean, midpoint_median in zip(window_ids, variations, midpoint_means, midpoint_medians):
        #        pi_1, pi_2, d_xy, f_st, fgv = calculate_popgen_from_array(variation, (meta['blocks_length'] * variation.shape[0]))
        #        blocks = variation.shape[0]
        #        window_info_rows.append([window_id, seq_name, midpoint_mean, midpoint_median, pi_1, pi_2, d_xy, f_st, fgv, blocks])
        #        # bsfs for each window
        #        mutuple, counts = np.unique(variation, return_counts=True, axis=0)
        #        tally = np.concatenate([counts[:, np.newaxis], mutuple], axis =-1)
        #        windows = np.array([window_id] * tally.shape[0])
        #        window_mutuple_tally.append(np.concatenate([windows[:, np.newaxis], tally], axis =-1))
        #window_bsfs_cols = ['window_id', 'count'] + [x+1 for x in range(meta['mutypes_count'])]
        #print(window_bsfs_cols)
        #print(window_mutuple_tally)
        #window_bsfs_df = pd.DataFrame(np.vstack(window_mutuple_tally), columns=window_bsfs_cols)
        #print("[+] Made %s windows" % window_bsfs_df['window_id'].nunique()) 
        #window_bsfs_f = "%s.window_bsfs.tsv" % self.prefix
        #window_bsfs_df.to_csv(window_bsfs_f, sep='\t', index=False)
        #print("[>] Created: %r" % str(window_bsfs_f))
        #window_info_cols = ['window_id', 'seq_name', 'midpoint_mean', 'midpoint_median', 'pi_%s' % meta['population_ids'][0], 'pi_%s' % meta['population_ids'][1], 'd_xy', 'f_st', 'f']
        #window_info_df = pd.DataFrame(window_info_rows, columns=window_info_cols)
        #window_info_f = "%s.window_info.tsv" % self.prefix
        #window_info_df.to_csv(window_info_f, sep='\t', index=False)
        #print("[>] Created: %r" % str(window_info_f))
        #self.plot_fst_genome_scan(window_info_df)
        #self.plot_pi_genome_scan(window_info_df)
        ##     plot_pi_scatter(window_df, '%s.pi_scatter.png' % parameterObj.dataset)

    def _get_interval_coordinates_for_sample_set(self, seq_name='', sample_set=[]):
        meta = self.data['seqs'].attrs
        sample_set_key = np.array([meta['intervals_idx_by_sample'][sample] for sample in sample_set]) # order should not matter...
        matrix_key = 'seqs/%s/intervals/matrix' % seq_name
        start_key = 'seqs/%s/intervals/starts' % seq_name
        end_key = 'seqs/%s/intervals/ends' % seq_name
        if matrix_key in self.data:
            mask = np.all(np.array(self.data[matrix_key])[:,sample_set_key], axis=1)
            return (np.array(self.data[start_key])[mask], np.array(self.data[end_key])[mask])
        return None

    #def _make_blocks(self, parameterObj, debug=False):
    #    meta = self.data['seqs'].attrs
    #    meta['blocks_length'] = parameterObj.block_length
    #    meta['blocks_span'] = parameterObj.block_span
    #    meta['blocks_gap_run'] = parameterObj.block_gap_run
    #    meta['blocks_max_missing'] = parameterObj.block_max_missing
    #    meta['blocks_max_multiallelic'] = parameterObj.block_max_multiallelic
    #    blocks_raw_per_sample_set_idx = collections.Counter()   # all possible blocks
    #    blocks_per_sample_set_idx = collections.Counter()       # all valid blocks => only these get saved to store
#
    #    with tqdm(total=(len(meta['seq_names']) * len(meta['sample_sets'])), desc="[%] Building blocks ", ncols=100, unit_scale=True) as pbar:        
    #        for seq_name in meta['seq_names']:        
    #            block_site_arrays = []
    #            for sample_set_idx, sample_set in enumerate(meta['sample_sets']):
    #                starts, ends = self._get_interval_coordinates_for_sample_set(seq_name=seq_name, sample_set=sample_set)
    #                # Cut sample-set specific blocks based on intervals and block-algoritm parameters
    #                block_site_array = cut_blocks(starts, ends, meta['blocks_length'], meta['blocks_span'], meta['blocks_gap_run'])
    #                print('block_site_array', block_site_array.shape)
    #                block_site_arrays.append(block_site_array)
    #                pbar.update(1)
    #            block_sites = np.concatenate(block_site_arrays, axis=0)
    #            print('block_sites', block_sites.shape)
    #            print(block_sites)
    #    sys.exit('.')

    def _make_blocks(self, parameterObj, debug=False):
        meta = self.data['seqs'].attrs
        meta['blocks_length'] = parameterObj.block_length
        meta['blocks_span'] = parameterObj.block_span
        meta['blocks_gap_run'] = parameterObj.block_gap_run
        meta['blocks_max_missing'] = parameterObj.block_max_missing
        meta['blocks_max_multiallelic'] = parameterObj.block_max_multiallelic
        blocks_raw_per_sample_set_idx = collections.Counter()   # all possible blocks
        blocks_per_sample_set_idx = collections.Counter()       # all valid blocks => only these get saved to store
        parse = True
        with tqdm(total=(len(meta['seq_names']) * len(meta['sample_sets'])), desc="[%] Building blocks ", ncols=100, unit_scale=True) as pbar:        
            for seq_name in meta['seq_names']:
                #if seq_name == 'monCan3F9.scf.00112':
                #    parse=True
                if parse == True:
                    pos_key = "seqs/%s/variants/pos" % (seq_name)
                    gt_key = "seqs/%s/variants/matrix" % (seq_name)
                    pos = np.array(self.data[pos_key], dtype=np.int64) if pos_key in self.data else None
                    sa_genotype_array = allel.GenotypeArray(self.data[gt_key].view(read_only=True)) if gt_key in self.data else None
                    #print('\nseq_name', seq_name)
                    #print('sa_genotype_array', sa_genotype_array)
                    for sample_set_idx, sample_set in enumerate(meta['sample_sets']):
                        print("\n", seq_name, sample_set_idx, sample_set)
                        start_end = self._get_interval_coordinates_for_sample_set(seq_name=seq_name, sample_set=sample_set)
                        if not start_end is None:
                            # Cut sample-set specific blocks based on intervals and block-algoritm parameters
                            starts, ends = start_end
                            block_sites = cut_blocks(starts, ends, meta['blocks_length'], meta['blocks_span'], meta['blocks_gap_run'])
                            #print('block_sites', block_sites.shape, block_sites)
                            #print('block_sites', type(block_sites), block_sites)
                            #print('np.any(block_sites)', np.any(block_sites))
                            #print('not block_sites is None', not block_sites is None)
                            if not block_sites is None and np.any(block_sites):
                                # Allocate starts/ends before overwriting position ints
                                block_starts = np.array(block_sites[:, 0], dtype=np.int64)
                                block_ends = np.array(block_sites[:, -1] + 1, dtype=np.int64)
                                # variants take longer than blocking
                                if np.any(pos) or pos is not None:
                                    ##print('pos', pos.shape, pos)
                                    idx_pos_in_block_sites = np.isin(pos, block_sites, assume_unique=True)
                                    #print('idx_pos_in_block_sites', idx_pos_in_block_sites)
                                    if np.any(idx_pos_in_block_sites):
                                        sample_set_vcf_idxs = [meta['variants_idx_by_sample'][sample] for sample in sample_set]
                                        idx_block_sites_in_pos = np.isin(block_sites, pos, assume_unique=True) 
                                        sa_sample_set_genotype_array = sa_genotype_array.subset(idx_pos_in_block_sites, sample_set_vcf_idxs)
                                        block_sites = genotype_to_mutype_array(sa_sample_set_genotype_array, idx_block_sites_in_pos, block_sites, debug)
                                    else:
                                        block_sites[:] = 2 # if no variants, set all to invariant    
                                else:
                                    block_sites[:] = 2 # if no variants, set all to invariant
                                multiallelic, missing, monomorphic, variation = block_sites_to_variation_arrays(block_sites)
                                valid = (np.less_equal(missing, meta['blocks_max_missing']) & np.less_equal(multiallelic, meta['blocks_max_multiallelic'])).flatten()
                                blocks_raw_per_sample_set_idx[sample_set_idx] += valid.shape[0]
                                blocks_per_sample_set_idx[sample_set_idx] += valid[valid==True].shape[0]
                                blocks_starts_key = 'seqs/%s/blocks/%s/starts' % (seq_name, sample_set_idx)
                                self.data.create_dataset(blocks_starts_key, data=block_starts[valid], overwrite=True)
                                blocks_ends_key = 'seqs/%s/blocks/%s/ends' % (seq_name, sample_set_idx)
                                self.data.create_dataset(blocks_ends_key, data=block_ends[valid], overwrite=True)
                                blocks_variation_key = 'seqs/%s/blocks/%s/variation' % (seq_name, sample_set_idx)
                                self.data.create_dataset(blocks_variation_key, data=variation[valid], overwrite=True)
                                blocks_missing_key = 'seqs/%s/blocks/%s/missing' % (seq_name, sample_set_idx)
                                self.data.create_dataset(blocks_missing_key, data=missing[valid], overwrite=True)
                                blocks_multiallelic_key = 'seqs/%s/blocks/%s/multiallelic' % (seq_name, sample_set_idx)
                                self.data.create_dataset(blocks_multiallelic_key, data=multiallelic[valid], overwrite=True)
                                #print("[*] Sample_set runtime: %.3fs" % (timer() - sample_set_start_time))
                        pbar.update(1)
        meta['blocks_by_sample_set_idx'] = dict(blocks_per_sample_set_idx) # keys are strings
        meta['blocks_raw_per_sample_set_idx'] = dict(blocks_raw_per_sample_set_idx) # keys are strings

    def _make_blocks_threaded(self, parameterObj, debug=False, threads=2):
        '''there might be some speed improvement here ... has to be finished...'''
        meta = self.data['seqs'].attrs
        meta['blocks_length'] = parameterObj.block_length
        meta['blocks_span'] = parameterObj.block_span
        meta['blocks_gap_run'] = parameterObj.block_gap_run
        meta['blocks_max_missing'] = parameterObj.block_max_missing
        meta['blocks_max_multiallelic'] = parameterObj.block_max_multiallelic
        blocks_raw_per_sample_set_idx = collections.Counter()   # all possible blocks
        blocks_per_sample_set_idx = collections.Counter()       # all valid blocks => only these get saved to store
        params = [(seqs, str(sample_seq_idx)) for seqs, sample_seq_idx in itertools.product(meta['seq_names'], range(0, len(meta['sample_sets'])))]
        print(params)
        results = []
        with poolcontext(processes=threads) as pool:
            with tqdm(total=len(params), desc="[%] ", ncols=100, unit_scale=True) as pbar:
                for blockObjs in pool.imap_unordered(block_algorithm, params):
                    results.append(blockObjs)
                    pbar.update()

        with tqdm(total=(len(meta['seq_names']) * len(meta['sample_sets'])), desc="[%] Building blocks ", ncols=100, unit_scale=True) as pbar:        
            for seq_name in meta['seq_names']:        
                pos_key = "seqs/%s/variants/pos" % (seq_name)
                gt_key = "seqs/%s/variants/matrix" % (seq_name)
                for sample_set_idx, (sample_set, sample_set_cartesian) in enumerate(zip(meta['sample_sets'], meta['sample_sets_inter'])):
                    params.append(seq_name, sample_set_idx, pos_key, gt_key)
                    pbar.update(1)
        meta['blocks_by_sample_set_idx'] = dict(blocks_per_sample_set_idx) # keys are strings
        meta['blocks_raw_per_sample_set_idx'] = dict(blocks_raw_per_sample_set_idx) # keys are strings

    def plot_bsfs_pcp(self, sample_set='X'):
        '''plots a bsfs pcp for a given sample_set. 
        '''
        # https://stackoverflow.com/questions/8230638/parallel-coordinates-plot-in-matplotlib
        # http://www.shengdongzhao.com/publication/tracing-tuples-across-dimensions-a-comparison-of-scatterplots-and-parallel-coordinate-plots/
        #meta = self.data['seqs'].attrs
        bsfs = bsfs_to_2d(self._get_block_bsfs(sample_sets='X'))
        print(bsfs)
        freq = bsfs[:,0] / np.sum(bsfs[:,0])
        meta = self.data['seqs'].attrs
        x = ['m_1', 'm_2', 'm_3', 'm_4']
        bins = np.linspace(0, 1, 9)
        freq = bins[np.digitize(freq, bins)]
        data = bsfs[:,1:] 
        cmap = plt.get_cmap("Greys")
        print('data', data.shape, data)
        print('freq', freq.shape, freq)
        fig, axes = plt.subplots(1, 3, sharey=False, figsize=(15,5))
        axes[0].plot(x,data.T, c=cmap(freq))
        axes[1].plot(x,data.T, c=cmap(freq))
        axes[2].plot(x,data.T, c=cmap(freq))
        plt.subplots_adjust(wspace=0)
        #     min_max_range[mutype] = [np.min(, df[col].max(), np.ptp(df[col])]
        #         df[col] = np.true_divide(df[col] - df[col].min(), np.ptp(df[col]))
        # ynames = ['m_%s' for idx in range(1, meta['mutypes_count'] + 1)]
        # import pandas as pd
        # from pandas.tools.plotting import parallel_coordinates
        # ax = pd.tools.plotting.parallel_coordinates(mutypes)

        # #ax.plot(window_df['rel_pos'], window_df[pi_A_key], color='orange', alpha=0.8, linestyle='-', linewidth=1, label=pi_A_key.replace('pi_', ''))
        # #ax.plot(window_df['rel_pos'], window_df[pi_B_key], color='dodgerblue', alpha=0.8, linestyle='-', linewidth=1, label=pi_B_key.replace('pi_', ''))
        # #y_lim = (min(window_df[pi_A_key].min(), window_df[pi_B_key].min()), max(window_df[pi_A_key].max(), window_df[pi_B_key].max()))
        # #ax.vlines(x_boundaries, y_lim[0], y_lim[1], colors=['lightgrey'], linestyles='dashed', linewidth=1)
        # #ax.set_ylim(y_lim)
        # #ax.spines['right'].set_visible(False)
        # #ax.spines['top'].set_visible(False)
        # #ax.legend(numpoints=1)
        # #plt.ylabel('Pi')
        # #plt.xlabel("Genome coordinate")
        out_f = '%s.pcp.png' % self.prefix
        plt.savefig(out_f, format="png")
        #print("[>] Created: %r" % str(out_f))
        #plt.close(fig)

        #fig, host = plt.subplots()
        #
        ## create some dummy data
        #ynames = ['P1', 'P2', 'P3', 'P4', 'P5']
        #N1, N2, N3 = 10, 5, 8
        #N = N1 + N2 + N3
        #category = np.concatenate([np.full(N1, 1), np.full(N2, 2), np.full(N3, 3)])
        #y1 = np.random.uniform(0, 10, N) + 7 * category
        #y2 = np.sin(np.random.uniform(0, np.pi, N)) ** category
        #y3 = np.random.binomial(300, 1 - category / 10, N)
        #y4 = np.random.binomial(200, (category / 6) ** 1/3, N)
        #y5 = np.random.uniform(0, 800, N)
        #
        ## organize the data
        #ys = np.dstack([y1, y2, y3, y4, y5])[0]
        #ymins = ys.min(axis=0)
        #ymaxs = ys.max(axis=0)
        #dys = ymaxs - ymins
        #ymins -= dys * 0.05  # add 5% padding below and above
        #ymaxs += dys * 0.05
        #dys = ymaxs - ymins
        #
        ## transform all data to be compatible with the main axis
        #zs = np.zeros_like(ys)
        #zs[:, 0] = ys[:, 0]
        #zs[:, 1:] = (ys[:, 1:] - ymins[1:]) / dys[1:] * dys[0] + ymins[0]
        #
        #
        #axes = [host] + [host.twinx() for i in range(ys.shape[1] - 1)]
        #for i, ax in enumerate(axes):
        #    ax.set_ylim(ymins[i], ymaxs[i])
        #    ax.spines['top'].set_visible(False)
        #    ax.spines['bottom'].set_visible(False)
        #    if ax != host:
        #        ax.spines['left'].set_visible(False)
        #        ax.yaxis.set_ticks_position('right')
        #        ax.spines["right"].set_position(("axes", i / (ys.shape[1] - 1)))
        #
        #host.set_xlim(0, ys.shape[1] - 1)
        #host.set_xticks(range(ys.shape[1]))
        #host.set_xticklabels(ynames, fontsize=14)
        #host.tick_params(axis='x', which='major', pad=7)
        #host.spines['right'].set_visible(False)
        #host.xaxis.tick_top()
        #host.set_title('Parallel Coordinates Plot', fontsize=18)
        #
        #colors = plt.cm.tab10.colors
        #for j in range(N):
        #    # to just draw straight lines between the axes:
        #    # host.plot(range(ys.shape[1]), zs[j,:], c=colors[(category[j] - 1) % len(colors) ])
        #
        #    # create bezier curves
        #    # for each axis, there will a control vertex at the point itself, one at 1/3rd towards the previous and one
        #    #   at one third towards the next axis; the first and last axis have one less control vertex
        #    # x-coordinate of the control vertices: at each integer (for the axes) and two inbetween
        #    # y-coordinate: repeat every point three times, except the first and last only twice
        #    verts = list(zip([x for x in np.linspace(0, len(ys) - 1, len(ys) * 3 - 2, endpoint=True)],
        #                     np.repeat(zs[j, :], 3)[1:-1]))
        #    # for x,y in verts: host.plot(x, y, 'go') # to show the control points of the beziers
        #    codes = [Path.MOVETO] + [Path.CURVE4 for _ in range(len(verts) - 1)]
        #    path = Path(verts, codes)
        #    patch = patches.PathPatch(path, facecolor='none', lw=1, edgecolor=colors[category[j] - 1])
        #    host.add_patch(patch)
        #plt.tight_layout()
        #plt.show()

# class Store(object):

#     def _from_zarr(self, parameterObj):
#         '''should be made more explicit'''
#         self.path = str(parameterObj.zstore)
#         self.prefix = str(parameterObj.zstore).rstrip(".z")
#         self.data = zarr.open(self.path, mode='r+')

#     def _from_input(self, parameterObj):
#         self.prefix = parameterObj.outprefix
#         self.path = self._get_path(parameterObj.outprefix)
#         self.data = zarr.open(self.path, mode='w')

#         self._parse_genome_file(parameterObj.genome_file)
#         self._parse_sample_file(
#             str(parameterObj.sample_file), 
#             parameterObj._pairedness
#             )
#         self._parse_vcf_file(str(parameterObj.vcf_file))
#         self.data.attrs['bed_f'] = str(parameterObj.bed_file)
#         self.add_stage(parameterObj)
        
#     def get_base_count(self, sequence_id=None, kind='total'):
#         if kind == 'total':
#             if sequence_id is None:
#                 return sum(self.data.attrs['sequence_length'])
#             else:
#                 if isinstance(sequence_id, list):
#                     sequence_ids = set(sequence_id)
#                 elif isinstance(sequence_id, str):
#                     sequence_ids = set([sequence_id])
#                 else:
#                     return None
#                 return sum([s_length for s_id, s_length in zip(
#                     self.data.attrs['sequence_length'], 
#                     self.data.attrs['sequence_ids']) if s_id in sequence_ids])

#     def get_pop_ids_by_sample_id(self):
#         return self.data.attrs['pop_ids_by_sample_id']

#     def info(self, verbose=False):
#         pop_ids_by_sample_id = self.get_pop_ids_by_sample_id()
#         pop_ids = self.data.attrs['pop_ids']
#         print("[=] ==========================================")
#         print("[+] [DataStore   ] %s" % self.path)
#         print("[+] [Genome      ] %s in %s sequence(s)" % (format_bases(self.get_base_count(kind='total')), len(self.data.attrs['sequence_ids'])))
#         print("[+] [Samples     ] %s in %s populations" % (len(pop_ids_by_sample_id), len(pop_ids)))
        
#         if verbose:
#             print("%s" % "\n".join(["[+] \t %s [%s]" % (sample_id, pop_id) for sample_id, pop_id in pop_ids_by_sample_id.items()]))
#         print("[+] [Sample-sets ] %s" % (len(self.data.attrs['idx_cartesian_sample_sets'])))
        
#         #self.data.attrs['block_length']
#         #self.data.attrs['block_span']
#         #print("[+] [  Variants  ] : %s variants" % (sum(self.data.attrs['variants'])))
#         #interval_counts = []
#         #interval_sites = []
#         #block_sites = []
#         #for sequence_id in self.data.attrs['sequence_ids']:
#         #    for sample_set_idx in self.data.attrs['sequence_ids']
#         #    interval_counts.append(self.data.attrs['block_span']) 
#         #    interval_sites.append(self.data.attrs['block_span'])
#         #    block_sites.append()

#         #print("[+] [  Intervals ] : %s in %sb" % (self.data.attrs['pairedness'], len(self.data.attrs['idx_cartesian_sample_sets'])))
#         #print("[+] [  Blocks    ] l=%s s=%s : %s" % (self.data.attrs['pairedness'], len(self.data.attrs['idx_cartesian_sample_sets'])))

#         #self.data.attrs['window_size']
#         #self.data.attrs['window_step']
# #
#         #for sequence_id in self.data.attrs['sequence_ids']:
#         #    window = self.data["%s/windows/window_id" % sequence_id]
#         #self.data.create_dataset("%s/windows/window_id")
#         #self.data.create_dataset("%s/windows/midpoint_mean")
#         #self.data.create_dataset("%s/windows/midpoint_median")
#         #self.data.create_dataset("%s/%s/interval_sites")
#         #self.data.create_dataset("%s/%s/block_sites")
#         #self.data.create_dataset("%s/%s/blocks/starts")
#         #self.data.create_dataset("%s/%s/blocks/ends")
#         #self.data.create_dataset("%s/%s/blocks/variation")
#         #self.data.create_dataset("%s/%s/blocks/multiallelic")
#         #self.data.create_dataset("%s/%s/blocks/missing")
#         #for sequence_id in self.data.attrs['sequence_ids']:

#         #print("[+] [Blocks] = %s: %s" % (self.data.attrs['pairedness'], len(self.data.attrs['idx_cartesian_sample_sets'])))
#         #print("[+] [Windows] = %s: %s" % (self.data.attrs['pairedness'], len(self.data.attrs['idx_cartesian_sample_sets'])))
#         #print("[+] [Grids] = %s: %s" % (self.data.attrs['pairedness'], len(self.data.attrs['idx_cartesian_sample_sets'])))
#         print("[=] ==========================================")

#     def tree(self):
#         return self.data.tree()

#     def attrs(self):
#         print([self.data.attrs['sample_sets'][idx] for idx in self.data.attrs['idx_cartesian_sample_sets']])
#         return "\n".join(
#             ["\t".join([k, str(len(v)), str(type(v)), str(v)]) if isinstance(v, list) else "\t".join([k, str(v), str(type(v)), str(v)]) for k, v in self.data.attrs.asdict().items()])


#     def _parse_vcf_file(self, vcf_file):
#         sample_ids = self.data.attrs['sample_ids']
#         sequence_ids = self.data.attrs['sequence_ids']
#         variant_counts = []
#         with warnings.catch_warnings():
#             warnings.simplefilter("ignore")
#             for seq_id in tqdm(sequence_ids, total=len(sequence_ids), desc="[%] Parsing input files ", ncols=100):
#                 zarr_key = ''
#                 #print(vcf_file)
#                 data_by_key = allel.read_vcf(vcf_file, region=seq_id, samples=sample_ids, fields=['samples', 'calldata/GT', 'variants/POS'])
#                 if data_by_key:
#                     for key, data in data_by_key.items():
#                         if key == 'samples':
#                             self.data.attrs['sample_ids_vcf'] = list(data)
#                             #print("vcf samples", list(data))
#                             self.data.attrs['sample_ids_to_vcf_idx'] = {sample_id: idx for idx, sample_id in enumerate(data)}
#                         elif key == 'calldata/GT':
#                             zarr_key = "%s/gt" % seq_id
#                             self.data.create_dataset(zarr_key, data=data) 
#                         elif key == 'variants/POS':
#                             zarr_key = "%s/pos" % seq_id
#                             data = np.array(data)
#                             variant_counts.append(data.shape[0])
#                             #print('np.array(data).shape', np.array(data).shape)
#                             unique_pos, counts_pos = np.unique(data, return_counts=True)
#                             duplicates = unique_pos[counts_pos > 1]
#                             if duplicates.any():
#                                 print("\n[-] Sequence %r: %s VCF records with non-unique positions found. Rescuing records by shifting position... (abort if this is not desired)" % (seq_id, len(duplicates)))
#                                 data = fix_pos_array(data)
#                             self.data.create_dataset(zarr_key, data=data - 1) # port to BED (0-based) coordinates
#                         else:
#                             logging.error("[X] Unknown key %r" % key)
#                             sys.exit()
#             self.data.attrs['variants'] = variant_counts
#             self.data.attrs['vcf_f'] = str(vcf_file)

#     def has_data(self, datatype):
#         if datatype in self.data.attrs:
#             return True
#         return False

#     def make_blocks(self, parameterObj, debug=False):
#         if self.has_data('blocks'):
#             if not parameterObj.overwrite:
#                 sys.exit("[X] Store %r already contains blocks.\n[X] These blocks => %r\n[X] Please specify '--force' to overwrite." % (self.path, self.get_stage_cmd('blocks')))
#             print('[-] Store %r already contains blocks. But these will be overwritten...' % (self.path))
#         self.data.attrs['block_length'] = parameterObj.block_length
#         self.data.attrs['block_gap_run'] = parameterObj.block_gap_run
#         self.data.attrs['block_span'] = parameterObj.block_span
#         self.data.attrs['block_max_missing'] = parameterObj.block_max_missing
#         self.data.attrs['block_max_multiallelic'] = parameterObj.block_max_multiallelic
#         self.data.attrs['mutypes_count'] = 4 # should be calculated from possible folded genotypes and pairedness
#         sample_ids = self.data.attrs['sample_ids'] # from BED file
#         sequence_ids = self.data.attrs['sequence_ids']
#         sample_sets = self.data.attrs['sample_sets']
#         #print(self.attrs())
#         logging.info("[#] Processing BED file %r ..." % self.data.attrs['bed_f'])
#         df = pd.read_csv(self.data.attrs['bed_f'], sep="\t", usecols=[0, 1, 2, 4], names=['sequence_id', 'start', 'end', 'samples'], dtype={'sequence_id': str, 'start': np.int, 'end': np.int, 'samples': str})
#         # remove sequence_ids that are not in sequence_names_array, sort, reset index
#         intervals_df = df[df['sequence_id'].isin(sequence_ids)].sort_values(['sequence_id', 'start'], ascending=[True, True]).reset_index(drop=True)
#         # get length column
#         intervals_df['length'] = intervals_df['end'] - intervals_df['start'] 
#         # get coverage matrix and drop columns of samples that are not in sample_ids_array
#         intervals_df = pd.concat([intervals_df, intervals_df.samples.str.get_dummies(sep=',').filter(sample_ids)], axis=1).drop(columns=['samples'])
#         # remove those intervals including less than two sample_ids (query sample_id columns : intervals_df[intervals_df.columns.intersection(sample_ids)])
#         intervals_df = intervals_df.loc[(intervals_df[intervals_df.columns.intersection(sample_ids)].sum(axis=1) > 1)]
#         #interval_sites_by_sample_set_idx = collections.defaultdict(list)
#         #block_sites_by_sample_set_idx = collections.defaultdict(list)
#         with tqdm(total=(len(sequence_ids) * len(sample_sets)), desc="[%] Calculating bSFSs ", ncols=100, unit_scale=True) as pbar:        
#             for seq_id in sequence_ids:        
#                 _intervals_df = intervals_df[intervals_df['sequence_id'] == seq_id]
#                 # To Do: put a check in so that it only progresses if there ARE intervals in the BED file for that sequence_id ... 
#                 seq_id_key = "%s/pos" % seq_id
#                 _pos = np.array([])
#                 if seq_id_key in self.data: 
#                     _pos = self.data["%s/pos" % seq_id] # zarr.core.array
#                     sa_genotype_array = allel.GenotypeArray(self.data["%s/gt" % seq_id])
#                 for sample_set_idx, sample_set in enumerate(sample_sets):
#                     if sample_set_idx in self.data.attrs['idx_cartesian_sample_sets']:
#                         try:
#                             sample_set_intervals_df = _intervals_df[_intervals_df[sample_set].all(axis='columns')]
#                         except KeyError:
#                             sys.exit("[X] Sample set %r not found in BED file" % sample_set)
#                         # Cut blocks based on intervals and block-algoritm parameters
#                         block_sites = cut_blocks(
#                             np.array(sample_set_intervals_df.start), 
#                             np.array(sample_set_intervals_df.end), 
#                             parameterObj.block_length, 
#                             parameterObj.block_span, 
#                             parameterObj.block_gap_run
#                             )
#                         block_starts = np.array(block_sites[:,0])           # has to be np.array() !
#                         block_ends = np.array(block_sites[:,-1] + 1)        # has to be np.array() !
#                         interval_space = sample_set_intervals_df['length'].sum()
#                         block_space = np.sum(((block_sites[:,-1] - block_sites[:,0]) + 1)) # block_space == span !
#                         if debug:
#                             print("#", seq_id, sample_set_idx, sample_set)
#                             print("# Block_sites 1", block_sites.shape)
#                             print(block_sites)
#                         # Variation
#                         if seq_id_key in self.data: 
#                             # positions inside block_sites that are variant     
#                             idx_block_sites_in_pos = np.isin(block_sites, _pos, assume_unique=True) # will crash if non-unique pos
#                             # variant positions in _pos that are within block_sites 
#                             idx_pos_in_block_sites = np.isin(_pos, block_sites, assume_unique=True) # will crash if non-unique pos
#                             # make room for reading in array
#                             sample_set_vcf_idxs = [self.data.attrs['sample_ids_to_vcf_idx'][sample_id] for sample_id in sample_set] # list of indices
#                             sa_sample_set_genotype_array = sa_genotype_array.subset(idx_pos_in_block_sites, sample_set_vcf_idxs)
#                             block_sites = genotype_to_mutype_array(sa_sample_set_genotype_array, idx_block_sites_in_pos, block_sites, debug)
#                         else:
#                             block_sites[:] = 2 # if no variants, all invariant
#                         multiallelic, missing, monomorphic, variation = block_sites_to_variation_arrays(block_sites, self.data.attrs['mutypes_count'])
#                         pi_1, pi_2, d_xy, f_st, fgv = calculate_popgen_from_array(variation, (self.data.attrs['block_length'] * variation.shape[0])) 
#                         if debug:
#                             print("# Block_sites 2")
#                             print(block_sites)
#                             print('# Variation: 0=HetB, 1=HetA, 2=HetAB, 3=Fixed')
#                             print(variation)
#                             print("[+] Pi_%s = %s; Pi_%s = %s; D_xy = %s; F_st = %s; FGV = %s" % (self.data.attrs['pop_ids'][0], pi_1, self.data.attrs['pop_ids'][1], pi_2, d_xy, f_st, fgv)) 
#                         # To Do: put explicit datatypes when saving stuff in ZARR (it might otherwise estimate it wrongly!)
#                         self.data.create_dataset("%s/%s/interval_sites" % (seq_id, sample_set_idx), data=np.array([interval_space]), overwrite=True)
#                         self.data.create_dataset("%s/%s/block_sites" % (seq_id, sample_set_idx), data=np.array([block_space]), overwrite=True)
#                         self.data.create_dataset("%s/%s/blocks/starts" % (seq_id, sample_set_idx), data=block_starts, overwrite=True)
#                         self.data.create_dataset("%s/%s/blocks/ends" % (seq_id, sample_set_idx), data=block_ends, overwrite=True)
#                         self.data.create_dataset("%s/%s/blocks/variation" % (seq_id, sample_set_idx), data=variation, overwrite=True)
#                         self.data.create_dataset("%s/%s/blocks/multiallelic" % (seq_id, sample_set_idx), data=multiallelic.flatten(), overwrite=True)
#                         self.data.create_dataset("%s/%s/blocks/missing" % (seq_id, sample_set_idx), data=missing.flatten(), overwrite=True)
#                     pbar.update(1)
#         self.add_stage(parameterObj)

#     def make_windows(self, parameterObj):
#         block_status = self.has_data('blocks')
#         if not block_status:
#             sys.exit("[X] No blocks found. Please make blocks first.")
#         if self.has_data('windows'):
#             if not parameterObj.overwrite:
#                 sys.exit("[X] Store %r already contains windows.\n[X] These windows => %r\n[X] Please specify '--force' to overwrite." % (self.path, self.get_stage_cmd('windows')))
#             print('[-] Store %r already contains windows. But these will be overwritten...' % (self.path))
#         sample_sets_idxs = self.data.attrs['idx_cartesian_sample_sets']
#         with tqdm(total=(len(self.data.attrs['sequence_ids']) * len(sample_sets_idxs)), desc="[%] Making windows ", ncols=100, unit_scale=True) as pbar: 
#             for seq_id in self.data.attrs['sequence_ids']: 
#                 mutypes, sample_set_covs, starts, ends = [], [], [], []
#                 for sample_set_idx in sample_sets_idxs:
#                     # first determine valid block mask
#                     missing = np.array(self.data["%s/%s/blocks/missing" % (seq_id, sample_set_idx)])
#                     multiallelic = np.array(self.data["%s/%s/blocks/multiallelic" % (seq_id, sample_set_idx)])
#                     valid = np.less_equal(missing, self.data.attrs['block_max_missing']) & np.less_equal(multiallelic, self.data.attrs['block_max_multiallelic'])
#                     mutype = np.array(self.data["%s/%s/blocks/variation" % (seq_id, sample_set_idx)])[valid]
#                     mutypes.append(mutype)
#                     start = np.array(self.data["%s/%s/blocks/starts" % (seq_id, sample_set_idx)])[valid]
#                     starts.append(start)
#                     end = np.array(self.data["%s/%s/blocks/ends" % (seq_id, sample_set_idx)])[valid]
#                     ends.append(end)
#                     sample_set_covs.append(np.full_like(end, sample_set_idx)) 
#                     pbar.update()
#                 mutype_array = np.concatenate(mutypes, axis=0)
#                 #print('mutype_array', mutype_array.shape, mutype_array[:])
#                 start_array = np.concatenate(starts[:], axis=0)
#                 #print('starts', starts)
#                 end_array = np.concatenate(ends, axis=0)
#                 #print('end_array', end_array.shape, end_array[:])
#                 sample_set_array = np.concatenate(sample_set_covs, axis=0)
#                 window_variation, window_starts, window_ends, window_pos_mean, window_pos_median = cut_windows(mutype_array, sample_sets_idxs, start_array, end_array, sample_set_array, num_blocks=parameterObj.window_size, num_steps=parameterObj.window_step)
#                 self.data.create_dataset("%s/windows/variation" % seq_id, data=window_variation, overwrite=True)
#                 self.data.create_dataset("%s/windows/starts" % seq_id, data=window_starts, overwrite=True)
#                 self.data.create_dataset("%s/windows/ends" % seq_id, data=window_ends, overwrite=True)
#                 self.data.create_dataset("%s/windows/pos_mean" % seq_id, data=window_pos_mean, overwrite=True)
#                 self.data.create_dataset("%s/windows/pos_median" % seq_id, data=window_pos_median, overwrite=True)
#         self.add_stage(parameterObj)
#         self.data.attrs['window_size'] = parameterObj.window_size
#         self.data.attrs['window_step'] = parameterObj.window_step

#     def dump_blocks(self, parameterObj):
#         # applies the missing & multiallelic thresholds
#         sample_sets_idxs = self.data.attrs['idx_cartesian_sample_sets']
#         data_by_key_by_sample_set_idx = collections.defaultdict(lambda: collections.defaultdict(list))
#         variation_global = []
#         with tqdm(total=(len(self.data.attrs['sequence_ids']) * len(sample_sets_idxs)), desc="[%] Writing bSFSs ", ncols=100, unit_scale=True) as pbar: 
#             for seq_id in self.data.attrs['sequence_ids']: 
#                 for sample_set_idx in sample_sets_idxs:
#                     # aggregate interval/block sites as determined from BED file and blocking algorithm
#                     data_by_key_by_sample_set_idx[sample_set_idx]['interval_sites'].append(np.array(self.data["%s/%s/interval_sites" % (seq_id, sample_set_idx)]))
#                     data_by_key_by_sample_set_idx[sample_set_idx]['block_sites'].append(np.array(self.data["%s/%s/block_sites" % (seq_id, sample_set_idx)]))
#                     # missing & multiallelic thresholds determine validity of blocks ...
#                     missing = np.array(self.data["%s/%s/blocks/missing" % (seq_id, sample_set_idx)])
#                     multiallelic = np.array(self.data["%s/%s/blocks/multiallelic" % (seq_id, sample_set_idx)])
#                     valid = np.less_equal(missing, parameterObj.block_max_missing) & np.less_equal(multiallelic, parameterObj.block_max_multiallelic)
#                     data_by_key_by_sample_set_idx[sample_set_idx]['block_sites_valid'].append(np.array([valid[valid == True].shape[0] * self.data.attrs['block_length']]))
#                     # aggregate variation/location data of valid blocks 
#                     data_by_key_by_sample_set_idx[sample_set_idx]['missing'].append(missing[valid])
#                     data_by_key_by_sample_set_idx[sample_set_idx]['multiallelic'].append(multiallelic[valid])
#                     variation = np.array(self.data["%s/%s/blocks/variation" % (seq_id, sample_set_idx)])[valid]
#                     data_by_key_by_sample_set_idx[sample_set_idx]['variation'].append(variation)
#                     variation_global.append(variation)
#                     data_by_key_by_sample_set_idx[sample_set_idx]['starts'].append(np.array(self.data["%s/%s/blocks/starts" % (seq_id, sample_set_idx)])[valid])
#                     data_by_key_by_sample_set_idx[sample_set_idx]['ends'].append(np.array(self.data["%s/%s/blocks/ends" % (seq_id, sample_set_idx)])[valid])
#                     #sample_set.append(np.full_like(end[valid], sample_set_idx)) 
#                     pbar.update()
#         variation_global_array = np.concatenate(variation_global, axis=0)
#         # popgen
#         variation_global = []
#         metrics_rows = []
#         # is order (pi_1, pi_2, d_xy, f_st, fgv) correct?
#         for sample_set_idx in data_by_key_by_sample_set_idx:
#             sample_set_ids = self.data.attrs['sample_sets'][sample_set_idx]
#             #print(data_by_key_by_sample_set_idx)
#             block_sites = np.sum(np.concatenate(data_by_key_by_sample_set_idx[sample_set_idx]['block_sites'], axis=0))
#             interval_sites = np.sum(np.concatenate(data_by_key_by_sample_set_idx[sample_set_idx]['interval_sites'], axis=0))
#             block_sites_valid = np.sum(np.concatenate(data_by_key_by_sample_set_idx[sample_set_idx]['block_sites_valid'], axis=0))
#             variation_array = np.concatenate(data_by_key_by_sample_set_idx[sample_set_idx]['variation'], axis=0)
#             missing_count = np.sum(np.concatenate(data_by_key_by_sample_set_idx[sample_set_idx]['missing'], axis=0))
#             multiallelic_count = np.sum(np.concatenate(data_by_key_by_sample_set_idx[sample_set_idx]['multiallelic'], axis=0))
#             hetB_count, hetA_count, hetAB_count, fixed_count = np.sum(variation_array, axis=0)
#             #pi_1, pi_2, d_xy, f_st, fgv = calculate_popgen_from_array(variation_array, (self.data.attrs['block_length'] * variation_array.shape[0]))    
#             pi_1, pi_2, d_xy, f_st, fgv = calculate_popgen_from_array(variation_array, block_sites_valid)    
#             metrics_rows.append([
#                 sample_set_ids[0], 
#                 sample_set_ids[1],     
#                 block_sites,
#                 interval_sites,
#                 block_sites_valid,
#                 np.divide(block_sites_valid, self.data.attrs['block_length']),
#                 fgv,
#                 missing_count,
#                 multiallelic_count,
#                 hetA_count, 
#                 hetB_count, 
#                 hetAB_count, 
#                 fixed_count,
#                 pi_1, 
#                 pi_2, 
#                 d_xy, 
#                 f_st
#                 ])
#         # output metrics 
#         header = [
#             self.data.attrs['pop_ids'][0], 
#             self.data.attrs['pop_ids'][1], 
#             'block_sites', 
#             'interval_sites', 
#             'block_sites_valid', 
#             'blocks', 
#             'fgv', 
#             'missing', 
#             'multiallelic', 
#             'hetA', 
#             'hetB', 
#             'hetAB', 
#             'fixed', 
#             'piA', 
#             'piB', 
#             'dxy', 
#             'fst'
#             ]
#         pd.DataFrame(data=metrics_rows, columns=header, dtype='int64').to_hdf("%s.block_stats.h5" % self.prefix, 'bsfs', format='table')

#         pi_1, pi_2, d_xy, f_st, fgv = calculate_popgen_from_array(variation_global_array, (self.data.attrs['block_length'] * variation_global_array.shape[0]))
#         print("[+] Pi_%s = %s; Pi_%s = %s; D_xy = %s; F_st = %s; FGVs = %s / %s blocks (%s)" % (self.data.attrs['pop_ids'][0], pi_1, self.data.attrs['pop_ids'][1], pi_2, d_xy, f_st, fgv, variation_global_array.shape[0], format_percentage(fgv / variation_global_array.shape[0]))) 
        
#         # mutuple barchart
#         mutypes, counts = np.unique(variation_global_array, return_counts=True, axis=0)
#         mutype_counter = collections.Counter({tuple(i):j for i,j in zip(mutypes, counts)})
#         plot_mutuple_barchart('%s.mutuple_barchart.png' % self.prefix, mutype_counter)

#         # mutuple tally
#         bsfs = np.concatenate([counts[:, np.newaxis], mutypes], axis =-1)
#         header = ['count'] + [x+1 for x in range(self.data.attrs['mutypes_count'])]
#         pd.DataFrame(data=bsfs, columns=header, dtype='int64').to_hdf("%s.blocks.h5" % self.prefix, 'tally', format='table')

#         # block coordinates (BED format)
#         #header = ['block_id', 'start', 'end', 'sample_set', 'multiallelic', 'missing']
#         #pd.DataFrame(data=bsfs, columns=header, dtype='int64').to_hdf("%s.blocks.h5" % self.prefix, 'bed', format='table')

#     def dump_windows(self, parameterObj):
#         window_info_rows = []
#         window_mutuple_tally = []
#         for sequence_id in tqdm(self.data.attrs['sequence_ids'], total=len(self.data.attrs['sequence_ids']), desc="[%] Generating output ", ncols=100):
#             variations = self.data["%s/windows/variation" % sequence_id]
#             #print(self.data["%s/windows/starts" % sequence_id][:])
#             #print(self.data["%s/windows/pos_mean" % sequence_id][:])
#             #print(self.data["%s/windows/pos_median" % sequence_id][:])
#             window_ids = np.array(["_".join([sequence_id, _start, _end]) for (_start, _end) in zip(
#                 np.array(self.data["%s/windows/starts" % sequence_id]).astype(str), 
#                 np.array(self.data["%s/windows/ends" % sequence_id]).astype(str))])
#             #window_ids = self.data["%s/windows/window_id" % sequence_id]
#             midpoint_means = self.data["%s/windows/pos_mean" % sequence_id]
#             midpoint_medians = self.data["%s/windows/pos_median" % sequence_id]
#             for window_id, variation, midpoint_mean, midpoint_median in zip(window_ids, variations, midpoint_means, midpoint_medians):
#                 pi_1, pi_2, d_xy, f_st, fgv = calculate_popgen_from_array(variation, (self.data.attrs['block_length'] * variation.shape[0]))
#                 window_info_rows.append([window_id, sequence_id, midpoint_mean, midpoint_median, pi_1, pi_2, d_xy, f_st, fgv/variation.shape[0]])
#                 # mutuple barchart
#                 mutypes, counts = np.unique(variation, return_counts=True, axis=0)
#                 tally = np.concatenate([counts[:, np.newaxis], mutypes], axis =-1)
#                 windows = np.array([window_id] * tally.shape[0])
#                 window_mutuple_tally.append(np.concatenate([windows[:, np.newaxis], tally], axis =-1))
        
#         window_bsfs_cols = ['window_id', 'count'] + [x+1 for x in range(self.data.attrs['mutypes_count'])]
#         window_bsfs_df = pd.DataFrame(np.vstack(window_mutuple_tally), columns=window_bsfs_cols)
#         print("[+] Made %s windows" % window_bsfs_df['window_id'].nunique()) 
#         window_bsfs_f = "%s.window_bsfs.tsv" % self.prefix
#         window_bsfs_df.to_csv(window_bsfs_f, sep='\t', index=False)
#         print("[>] Created: %r" % str(window_bsfs_f))

#         window_info_cols = ['window_id', 'sequence_id', 'midpoint_mean', 'midpoint_median', 'pi_%s' % self.data.attrs['pop_ids'][0], 'pi_%s' % self.data.attrs['pop_ids'][1], 'd_xy', 'f_st', 'fgv']
#         window_info_df = pd.DataFrame(window_info_rows, columns=window_info_cols)
#         window_info_f = "%s.window_info.tsv" % self.prefix
#         window_info_df.to_csv(window_info_f, sep='\t', index=False)
#         print("[>] Created: %r" % str(window_info_f))        
#         self.plot_fst_genome_scan(window_info_df)
#         self.plot_pi_genome_scan(window_info_df)
#         #     plot_pi_scatter(window_df, '%s.pi_scatter.png' % parameterObj.dataset)
    
#     def plot_pi_genome_scan(self, window_df):
#         offset_by_sequence_id = {}
#         offset = 0
#         x_boundaries = []
#         for sequence_id, sequence_length in zip(self.data.attrs['sequence_ids'], self.data.attrs['sequence_length']):
#             offset_by_sequence_id[sequence_id] = offset
#             x_boundaries.append(offset)
#             offset += sequence_length
#         x_boundaries.append(offset)
#         #print([(sequenceObj.id, sequenceObj.length) for sequenceObj in sequenceObjs])
#         #print(x_boundaries)
#         fig = plt.figure(figsize=(18,4), dpi=200, frameon=True)
#         #connecting dots
#         ax = fig.add_subplot(111)  
#         window_df['rel_pos'] = window_df['midpoint_median'] + window_df['sequence_id'].map(offset_by_sequence_id)
#         window_df.sort_values(['rel_pos'], inplace=True)
#         #print(window_df)
#         pi_A_key = list(window_df.columns)[4]
#         pi_B_key = list(window_df.columns)[5]
#         ax.plot(window_df['rel_pos'], window_df[pi_A_key], color='orange', alpha=0.8, linestyle='-', linewidth=1, label=pi_A_key.replace('pi_', ''))
#         ax.plot(window_df['rel_pos'], window_df[pi_B_key], color='dodgerblue', alpha=0.8, linestyle='-', linewidth=1, label=pi_B_key.replace('pi_', ''))
#         y_lim = (min(window_df[pi_A_key].min(), window_df[pi_B_key].min()), max(window_df[pi_A_key].max(), window_df[pi_B_key].max()))
#         ax.vlines(x_boundaries, y_lim[0], y_lim[1], colors=['lightgrey'], linestyles='dashed', linewidth=1)
#         ax.set_ylim(y_lim)
#         ax.spines['right'].set_visible(False)
#         ax.spines['top'].set_visible(False)
#         ax.legend(numpoints=1)
#         plt.ylabel('Pi')
#         plt.xlabel("Genome coordinate")
#         out_f = '%s.pi_genome_scan.png' % self.prefix
#         plt.tight_layout()
#         fig.savefig(out_f, format="png")
#         print("[>] Created: %r" % str(out_f))
#         plt.close(fig)

#     def plot_fst_genome_scan(self, window_df):
#         offset_by_sequence_id = {}
#         offset = 0
#         x_boundaries = []
#         for sequence_id, sequence_length in zip(self.data.attrs['sequence_ids'], self.data.attrs['sequence_length']):
#             offset_by_sequence_id[sequence_id] = offset
#             x_boundaries.append(offset)
#             offset += sequence_length
#         x_boundaries.append(offset)
#         fig = plt.figure(figsize=(18,4), dpi=200, frameon=True)
#         #connecting dots
#         ax = fig.add_subplot(111)  
#         y_lim = (0.0, 1.0)
#         window_df['rel_pos'] = window_df['midpoint_median'] + window_df['sequence_id'].map(offset_by_sequence_id)
#         window_df.sort_values(['rel_pos'], inplace=True)
#         ax.plot(window_df['rel_pos'], window_df['f_st'], color='lightgrey', alpha=0.8, linestyle='-', linewidth=1)
#         scatter = ax.scatter(window_df['rel_pos'], window_df['f_st'], c=window_df['d_xy'], alpha=1.0, cmap='PiYG_r', edgecolors='white', marker='o', s=40, linewidth=0.2)
#         cbar = fig.colorbar(scatter, ax=ax)
#         cbar.ax.set_title('D_xy')
#         ax.vlines(x_boundaries, 0.0, 1.0, colors=['lightgrey'], linestyles='dashed', linewidth=1)
#         ax.set_ylim(y_lim)
#         ax.spines['right'].set_visible(False)
#         ax.spines['top'].set_visible(False)
#         plt.ylabel('F_st')
#         plt.xlabel("Genome coordinate")
#         ax.autoscale_view(tight=None, scalex=True, scaley=True)
#         out_f = '%s.fst_genome_scan.png' % self.prefix
#         fig.savefig(out_f, format="png")
#         plt.close(fig)
#         print("[>] Created: %r" % str(out_f))

    #def dump_blocks(self, parameterObj, cartesian_only=True):
    #    meta = self.data['seqs'].attrs
    #    sample_set_idxs = [idx for (idx, is_cartesian) in enumerate(meta['sample_sets_inter']) if is_cartesian] if cartesian_only else range(len(meta['sample_sets']))
    #    variation_global = []
    #    with tqdm(total=(len(meta['seq_names']) * len(sample_set_idxs)), desc="[%] Writing bSFSs ", ncols=100, unit_scale=True) as pbar: 
    #        for seq_name in meta['seq_names']: 
    #            for sample_set_idx in sample_set_idxs:
    #                variation_key = 'seqs/%s/blocks/%s/variation' % (seq_name, sample_set_idx)
    #                variation_global.append(np.array(self.data[variation_key]))#[valid]
    #                pbar.update()
    #    variation_global_array = np.concatenate(variation_global, axis=0)
    #    # popgen
    #    variation_global = []
        #metrics_rows = []
        # is order (pi_1, pi_2, d_xy, f_st, fgv) correct?
        # for sample_set_idx in data_by_key_by_sample_set_idx:
        #     sample_set_ids = self.data.attrs['sample_sets'][sample_set_idx]
        #     #print(data_by_key_by_sample_set_idx)
        #     block_sites = np.sum(np.concatenate(data_by_key_by_sample_set_idx[sample_set_idx]['block_sites'], axis=0))
        #     interval_sites = np.sum(np.concatenate(data_by_key_by_sample_set_idx[sample_set_idx]['interval_sites'], axis=0))
        #     block_sites_valid = np.sum(np.concatenate(data_by_key_by_sample_set_idx[sample_set_idx]['block_sites_valid'], axis=0))
        #     variation_array = np.concatenate(data_by_key_by_sample_set_idx[sample_set_idx]['variation'], axis=0)
        #     missing_count = np.sum(np.concatenate(data_by_key_by_sample_set_idx[sample_set_idx]['missing'], axis=0))
        #     multiallelic_count = np.sum(np.concatenate(data_by_key_by_sample_set_idx[sample_set_idx]['multiallelic'], axis=0))
        #     hetB_count, hetA_count, hetAB_count, fixed_count = np.sum(variation_array, axis=0)
        #     #pi_1, pi_2, d_xy, f_st, fgv = calculate_popgen_from_array(variation_array, (self.data.attrs['block_length'] * variation_array.shape[0]))    
        #     pi_1, pi_2, d_xy, f_st, fgv = calculate_popgen_from_array(variation_array, block_sites_valid)    
        #     metrics_rows.append([
        #         sample_set_ids[0], 
        #         sample_set_ids[1],     
        #         block_sites,
        #         interval_sites,
        #         block_sites_valid,
        #         np.divide(block_sites_valid, self.data.attrs['block_length']),
        #         fgv,
        #         missing_count,
        #         multiallelic_count,
        #         hetA_count, 
        #         hetB_count, 
        #         hetAB_count, 
        #         fixed_count,
        #         pi_1, 
        #         pi_2, 
        #         d_xy, 
        #         f_st
        #         ])
        # # output metrics 
        # header = [
        #     self.data.attrs['pop_ids'][0], 
        #     self.data.attrs['pop_ids'][1], 
        #     'block_sites', 
        #     'interval_sites', 
        #     'block_sites_valid', 
        #     'blocks', 
        #     'fgv', 
        #     'missing', 
        #     'multiallelic', 
        #     'hetA', 
        #     'hetB', 
        #     'hetAB', 
        #     'fixed', 
        #     'piA', 
        #     'piB', 
        #     'dxy', 
        #     'fst'
        #     ]
        # pd.DataFrame(data=metrics_rows, columns=header, dtype='int64').to_hdf("%s.block_stats.h5" % self.prefix, 'bsfs', format='table')

        #pi_1, pi_2, d_xy, f_st, fgv = calculate_popgen_from_array(variation_global_array, (self.data.attrs['block_length'] * variation_global_array.shape[0]))
        #print("[+] Pi_%s = %s; Pi_%s = %s; D_xy = %s; F_st = %s; FGVs = %s / %s blocks (%s)" % (self.data.attrs['pop_ids'][0], pi_1, self.data.attrs['pop_ids'][1], pi_2, d_xy, f_st, fgv, variation_global_array.shape[0], format_percentage(fgv / variation_global_array.shape[0]))) 

