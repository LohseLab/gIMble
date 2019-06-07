from pandas import read_csv 
from timeit import default_timer as timer
from lib.classes import BlockObj, BedObj, EntityCollection

import collections

from more_itertools import flatten
import numpy as np 
from tqdm import tqdm
#### Entities

import sys
from lib.functions import memory_usage_psutil, check_file, poolcontext, check_path, check_prefix, format_bases, format_count, format_percentage

'''
- cli simple
- setup.py (example https://github.com/tskit-dev/tszip/blob/master/setup.py)
 => entry_points={
        'console_scripts': [
            'tszip=tszip.cli:tszip_main',
            'tsunzip=tszip.cli:tsunzip_main',
'''

class ParameterObj(object):
    def __init__(self, args):
        self.args = args
        self.sample_file = check_file(args.get('--sample_file', None))
        self.genome_file = check_file(args.get('--genome_file', None))
        self.bed_file = check_file(args.get('--bed_file', None))
        self.path = check_path(args.get('--prefix', None))
        self.prefix = check_prefix(args.get('--prefix', None))
        self.dataset = self.prefix if self.path is None else (self.path / self.prefix)

        self.block_length = int(args['--block_length'])
        self.min_interval_length = int(args['--min_interval_length'])
        self.max_block_length = int(args['--max_block_span'])
        if self.max_block_length < self.block_length:
            sys.exit("[X] '--max_block_length' must be greater or equal to '--block_length'")
        self.max_interval_distance = (self.max_block_length - self.block_length) + 1 
        self.threads = int(args.get('--threads'))
        self.bed_length_total = 0

def generate_region_dfs(parameterObj, entityCollection):
    '''
    Generates a list of region_dfs 
    '''
    print("[#] Splitting bed file into chunks for downstream processing (this might take a while) ...")
    df = read_csv( 
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
    df = df[df['sequence_id'].isin(entityCollection.sequence_idx_by_sequence_id)] # filter rows based on sequence_ids
    df['length'] =  df['end'] - df['start'] # compute length column
    parameterObj.bed_length_total = int(df['length'].sum())
    genome_bases_percentage = format_percentage(parameterObj.bed_length_total / entityCollection.count('bases'))
    print("[+] Found %s BED intervals adding up to %s (%s of genome) ..." % (
        format_count(len(df)), 
        format_bases(parameterObj.bed_length_total), 
        genome_bases_percentage
        ))
    df = df[df['length'] >= parameterObj.min_interval_length] # filter intervals shorter than MIN_INTERVAL_LEN
    #df['samples_ids'] = df['samples'].apply(entityCollection.sample_string_to_sample_ids) # samples to frozenset
    df['pair_idxs'] = df['samples'].apply(entityCollection.sample_string_to_pair_idxs) # samples to frozenset pairs
    df = df.dropna() # Drop intervals that don't affect pairs
    df['distance'] = np.where((df['sequence_id'] == df['sequence_id'].shift(-1)), df['start'].shift(-1) - df['end'], parameterObj.max_interval_distance + 1) # compute distance to next interval
    # region_dfs
    region_ids = (df["distance"].shift() > float(parameterObj.max_interval_distance)).cumsum() # generate indices for splitting
    region_dfs = []
    for idx, region_df in df.groupby(region_ids):
        if region_df.length.sum() > parameterObj.block_length: # remove regions below block_length
            region_df = region_df.drop(columns=['distance', 'samples'], axis=1) # remove distance/sample_idsx columns
            region_dfs.append(region_df)
    return region_dfs

def make_blocks(parameterObj, region_dfs, entityCollection):
    results = []
    print("[+] Analysing %s BED regions using %s threads ..." % (len(region_dfs), parameterObj.threads))
    if parameterObj.threads < 2:
        for region_idx, region_df in enumerate(tqdm(region_dfs, total=len(region_dfs), desc="[%] ", ncols=100, unit_scale=True)):
            results.append(block_algorithm((region_df, region_idx, parameterObj, entityCollection)))
    else:
        # if multiple threads then arguments have to be passed to algorithm
        params = [(region_df, region_idx, parameterObj, entityCollection) for region_idx, region_df in enumerate(region_dfs)]
        with poolcontext(processes=parameterObj.threads) as pool:
            with tqdm(total=len(params), desc="[%] ", ncols=100, unit_scale=True) as pbar:
                for blockObjs in pool.imap_unordered(block_algorithm, params):
                    results.append(blockObjs)
                    pbar.update()
    entityCollection.add_blockObjs(list(flatten(results)))

def block_algorithm(params):
    region_df, region_idx, parameterObj, entityCollection = params
    contig_id = region_df['sequence_id'].iloc[0]
    block_idx = 0
    blockObjs = []
    block_id = "%s.i%s.b%s" % (contig_id, region_idx, block_idx)
    blockObj = BlockObj(block_id, parameterObj.block_length)
    bedObjs = collections.deque([BedObj(sequence_id, start, end, length, pair_idxs) for sequence_id, start, end, pair_idxs, length in region_df.values.tolist()])
    while 1:
        try:
            bedObj = bedObjs.popleft()
            bedObj_score = (len(bedObj.pair_idxs) / len(entityCollection.pairObjs)) * (min(bedObj.length, blockObj.needed) / parameterObj.block_length)
            remainder_bedObj = blockObj.add_bedObj(bedObj, parameterObj, entityCollection)
            if blockObj.score >= bedObj_score:
                if not blockObj.needed:
                    blockObjs.append(blockObj)
                    block_idx += 1
                    block_id = "%s.i%s.b%s" % (contig_id, region_idx, block_idx)
                    blockObj = BlockObj(block_id, parameterObj.block_length)
                    if remainder_bedObj:
                        bedObjs.appendleft(remainder_bedObj)
            else: # span violation
                bedObjs.appendleft(bedObj)
                blockObj = BlockObj(block_id, parameterObj.block_length)
        except IndexError:
            break
    return blockObjs

### TASKS

def task_generate_parameterObj(args):
    start = timer()
    parameterObj = ParameterObj(args)
    print("[+] Read parameters in %.3fs (%.2fMB)" % (timer() - start, memory_usage_psutil()))
    return parameterObj

def task_generate_entityCollection(parameterObj):
    start = timer()
    entityCollection = EntityCollection()
    entityCollection.parse_sample_file(parameterObj)
    print("[+] Read %s samples from %s populations and generated %s pairs in %.3fs." % (\
        entityCollection.count('samples'), \
        entityCollection.count('populations'), \
        entityCollection.count('pairs'), \
        timer() - start))
    entityCollection.parse_genome_file(parameterObj)
    print("[+] Read %s sequences with total length of %s in %.3fs" % (\
        entityCollection.count('sequences'), \
        format_bases(entityCollection.count('bases')), \
        timer() - start))
    return entityCollection

def task_generate_region_dfs(parameterObj, entityCollection):
    print("[#] Processing BED file ...")
    region_dfs = generate_region_dfs(parameterObj, entityCollection)
    print("[+] BED intervals processed (%.2fMB)" % (memory_usage_psutil()))
    return region_dfs

def task_make_blocks(parameterObj, region_dfs, entityCollection):
    print("[#] Generating blocks ...")
    make_blocks(parameterObj, region_dfs, entityCollection)
    block_count = format_count(entityCollection.count('blocks'))
    blocked_bases = entityCollection.count('blocks') * parameterObj.block_length
    bed_bases_percentage = format_percentage(blocked_bases / parameterObj.bed_length_total)
    total_bases_percentage = format_percentage(blocked_bases / entityCollection.count('bases')) 
    print("[+] Made %s blocks covering %s (%s of BED intervals, %s of genome) (%.2fMB)" % (\
        block_count, 
        format_bases(blocked_bases),
        bed_bases_percentage,
        total_bases_percentage,
        memory_usage_psutil()
        ))

def task_write_block_output(parameterObj, entityCollection):
    print("[#] Generating output ...")
    entityCollection.generate_block_output(parameterObj)
