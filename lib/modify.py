from timeit import default_timer as timer
from lib.functions import plot_mutuple_barchart, check_path, check_file, check_prefix, memory_usage_psutil, format_count, format_bases, format_percentage, create_hdf5_store
from lib.classes import CoordinateTransformObj, EntityCollection
import numpy as np
import pandas as pd
from tqdm import tqdm
import collections

FULL_MUTYPE_ORDER = ['hetA', 'fixed', 'hetB', 'hetAB', 'missing', 'multiallelic']
MUTYPE_ORDER = ['hetA', 'fixed', 'hetB', 'hetAB']
MUTYPE_OTHER = ['missing', 'multiallelic']
GT_ORDER = ['TOTAL', 'MISS', 'HOM', 'HET']

class BlockParameterObj(object):
    def __init__(self, args):
        self.args = args
        self.blocks_file = check_file(args.get('--blocks', None))
        self.coordinates_file = check_file(args.get('--coordinates', None))
        self.path = check_path(args.get('--prefix', None))
        self.prefix = check_prefix(args.get('--prefix', None))
        self.dataset = "%s.modified" % self.prefix if self.path is None else "%s.modified" % (self.path / self.prefix) 
        self.threads = int(args.get('--threads'))
        self.sample_file = check_file(args.get('--sample_file', None))
        self.genome_file = check_file(args.get('--genome_file', None))
        self.max_block_length = None

class VariantParameterObj(object):
    def __init__(self, args):
        self.args = args
        self.path = check_path(args.get('--prefix', None))
        self.prefix = check_prefix(args.get('--prefix', None))
        self.dataset = "%s.modified" % self.prefix if self.path is None else "%s.modified" % (self.path / self.prefix) 
        self.threads = int(args.get('--threads'))
        self.variants_file = check_file(args.get('--variants_h5', None))
        self.max_missing = int(args.get('--max_missing'))
        self.max_multiallelic = int(args.get('--max_multiallelic'))


def parse_coordinates(parameterObj):
    mapping_df = pd.read_csv(
        parameterObj.coordinates_file, 
        sep="\t", 
        usecols=[0, 1, 2, 3, 4, 5, 6], 
        names=['chrom_id', 'chrom_start', 'chrom_end', 'sequence_id', 'sequence_start', 'sequence_end', 'sequence_orientation'], 
        skiprows=1, 
        header=None, 
        dtype={ 
            'chrom_id': 'category', 
            'chrom_start': np.int, 
            'chrom_end': np.int, 
            'sequence_id': 'category', 
            'sequence_start': np.int, 
            'sequence_end': np.int, 
            'sequence_orientation': 'category' 
            })
    # convert coordinates to 0-based (BED)
    mapping_df['chrom_start'] = mapping_df['chrom_start'] - 1
    mapping_df['sequence_start'] = mapping_df['sequence_start'] - 1
    coordinateTransformObj = CoordinateTransformObj()
    for chrom_id, chrom_start, chrom_end, sequence_id, sequence_start, sequence_end, sequence_orientation in tqdm(mapping_df.values.tolist(), total=len(mapping_df.index), desc="[%] ", ncols=100):
        coordinateTransformObj.add_tuple(sequence_id, int(sequence_start), int(sequence_end), sequence_orientation, chrom_id, int(chrom_start), int(chrom_end))
    return coordinateTransformObj
    

def task_generate_entityCollection(parameterObj):
    start = timer()
    print("[#] Building entities based on samples and sequences...")
    entityCollection = EntityCollection()
    entityCollection.parse_sample_file(parameterObj)
    print("[+] Read %s samples from %s populations and generated %s pairs in %.3fs." % (\
        entityCollection.count('samples'), 
        entityCollection.count('populations'), 
        entityCollection.count('pairs'), 
        timer() - start))
    entityCollection.parse_genome_file(parameterObj)
    print("[+] Read %s sequences with total length of %s in %.3fs" % (
        entityCollection.count('sequences'), 
        format_bases(entityCollection.count('bases')), 
        timer() - start))
    return entityCollection

def task_parse_parameters(args, mode='blocks'):
    start = timer()
    if mode == 'blocks':
        parameterObj = BlockParameterObj(args)
    elif mode == 'variants':
        parameterObj = VariantParameterObj(args)
    print("[+] Read parameters in %.3fs (%.2fMB)" % (timer() - start, memory_usage_psutil()))
    return parameterObj

def task_filter_variants(parameterObj):
    print("[#] Loading variants ...")
    variants_hdf5_store = pd.HDFStore(parameterObj.variants_file)    
    variant_cols = ['block_id', 'pair_idx', 'hetA', 'fixed', 'hetB', 'hetAB', 'missing', 'multiallelic']
    variants_block_df = pd.read_hdf(variants_hdf5_store, key='blocks').reindex(columns=variant_cols)
    variants_hdf5_store.close()
    filtered_variants_block_df = variants_block_df[(
        (variants_block_df['missing'] <= parameterObj.max_missing) & 
        (variants_block_df['multiallelic'] <= parameterObj.max_multiallelic)
        )]
    print("[+] Excluded %s out of %s mutuples (%s)" % (
        (len(variants_block_df) - len(filtered_variants_block_df)),
        len(variants_block_df),
        format_percentage(1 - len(filtered_variants_block_df) / len(variants_block_df))
        ))
    
    mutuples = []
    print("[#] Analysing variants in blocks...")
    for block_id, pair_idx, hetA, fixed, hetB, hetAB, missing, multiallelic in tqdm(filtered_variants_block_df.values.tolist(), total=len(filtered_variants_block_df), desc="[%] ", ncols=100, unit_scale=True):
        mutuples.append((hetA, fixed, hetB, hetAB))

    print("[#] Creating dataframe of global mutuple tallies...")
    global_mutuple_tally_cols = ['count'] + MUTYPE_ORDER 
    global_mutuple_tally_vals = []
    global_mutype_counter = collections.Counter(mutuples)
    plot_mutuple_barchart(
            '%s.mutuple_barchart.png' % parameterObj.dataset,
            global_mutype_counter
            )
    for mutype, count in global_mutype_counter.most_common():
        global_mutuple_tally_vals.append([count] + list(mutype))
    out_f = '%s.variants.h5' % parameterObj.dataset
    filtered_variants_hdf5_store = create_hdf5_store(
        out_f=out_f, 
        path=parameterObj.path
        )
    filtered_variants_block_df.to_hdf(filtered_variants_hdf5_store, 'blocks', append=True)
    global_mutuple_tally_df = pd.DataFrame(global_mutuple_tally_vals, columns=global_mutuple_tally_cols)
    global_mutuple_tally_df.to_hdf(filtered_variants_hdf5_store, 'mutypes', append=True)
    filtered_variants_hdf5_store.close()

def task_load_blockObjs(parameterObj, entityCollection):
    print("[#] Loading blocks ...")
    entityCollection.load_blocks(parameterObj)
    block_count = format_count(entityCollection.count('blocks'))
    blocked_bases = entityCollection.count('blocks') * parameterObj.block_length
    total_bases_percentage = format_percentage(blocked_bases / entityCollection.count('bases')) 
    print("[+] Read %s blocks covering %s (%s of genome) (%.2fMB)" % (
        block_count, 
        format_bases(blocked_bases),
        total_bases_percentage,
        memory_usage_psutil()
        ))

def task_parse_coordinates(parameterObj):
    start = timer()
    print("[#] Parse coordinate transformation file ...")
    coordinateTransformObj = parse_coordinates(parameterObj)
    print("[+] Read file in %.3fs (%.2fMB)" % (timer() - start, memory_usage_psutil()))
    return coordinateTransformObj

def task_transform_coordinates(parameterObj, entityCollection, coordinateTransformObj):
    start = timer()
    print("[#] Transforming coordinates ...")
    entityCollection.transform_coordinates(parameterObj, coordinateTransformObj)
    print("[+] Transformed coordinates in %.3fs (%.2fMB)" % (timer() - start, memory_usage_psutil()))

def task_write_modify_output(parameterObj, entityCollection):
    entityCollection.generate_modify_output(parameterObj)
