from timeit import default_timer as timer 
from lib.functions import check_file, memory_usage_psutil, check_path, check_prefix, format_count, format_percentage, format_bases
from lib.classes import EntityCollection, WindowObj

from tqdm import tqdm
import more_itertools
#import collections


'''
[To Do]
- read variant json vs tsv
    
    # json (new dependency: simplejson)
        - size      : 408.3 Mb
        - reading   : 12.827509473000006 s
        - writing   : 100.27862279300001 s
    
    # tsv 
        - size      : 163.8 Mb
        - reading   : 48.378815994 s (w/ progressbar)
        - writing   : 3.968107513999996 s

    # hd5 (new dependency: pytables)
        - size      : 167.1 Mb
        - reading   : 44.776825606 (w/ progressbar)
        - writing   : 1.1913998750000019 s
        * lot of possibility for improvement by doing multidimensional dataframes
            - https://stackoverflow.com/questions/13575090/construct-pandas-dataframe-from-items-in-nested-dictionary
    import pandas as pd
    a = pd.read_hdf('/Users/dlaetsch/git/gIMble/hmel.chr18/win_test.variants.blocks.hd5')
    b = pd.read_hdf('/Users/dlaetsch/git/gIMble/hmel.chr18/var_test.variants.blocks.hd5')
    a == b

- re-use block parsing
- split blocks based on max_block_distance
- use itertools windows to make windows
    - re-use windownames
- generate output
    - window-based mutype counters (use faster format)

'''

class ParameterObj(object):
    def __init__(self, args):
        self.args = args
        self.sample_file = check_file(args.get('--sample_file', None))
        self.genome_file = check_file(args.get('--genome_file', None))
        self.blocks_file = check_file(args.get('--blocks_h5', None))
        self.variants_file = check_file(args.get('--variants_h5', None))
        self.path = check_path(args.get('--prefix', None))
        self.prefix = check_prefix(args.get('--prefix', None))
        self.dataset = self.prefix if self.path is None else (self.path / self.prefix)
        self.window_size = int(args.get('--window_size'))
        self.overlap = int(args.get('--overlap'))
        self.max_block_distance = int(args.get('--max_block_distance'))
        self.threads = int(args.get('--threads'))
        self.block_length = None
        self.max_block_length = None

def make_windows(parameterObj, entityCollection, block_id_batches):
    windowObjs = []
    for block_id_batch in tqdm(block_id_batches, total=len(block_id_batches), desc="[%] ", ncols=100, unit_scale=True):
        for block_ids in more_itertools.windowed(block_id_batch, parameterObj.window_size, step=parameterObj.overlap):
            if not block_ids[-1] is None:
                blockObjs = [entityCollection.blockObjs[entityCollection.block_idx_by_block_id[block_id]] for block_id in block_ids]
                windowObjs.append(WindowObj(blockObjs))
    entityCollection.add_windowObjs(windowObjs)

def task_generate_entityCollection(parameterObj):
    start = timer()
    print("[#] Building entities based on samples and sequences...")
    entityCollection = EntityCollection()
    entityCollection.parse_sample_file(parameterObj)
    print("[+] Read %s samples from %s populations and generated %s pairs in %.3fs." % (\
        entityCollection.count('samples'), \
        entityCollection.count('populations'), \
        entityCollection.count('pairs'), \
        timer() - start))
    entityCollection.parse_genome_file(parameterObj)
    print("[+] Read %s sequences with total length of %s b in %.3fs" % (\
        entityCollection.count('sequences'), \
        entityCollection.count('bases'), \
        timer() - start))
    return entityCollection

def task_parse_parameters(args):
    start = timer()
    print("[#] Parsing parameters ...")
    parameterObj = ParameterObj(args)
    print("[+] Read parameters in %.3fs (%.2fMB)" % (timer() - start, memory_usage_psutil()))
    return parameterObj

def task_get_block_regions(parameterObj, entityCollection):
    print("[#] Loading blocks ...")
    block_dfs = entityCollection.load_blocks(parameterObj, purpose='windows')
    block_count = format_count(entityCollection.count('blocks'))
    blocked_bases = entityCollection.count('blocks') * parameterObj.block_length
    total_bases_percentage = format_percentage(blocked_bases / entityCollection.count('bases')) 
    print("[+] Read %s blocks covering %s (%s of genome) (%.2fMB)" % (
        block_count, 
        format_bases(blocked_bases),
        total_bases_percentage,
        memory_usage_psutil()
        ))
    return block_dfs

def task_get_variants_df(parameterObj, entityCollection):
    print("[#] Loading variants ...")
    return entityCollection.get_mutype_counters_by_block_id(parameterObj)

def task_make_windows(parameterObj, entityCollection, block_dfs):
    print("[#] Making windows ...")
    make_windows(parameterObj, entityCollection, block_dfs)

def task_write_window_output(parameterObj, entityCollection, mutype_counters_by_block_id):
    print("[#] Generating output ...")
    entityCollection.generate_window_output(parameterObj, entityCollection, mutype_counters_by_block_id)