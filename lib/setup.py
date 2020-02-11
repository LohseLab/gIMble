# import pathlib
# import allel
import pandas as pd
import numpy as np
# import os
# import shutil
# import zarr
# import matplotlib as mpl
# import matplotlib.pyplot as plt
# import seaborn as sns
# sns.set_style('white')
# sns.set_style('ticks')
# sns.set_context('notebook')
# import re

# import itertools
# import collections
# from tqdm import tqdm
# import scipy
# import lib.functions
# import logging

class ParameterObj(object):
    def __init__(self, args):
        print(args)
        self.module = 'setup'
        self.vcf_file = args['--vcf']
        self.bed_file = args['--bed']
        self.genome_file = args['--genome']
        self.sample_file = args['--sample'] 
        self.outprefix = args['--outprefix']

        self.block_length = int(args['--block_length'])
        self.block_span = int(args['--block_span'])
        self.block_gap_run = int(args['--block_gap_run'])

        self.pairedness = int(args['--pairedness'])

#######################################################

#def get_metricObjs(storeObj, parameterObj):
#    metricObjs, fails = [], []
#    for seq_id, seq_length in tqdm(metadataObj.yield_seq_tuple(), total=len(metadataObj.seq_ids), desc="[%] ", ncols=100):
#        metricObj = MetricObj(seq_id, seq_length)
#        metricObj.valid = (seq_id in dataset)
#        if metricObj.valid:
#            metricObj.collect_metrics(dataset, metadataObj, indels_flag)
#            metricObjs.append(metricObj)
#        else:
#            fails.append(seq_id)
#    for seq_id in fails:
#        logging.info("[-] No variation data found for Seq ID %r" % seq_id)
#    return metricObjs

# import logging

COLOURS = ['orange', 'cornflowerblue', 'lightskyblue', 'gold']
SHADES = ['black', 'darkgrey', 'lightgrey', 'white'] 

'''
import numpy as np; 
from timeit import default_timer as timer; 
starts = np.arange(0, 10**8, 100); 
ends = starts + 50; 
a = np.array((starts, ends)).T; 
logging.info('# create_ranges()'); 
t1=timer(); 
cr = create_ranges(a); 
logging.info(timer()-t1) ; 
logging.info('# list_comprehension()'); 
t2=timer(); 
lc = np.concatenate([np.arange(x[0], x[1]) for x in a]); 
logging.info(timer()-t2); 
logging.info("Same = %s" % all(cr==lc))
# create_ranges()
0.46567643700001327
# list_comprehension()
2.9621128510000005
Same = True
'''
def get_accessability_matrix_from_bed_file(bed_file, samples):
    ''' samples -> [array_like, bool]'''
    def read_bed(bedfile):
        bed_df = pd.read_csv(
        bedfile,
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
        return bed_df

def create_ranges(aranges):
    # https://stackoverflow.com/a/47126435
    #print("aranges", aranges)
    l = aranges[:,1] - aranges[:,0]
    #print("l", l)
    clens = l.cumsum()
    #print("clens", clens)
    ids = np.ones(clens[-1], dtype=int)
    #print("ids", ids)
    ids[0] = aranges[0,0]
    #print("ids[0]", ids[0])
    ids[clens[:-1]] = aranges[1:, 0] - aranges[:-1, 1]+1
    #print("ids[clens[:-1]]", ids[clens[:-1]])
    #print("ids.cumsum()", ids.cumsum())
    return ids.cumsum()

# if __name__ == "__main__":
#     '''
#     -> GFF/ZARR : https://www.biostars.org/p/335077/
#     '''
#     start_time = timer()
#     args = docopt(__doc__)
#     logging.info(args)
#     # MAIN
#     MAX_NON_SNP_DISTANCE = 25
#     WINDOWSIZE = 10000
#     GT_LABELS = ['TOTAL', 'MISS', 'HET', 'HOMALT', 'HOMREF']
#     metadataObj, metricObjs = process_zarr(args['--zarr'], args['--indels'], args['--prefix'])
#     generate_results(metadataObj, metricObjs)