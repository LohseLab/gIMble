#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""usage: gIMble simulate                   [-z <DIR>] [-o <DIR> | -b <INT>] -c <FILE> [-r <INT>] [-t <INT>] [-h|--help] [-l <STR>] [-g]
                                            
    Options:
        -h --help                                   show this
        -z, --zarr DIR                              Path to zarr store
        -o, --outprefix DIR                         prefix to make new zarr store
        -c, --config_file FILE                      Config file with parameters
        -b, --blocks INT                            Number of blocks per replicate
        -r, --replicates INT                        Number of replicates per parametercombo
        -t, --threads INT                           Threads [default: 1]
        -l, --label STR                             Custom name for simulation run
        -g, --grid
        
"""
import pathlib
import collections
from timeit import default_timer as timer
from docopt import docopt
import sys, os
import lib.gimble
import lib.simulate
import numpy as np
import zarr
import pandas as pd
import itertools

"""
test command:
gIMble simulate -m 
./gIMble model -s A,B -p 2 -n 1,1 -m 'A>B' -j 'A,B' -o output/test
./gIMble model -s A,B -p 2 -n 2,1 -j 'A,B' -o output/test

./gIMble simulate -m output/s_A_B.p2.n_1_1.m_AtoB.j_A_B.model.tsv
./gIMble simulate -m output/test.model.tsv

./gIMble simulate -m output/s_A_B.p2.n_1_1.m_AtoB.j_A_B.model.tsv -c /Users/s1854903/git/gIMble/output/s_A_B.p2.n_1_1.m_AtoB.j_A_B.model.config.yaml -z /Users/s1854903/git/gIMble/test.z
./gIMble simulate -z output/sims.z -m output/test.model.tsv -c output/test.model.config.yaml
./gIMble simulate -z output/test.z -w output/windows.tsv 
"""


class SimulateParameterObj(lib.gimble.ParameterObj):
    """Sanitises command line arguments and stores parameters."""

    def __init__(self, params, args):
        super().__init__(params)
        self.config_file = self._get_path(args["--config_file"])
        self.zstore = self._get_path(args["--zarr"])
        self.prefix = self._get_prefix(args["--outprefix"])
        self.threads = self._get_int(args["--threads"])
        self.label = args["--label"]
        self.sim_grid = args["--grid"]
        self._set_or_write_config(args["--blocks"], args["--replicates"])
        #self._set_recombination_rate()  
        
    def _set_or_write_config(self, blocks, replicates):
        #in case no config file is provided
        if self.config_file is None:
            print("[-] No config file found.")
            sys.exit("[X] Use gimble model to generate a config file.")
        else:
            print("[+] Validating config file.")
            self._parse_config(self.config_file)
            if blocks:
                blocks = self._get_int(blocks)
                self.config['simulations']['blocks'] = blocks
            if replicates:
                replicates = self._get_int(replicates)
                self.config['simulations']['replicates'] = replicates

    def _set_recombination_rate(self):      
        self.recombination_map = None
        rmap_path = self.config["simulations"]["recombination_map"]
        rec = self.config["simulations"]["recombination_rate"]
        rec=0.0 if rec==None else rec
        if os.path.isfile(rmap_path):
            if not self.sim_grid:
                print("[-] A recombination map can only be used with the flag --sim_grid.")
                print(f"[-] Simulate will continue using r={rec}") 
                self.config["parameters"]["recombination"] = [rec,]    
            else:
                rbins = self.config["simulations"]["number_bins"]
                cutoff = self.config["simulations"]["cutoff"]
                scale = self.config["simulations"]["scale"]
                self.config["parameters"]["recombination"] = [None,]
                self.recombination_map = self._parse_recombination_map(rmap_path, cutoff, rbins, scale)
        else:
            self.config["parameters"]["recombination"] = [rec,]

    def _parse_recombination_map(self, path, cutoff, bins, scale):
        #load bedfile
        hapmap = pd.read_csv(path, sep='\t', 
            names=['sequence', 'start', 'end', 'rec'], header=None)
        #check validity: is number of columns correct?
        #from cM/Mb to rec/bp
        hapmap['rec_scaled'] = hapmap['rec']*1e-8
        return self._make_bins(hapmap, scale, cutoff, bins)

    def _make_bins(self, df, scale, cutoff=90, bins=10):
        clip_value = np.percentile(df['rec_scaled'], cutoff)
        df['rec_clipped'] = df['rec_scaled'].clip(lower=None, upper=clip_value)
        df['rec_clipped'].replace(0,np.nan,inplace=True)
        start, stop =  df['rec_clipped'].min(), df['rec_clipped'].max()
        #determine bins
        if scale.lower() == "log":
            start, stop = np.log10(start), np.log10(stop)
            bin_edges = np.logspace(start, stop, num=bins+1)
        elif scale.lower() == "lin": 
            bin_edges = np.linspace(start, stop, num=bins+1)
        else:
            sys.exit("[X] Scale of recombination values to be simulated should either be LINEAR or LOG")
        to_be_simulated  = [(bstop + bstart)/2 for bstart, bstop in zip(bin_edges[:-1],bin_edges[1:])]     
        df['rec_bins'] = pd.cut(df['rec_clipped'], bins, labels=to_be_simulated).astype(float)
        df['rec_bins'].replace(np.nan, 0.0, inplace=True)
        return df

    def _validate_recombination_map(self, store, df):
        meta_seqs = store._get_meta('seqs')
        meta_windows = store._get_meta('windows')
        MAX_SEQNAME_LENGTH = max([len(seq_name) for seq_name in meta_seqs['seq_names']])
        sequences = np.zeros(meta_windows['count'], dtype='<U%s' % MAX_SEQNAME_LENGTH)
        starts = np.zeros(meta_windows['count'], dtype=np.int64)
        ends = np.zeros(meta_windows['count'], dtype=np.int64)
        df_store = pd.DataFrame({'sequence':sequences, 'start':starts, 'end':ends})

        df_to_test = df[['sequence', 'start', 'end']]
        df_merge = df_to_test.merge(df_store, how='outer', on=['sequence', 'start', 'end'])
        if df_store.shape != df_merge.shape:
            sys.exit("[X] Recombination map coordinates do not match window coordinates. Use query to get window coordinates.")
        

def main(params):
    try:
        start_time = timer()
        args = docopt(__doc__)
        parameterObj = SimulateParameterObj(params, args)
        if parameterObj.zstore:
            gimbleStore = lib.gimble.Store(path=parameterObj.zstore)
        elif parameterObj.prefix:
            gimbleStore = lib.gimble.Store(prefix=parameterObj.prefix, create=True)
        else:
            sys.exit("[X] No config and no prefix specified. Should have been caught.")
        #perform recmap checks
        gimbleStore.simulate(parameterObj)

        print("[*] Total runtime: %.3fs" % (timer() - start_time))
    except KeyboardInterrupt:
        print("\n[X] Interrupted by user after %s seconds!\n" % (timer() - start_time))
        exit(-1)
