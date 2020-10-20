#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""usage: gIMble simulate                   [-z <DIR>] [-o <DIR> | -b <INT>] -c <FILE> [-r <INT>] [-t <INT>] [-h|--help] [-l <STR>]
                                            
    Options:
        -h --help                                   show this
        -z, --zarr DIR                              Path to zarr store
        -o, --outprefix DIR                         prefix to make new zarr store
        -c, --config_file FILE                      Config file with parameters
        -b, --blocks INT                            Number of blocks per replicate
        -r, --replicates INT                        Number of replicates per parametercombo
        -t, --threads INT                           Threads [default: 1]
        -l, --label STR                             Custom name for simulation run
        
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
        self._set_or_write_config(args["--blocks"], args["--replicates"])
        self._set_recombination_rate()  
        
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
        rmap_path = self.config["simulations"]["recombination_map"]
        if os.path.isfile(rmap_path):
            #check validity of recombinatin map in _parse_recombination_map()
            rbins = self.config["simulations"]["number_bins"]
            cutoff = self.config["simulations"]["cutoff"]
            scale = self.config["simulations"]["scale"]
            window_bin, to_be_simulated = self._parse_recombination_map(rmap_path, cutoff, rbins, scale)
            self.config["parameters"]["recombination"] = to_be_simulated
        else:
            self.config["parameters"]["recombination"] = [self.config["simulations"]["recombination_rate"]]

    def _parse_recombination_map(self, v, cutoff, bins, scale):
        #load bedfile
        hapmap = pd.read_csv(v, sep='\t', 
            names=['chrom', 'chromStart', 'chromEnd', 'rec'], header=None)
        #check validity
        #self._check_recombination_map(self.zstore, hapmap)
        #from cM/Mb to rec/bp
        hapmap['rec_scaled'] = hapmap['rec']*1e-8
        return self._make_bins(hapmap['rec_scaled'], scale, cutoff, bins)

    def _make_bins(self, x, scale, cutoff=90, bins=10):
        x=np.array(x)
        clip_value = np.percentile(x,cutoff)
        clip_array = np.clip(x, a_min=None, a_max=clip_value)
        with_zero = np.count_nonzero(x==0.0)
        if scale.lower() == "log":
            start, stop = np.min(clip_array), np.max(clip_array)
            if with_zero:
                zero_index = np.where(clip_array==0.0)
                start=np.min(np.delete(clip_array, zero_index))
            start, stop = np.log10(start), np.log10(stop)
            bin_edges = np.logspace(start, stop, num=bins)
            counts, bin_edges = np.histogram(clip_array, bins=bin_edges)
        
        elif scale.lower() == "lin": 
            counts, bin_edges = np.histogram(clip_array, bins=bins)
        else:
            sys.exit("[-] Scale of recombination values to be simulated should either be LINEAR or LOG")
        bin_with_counts = counts>0
        no_counts = bins-sum(1 for v in bin_with_counts if v)
        print(f"[+] There are {no_counts} recombination bins without counts")
        to_be_simulated  = [(stop + start)/2 for start, stop, c in zip(bin_edges[:-1],bin_edges[1:], bin_with_counts) if c]
        window_bin = np.digitize(clip_array, bin_edges, right=True)
        return (window_bin, list(to_be_simulated))

    def _validate_recombination_map(self, store, df):
        starts = store.data[chrom]['windows/starts'][:]
        ends = store.data[chrom]['windows/ends'][:]
        if set(starts) != set(df['starts']):
            sys.exit("[X] Starts recombination map do not match window coordinates")
        if set(ends) != set(df['ends']):
            sys.exit("[X] Ends recombination map do not match window coordinates")

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
        
        gimbleStore.simulate(parameterObj)

        print("[*] Total runtime: %.3fs" % (timer() - start_time))
    except KeyboardInterrupt:
        print("\n[X] Interrupted by user after %s seconds!\n" % (timer() - start_time))
        exit(-1)
