#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""usage: gimble gridsearch -z <FILE> -g <STR> (-t <STR>|-s <STR> [-w]) [-n <INT> -c <INT> -f] [-h|--help]
                                            
                                            
    Options:
        -h --help                                   show this
        
        -z, --zarr_file <FILE>                      Path to existing GimbleStore
        -g, --grid_key <STR>                        Makegrid key run in GimbleStore
        -t, --tally_key <STR>                       Tally key
        -s, --sim_key <STR>                         Simulation label
        -w, --windowsum                             Sum simulation windows [default: False]
        -n, --num_cores <INT>                       Number of cores [default: 1]
        -c, --chunksize <INT>                       Size of chunks to use in parallelisation 
                                                        (greater chunksize => greater RAM requirements)
                                                        [default: 500]
        -f, --overwrite                             Overwrite results in GimbleStore

"""
from timeit import default_timer as timer
from docopt import docopt
import lib.gimble
import lib.math

class GridsearchParameterObj(lib.gimble.ParameterObj):
    def __init__(self, params, args):
        super().__init__(params)
        self.zstore = self._get_path(args['--zarr_file'])
        self.tally_key = args['--tally_key']
        self.sim_key = args['--sim_key']
        self.windowsum = args['--windowsum']
        self.grid_key = args['--grid_key']
        self.overwrite = args['--overwrite']
        self.num_cores = self._get_int(args['--num_cores'])    # number of cores for independent processes
        self.chunksize = self._get_int(args['--chunksize'])    # size of chunks in first dimension of tally/grid dask array 


def main(params):
    try:
        start_time = timer()
        args = docopt(__doc__)
        parameterObj = GridsearchParameterObj(params, args)
        gimbleStore = lib.gimble.Store(path=parameterObj.zstore, create=False)
        gimbleStore.gridsearch(
            tally_key=parameterObj.tally_key,
            sim_key=parameterObj.sim_key,
            grid_key=parameterObj.grid_key,
            windowsum=parameterObj.windowsum,
            num_cores=parameterObj.num_cores,
            chunksize=parameterObj.chunksize,
            overwrite=parameterObj.overwrite,
            )
        print("[*] Total runtime was %s" % (lib.gimble.format_time(timer() - start_time)))
    except KeyboardInterrupt:
        print("\n[X] Interrupted by user after %s !\n" % (lib.gimble.format_time(timer() - start_time)))
        exit(-1)