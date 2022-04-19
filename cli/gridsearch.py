#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""usage: gimble gridsearch -z <FILE> -g <STR> [-t <STR>|-s <STR>] [-n <INT> -c <INT> -f] [-h|--help]
                                            
                                            
    Options:
        -h --help                                   show this
        
        -z, --zarr_file <FILE>                      Path to existing GimbleStore
        -g, --grid_label <STR>                      Label of makegrid run in GimbleStore
        -t, --tally_label <STR>                     Tally label
        -s, --sim_label <STR>                       Simulation label
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
        self.tally_label = args['--tally_label']
        self.sim_label = args['--sim_label']
        self.grid_label = args['--grid_label']
        self.overwrite = args['--overwrite']
        self.num_cores = self._get_int(args['--num_cores'])    # number of workers for independent processes
        self.chunksize = self._get_int(args['--chunksize'])    # size of chunks in first dimension of tally/grid dask array 
        
def main(params):
    try:
        start_time = timer()
        args = docopt(__doc__)
        parameterObj = GridsearchParameterObj(params, args)
        gimbleStore = lib.gimble.Store(path=parameterObj.zstore, create=False)
        gimbleStore.gridsearch(
            tally_label=parameterObj.tally_label,
            sim_label=parameterObj.sim_label,
            grid_label=parameterObj.grid_label,
            num_cores=parameterObj.num_cores,
            chunksize=parameterObj.chunksize,
            overwrite=parameterObj.overwrite,
            )
        print("[*] Total runtime was %s" % (lib.gimble.format_time(timer() - start_time)))
    except KeyboardInterrupt:
        print("\n[X] Interrupted by user after %s !\n" % (lib.gimble.format_time(timer() - start_time)))
        exit(-1)