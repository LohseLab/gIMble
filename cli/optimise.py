#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
usage: gimble optimize                  [-z FILE] -c FILE -m FILE (-b [-n INT]|-w| --simID <STR>) [--track] 
                                            [-x FLOAT -i INT] [-f FLOAT] [-p] [--inner_pool <INT> --outer_pool <INT>]
                                            [-h|--help]
                                            
                                            
    Options:
        -h --help                                   show this
        -z, --zarr_file FILE                        ZARR datastore
        -c, --config_file FILE                      INI config file
        -m, --model_file FILE                       gimble model TSV
        -b, --blocks                                Optimize based on blocks 
        -w, --windows                               Optimize based on windows (might take very long)
        --simID STR                                 Provide name of simulation run to optimize
        --inner_pool INT                            Number of processes used to optimize a single data point [default: 1] 
        --outer_pool INT                            Number of data points processed in parallel [default: 1]
        -n, --n_points INT                          Number of starting points [default: 1]
        -i, --iterations INT                        Number of iterations to perform when optimizing [default: 100]
        -x, --xtol_rel FLOAT                        Set relative tolerance on norm of vector of optimisation parameters [default: -1.0]
        -f, --ftol_rel FLOAT                        Set relative tolerance on lnCL [default: -1.0]
        --track                                     Track likelihood search                        
"""
from timeit import default_timer as timer
from docopt import docopt
import sys
import lib.gimble
import lib.math
import zarr
import numpy as np

class OptimizeParameterObj(lib.gimble.ParameterObj):
    '''Sanitises command line arguments and stores parameters.'''

    def __init__(self, params, args):
        super().__init__(params)
        self.zstore = self._get_path(args['--zarr_file'])
        self.data_type = self._get_datatype(args)
        self.config_file = self._get_path(args['--config_file'])
        self.model_file = self._get_path(args['--model_file'])
        #self.threads, self.gridThreads = [self._get_int(t) for t in args["--threads"].split(',')]
        self.threads = self._get_int(args['--inner_pool']) #number of workers for a single set of equations to be solved
        self.gridThreads = self._get_int(args['--outer_pool']) #number of workers for independent processes
        self.numPoints = self._get_int(args['--n_points'])
        self.max_eval = self._get_int(args['--iterations'])
        self.xtol_rel = self._get_float(args['--xtol_rel'])
        self.ftol_rel = self._get_float(args['--ftol_rel'])
        self.trackHistory = args['--track']
        self.config = None
        self.toBeSynced = None
        self.reference = None
        self._parse_config(self.config_file)

    def _get_datatype(self, args):
        choices = [args['--blocks'], args['--windows'], args['--simID']]
        if all(choices) or not any(choices):
            sys.exit("[X] Please specify either '--blocks', '--windows' or '--simID'.")
        if args['--blocks']:
            return 'blocks'
        if args['--windows']:
            return 'windows'
        if args['--simID']:
            self.label = args["--simID"]
            return 'simulate'

def main(params):
    try:
        start_time = timer()
        args = docopt(__doc__)
        parameterObj = OptimizeParameterObj(params, args)
        gimbleStore = lib.gimble.Store(path=parameterObj.zstore, create=False)
        gimbleStore.optimize(parameterObj)
        
        print("[*] Total runtime: %.3fs" % (timer() - start_time))
    except KeyboardInterrupt:
        print("\n[X] Interrupted by user after %s seconds!\n" % (timer() - start_time))
        exit(-1)