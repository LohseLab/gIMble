#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""usage: gimble optimise                  [-z FILE] -c FILE [-m FILE] [-b|-w|-e] [-t STR] [-n INT] 
                                            [-x FLOAT -i INT] [-f FLOAT] [-p]
                                            [-h|--help]
                                            
                                            
    Options:
        -h --help                                   show this
        -z, --zarr_file FILE                        ZARR datastore
        -c, --config_file FILE                      INI config file
        -m, --model_file FILE                       gimble model TSV
        -b, --blocks                                Optimise based on blocks 
        -w, --windows                               Optimise based on windows (might take very long)
        -e, --etp                                   Optimise based on model (debugging only)
        -t, --threads STR                           Threads [default: 1,1]
        -n, --n_points INT                          Number of starting points [default: 1]
        -i, --iterations INT                        Number of iterations to perform when optimizing [default: 100]
        -x, --xtol_rel FLOAT                        Set relative tolerance on optimisation parameters [default: 0.00000001]
        -f, --ftol_rel FLOAT                        Set relative tolerance on CL [default: 0.000000001]
        -p, --trackPath                             Track likelihood search                        
"""
from timeit import default_timer as timer
from docopt import docopt
import sys
import lib.gimble
import lib.math
import zarr
import numpy as np

class OptimiseParameterObj(lib.gimble.ParameterObj):
    '''Sanitises command line arguments and stores parameters.'''

    def __init__(self, params, args):
        super().__init__(params)
        self.zstore = self._get_path(args['--zarr_file'])
        self.data_type = self._get_datatype(args)
        self.config_file = self._get_path(args['--config_file'])
        self.model_file = self._get_path(args['--model_file'])
        self.threads, self.gridThreads = [self._get_int(t) for t in args["--threads"].split(',')]
        #self.threads, self.gridThreads = self._get_threads(args["--threads"])
        self.config = self._parse_config(self.config_file)
        self.numPoints = self._get_int(args['--n_points'])
        self.iterations = self._get_int(args['--iterations'])
        self.xtol_rel = self._get_float(args['--xtol_rel'])
        self.ftol_rel = self._get_float(args['--ftol_rel'])
        self.trackHistory = args['--trackPath']
        self._process_config()

    def _get_datatype(self, args):
        choices = [args['--blocks'], args['--windows'], args['--etp']]
        if all(choices) or not any(choices):
            sys.exit("[X] Please specify either '--blocks' or '--windows' or '--etp.")
        if args['--blocks']:
            return 'blocks'
        if args['--windows']:
            return 'windows'
        if args['--etp']:
            return 'etp'
        
def main(params):
    try:
        start_time = timer()
        args = docopt(__doc__)
        parameterObj = OptimiseParameterObj(params, args)
        print("[+] Generated all parameter combinations.")
        if parameterObj.data_type == 'etp':
            #depends on Dom will have resolved the zarr structuring
            #sys.exit("[X] this part of the code needs to be replaced with the path to the actual data.")
            z=zarr.open('output/test.z')
            data =np.array(z['grids/grid_7'][0])    
        else:
            gimbleStore = lib.gimble.Store(path=parameterObj.zstore, create=False)
            if not gimbleStore.has_stage(parameterObj.data_type):
                sys.exit("[X] GStore has no %r." % parameterObj.data_type)
            data = gimbleStore.get_bsfs(data_type=parameterObj.data_type, population_by_letter=parameterObj.config['populations'], sample_sets="X", kmax_by_mutype=parameterObj.config['k_max'])
            print(data)
        # load math.EquationSystemObj
        equationSystem = lib.math.EquationSystemObj(parameterObj)
        # initiate model equations
        equationSystem.initiate_model(parameterObj)
        optimizeResult=equationSystem.optimize_parameters(
            data, 
            maxeval=parameterObj.iterations, 
            xtol_rel=parameterObj.xtol_rel, 
            ftol_rel=parameterObj.ftol_rel, 
            numPoints=parameterObj.numPoints, 
            threads=parameterObj.threads, 
            gridThreads=parameterObj.gridThreads, 
            trackHistory=parameterObj.trackHistory)
        
        print("[*] Total runtime: %.3fs" % (timer() - start_time))
    except KeyboardInterrupt:
        print("\n[X] Interrupted by user after %s seconds!\n" % (timer() - start_time))
        exit(-1)