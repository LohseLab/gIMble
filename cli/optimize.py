#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""usage: gimble optimize                  [-z FILE] -c FILE [-m FILE] [-b|-w] [-t STR] [-h|--help] [-n INT]
                                            
                                            
    Options:
        -h --help                                   show this
        -z, --zarr_file FILE                        ZARR datastore
        -b, --blocks
        -w, --windows
        -c, --config_file FILE
        -m, --model_file FILE
        -t, --threads STR                           Threads [default: 1,1]
        -n, --n_points INT                          Number of starting points [default: 1]
"""
import pathlib
import collections
from timeit import default_timer as timer
from docopt import docopt
import sys
import lib.gimble
import lib.math
import zarr

class GridsearchParameterObj(lib.gimble.ParameterObj):
    '''Sanitises command line arguments and stores parameters.'''

    def __init__(self, params, args):
        super().__init__(params)
        self.zstore = self._get_path(args['--zarr_file'])
        self.data_type = self._get_datatype([args['--blocks'], args['--windows']])
        self.config_file = self._get_path(args['--config_file'])
        self.model_file = self._get_path(args['--model_file'])
        self.threads, self.gridThreads = [self._get_int(t) for t in args["--threads"].split(',')]
        self.config = self._parse_config(self.config_file)
        self.numPoints = self._get_int(args['--n_points'])
        self._process_config()

    def _get_datatype(self, args):
        if not any(args):
            return None
        elif args[0]:
            return 'blocks'
        elif args[1]:
            return 'windows'
        else:
            sys.exit("[X1] This should not have happened.")
        
def main(params):
    try:
        start_time = timer()
        args = docopt(__doc__)
        parameterObj = GridsearchParameterObj(params, args)
        print("[+] Generated all parameter combinations.")
        equationSystem = lib.math.EquationSystemObj(parameterObj)
        equationSystem.initiate_model(parameterObj)
        #load gimblestore and stored ETP numpy array
        #gimbleStore = lib.gimble.Store(path=parameterObj.zstore, create=False)
        #load ETP numpy array
        #if not gimbleStore._is_zarr_group(parameterObj.grid_name, 'grids'):
        #    sys.exit("[X] Specified grid not found in zarr store.")
        #grid = zarr.load(gimbleStore.data[f'grids/{parameterObj.grid_name}'])
        #data = gimbleStore.get_bsfs_matrix(
        #    data=parameterObj.data_type, 
        #    population_by_letter=parameterObj.config['population_ids'], 
        #    cartesian_only=True, 
        #    kmax_by_mutype=parameterObj.config['k_max'])
        #load data: in test.z in output folder
        z=zarr.open('../output/test.z')
        data =z['grids/grid_7'][0]
        #run optimisation
        equationSystem.optimize_parameters(data, maxeval=10, xtol_rel=0.01, numPoints=parameterObj.numPoints, threads=parameterObj.threads, gridThreads=parameterObj.gridThreads)
        
        print("[*] Total runtime: %.3fs" % (timer() - start_time))
    except KeyboardInterrupt:
        print("\n[X] Interrupted by user after %s seconds!\n" % (timer() - start_time))
        exit(-1)