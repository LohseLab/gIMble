#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""usage: gimble optimise                  [-z FILE] -c FILE [-m FILE] [-b|-w|-e] [-t INT] [-h|--help]
                                            
                                            
    Options:
        -h --help                                   show this
        -z, --zarr_file FILE                        ZARR datastore
        -c, --config_file FILE                      INI config file
        -m, --model_file FILE                       gimble model TSV
        -b, --blocks                                Optimise based on blocks 
        -w, --windows                               Optimise based on windows (might take very long)
        -e, --etp                                   Optimise based on model (debugging only)
        -t, --threads INT                           Threads [default: 1]
"""
from timeit import default_timer as timer
from docopt import docopt
import sys
import lib.gimble
import lib.math
import zarr

class OptimiseParameterObj(lib.gimble.ParameterObj):
    '''Sanitises command line arguments and stores parameters.'''

    def __init__(self, params, args):
        super().__init__(params)
        self.zstore = self._get_path(args['--zarr_file'])
        self.data_type = self._get_datatype(args)
        self.config_file = self._get_path(args['--config_file'])
        self.model_file = self._get_path(args['--model_file'])
        self.threads = self._get_int(args["--threads"])
        self.config = self._parse_config(self.config_file)
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
            # load math.EquationSystemObj
            equationSystem = lib.math.EquationSystemObj(parameterObj)
            # initiate model equations
            equationSystem.initiate_model(parameterObj)
            equationSystem.calculate_all_ETPs()
            equationSystem.optimise_parameters(equationSystem.ETPs, maxeval=50, localOptimum=False)
        else:
            gimbleStore = lib.gimble.Store(path=parameterObj.zstore, create=False)
            if not gimbleStore.has_stage(parameterObj.data_type):
                sys.exit("[X] GStore has no %r." % parameterObj.data_type)
            # load math.EquationSystemObj
            equationSystem = lib.math.EquationSystemObj(parameterObj)
            # initiate model equations
            equationSystem.initiate_model(parameterObj)
            data = gimbleStore.get_bsfs(data_type=parameterObj.data_type, sample_sets="X", kmax_by_mutype=parameterObj.config['k_max'])
            equationSystem.optimise_parameters(data, maxeval=50, localOptimum=False)
        #load ETP numpy array
        #if not gimbleStore._is_zarr_group(parameterObj.grid_name, 'grids'):
        #    sys.exit("[X] Specified grid not found in zarr store.")
        #grid = zarr.load(gimbleStore.data[f'grids/{parameterObj.grid_name}'])
        #data = gimbleStore.get_bsfs_matrix(
        #    data=parameterObj.data_type, 
        #    population_by_letter=parameterObj.config['population_ids'], 
        #    cartesian_only=True, 
        #    kmax_by_mutype=parameterObj.config['k_max'])
        #load data
        #z=zarr.open('/Users/s1854903/Documents/ongoing_projects/gIMble/output/new_configs/test.z')
        #data =z['grids/grid_7'][0]
        #run optimisation
        #equationSystem.optimise_parameters(data, maxeval=50, localOptimum=False)
        
        print("[*] Total runtime: %.3fs" % (timer() - start_time))
    except KeyboardInterrupt:
        print("\n[X] Interrupted by user after %s seconds!\n" % (timer() - start_time))
        exit(-1)