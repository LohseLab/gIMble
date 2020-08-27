#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""usage: gimble optimise                  -z FILE -c FILE (-b|-w) [-h|--help]
                                            
                                            
    Options:
        -h --help                                   show this
        -z, --zarr_file FILE                        ZARR datastore
        -b, --blocks
        -w, --windows
        -c, --config_file FILE
"""
import pathlib
import collections
from timeit import default_timer as timer
from docopt import docopt
import sys
import lib.gimble
import lib.math

class GridsearchParameterObj(lib.gimble.ParameterObj):
    '''Sanitises command line arguments and stores parameters.'''

    def __init__(self, params, args):
        super().__init__(params)
        self.zstore = self._get_path(args['--zarr_file'])
        self.data_type = self._get_datatype([args['--blocks'], args['--windows']])
        self.config_file = self._get_path(args['--config_file'])
        self.config = self._parse_config(self.config_file)

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
        gimbleStore = lib.gimble.Store(path=parameterObj.zstore, create=False)
        #load ETP numpy array
        if not gimbleStore._is_zarr_group(parameterObj.grid_name, 'grids'):
            sys.exit("[X] Specified grid not found in zarr store.")
        grid = zarr.load(gimbleStore.data[f'grids/{parameterObj.grid_name}'])
        data = gimbleStore.get_bsfs_matrix(
            data=parameterObj.data_type, 
            population_by_letter=parameterObj.config['population_ids'], 
            cartesian_only=True, 
            kmax_by_mutype=parameterObj.config['k_max'])

        composite_likelihoods = [lib.math.calculate_composite_likelihood(ETPs, data) for ETPs in grid]
        for idx, L in enumerate(composite_likelihoods):
            print('[+] parameter combination: %s: L=-%s' % (idx, L))
        print("[*] Total runtime: %.3fs" % (timer() - start_time))
    except KeyboardInterrupt:
        print("\n[X] Interrupted by user after %s seconds!\n" % (timer() - start_time))
        exit(-1)