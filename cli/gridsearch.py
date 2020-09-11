#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""usage: gimble gridsearch                  -z FILE (-g STR) -c FILE (-b|-w) [-h|--help]
                                            
                                            
    Options:
        -h --help                                   show this
        -z, --zarr_file FILE                        ZARR datastore
        -b, --blocks                                Using blocks
        -w, --windows                               Using windows
        -c, --config_file FILE
        -g, --grid_name STR                         Unique grid hash
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
        self.grid_name = args['--grid_name']
        self.data_type = self._get_datatype([args['--blocks'], args['--windows']])
        self.config_file = self._get_path(args['--config_file'])
        self.config = self._parse_config(self.config_file)
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
        unique_hash = parameterObj._get_unique_hash()
        gimbleStore = lib.gimble.Store(path=parameterObj.zstore, create=False)
        
        grid, meta = gimbleStore._get_grid(unique_hash)
        print(unique_hash) 
        if grid is None:
            sys.exit("[X] Please provide one of the grid(s): %s" % ",".join(gimbleStore.data['grids/']))
        #@Dom this needs to be verified whether the data part is working.
        data = gimbleStore.get_bsfs(
            data_type=parameterObj.data_type, 
            population_by_letter=parameterObj.config['populations'], 
            sample_sets='X', 
            kmax_by_mutype=parameterObj.config['k_max'])
        composite_likelihoods = [lib.math.calculate_composite_likelihood(ETPs, data) for ETPs in grid]
        for idx,L in enumerate(composite_likelihoods):
            print('[+] parameter combination: %s: L=-%s' % (meta[str(idx)], L))
        print("[*] Total runtime: %.3fs" % (timer() - start_time))
    except KeyboardInterrupt:
        print("\n[X] Interrupted by user after %s seconds!\n" % (timer() - start_time))
        exit(-1)