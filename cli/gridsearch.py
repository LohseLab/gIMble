#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""usage: gimble gridsearch                  -z FILE [-g STR] -c FILE (-b|-w) [-h|--help]
                                            
                                            
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
        grids, grid_meta_dict = gimbleStore._get_grid(unique_hash) 
        if grids is None:
            sys.exit("[X] Please provide one of the grid(s): %s" % ",".join(gimbleStore.data['grids/']))
        data = gimbleStore.get_bsfs(
            data_type=parameterObj.data_type, 
            population_by_letter=parameterObj.config['population_by_letter'], 
            sample_sets='X', 
            kmax_by_mutype=parameterObj.config['k_max'],
            as_dask=True)
        # numpy
        gridsearch_result = gimbleStore.gridsearch_np(data=data, grids=grids)
        
        # dask
        #gridsearch_result = gimbleStore.gridsearch(data=data, grids=grids)

        # test_dask
        #gridsearch_result = gimbleStore.test_dask(data=data, grids=grids)
    
        output_f = gimbleStore._write_gridsearch_bed(parameterObj=parameterObj, data=gridsearch_result, grid_meta_dict=grid_meta_dict)
        print("[+] Wrote %r." % output_f)
        print("[*] Total runtime: %.3fs" % (timer() - start_time))
        #print('gridsearch_result', gridsearch_result.shape)
        #output_f = gimbleStore._write_gridsearch_bed(parameterObj=parameterObj, data=gridsearch_result, grid_meta_dict=grid_meta_dict)
        #print("[+] Wrote %r." % output_f)
        #window_params = lib.math.get_window_params(meta=meta, gridsearch_results=gridsearch_results)
        #from dask.distributed import Client
        #client = Client()  # start local workers as processes
        #client = Client(processes=False)  # start local workers as threads
        #for parameters in itertools.product(grid, ):
        #    lazy_result = dask.delayed(costly_simulation)(parameters)
        #    lazy_results.append(lazy_result)
        #dask.array.apply_over_axes(np.sum, dask.array.from_array(data, chunks=(1,4,4,4,4)), axes = [1,2,3,4]).compute()
        #composite_likelihoods = [lib.math.calculate_composite_likelihood(ETPs, data) for ETPs in grid]
        #for idx,L in enumerate(composite_likelihoods):
        #    print('[+] parameter combination: %s: L=-%s' % (meta[str(idx)], L))
        
    except KeyboardInterrupt:
        print("\n[X] Interrupted by user after %s seconds!\n" % (timer() - start_time))
        exit(-1)