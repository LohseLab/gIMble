#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""usage: gIMble info                  -z FILE [-h|--help]
                                            
    Options:
        -h --help                                   show this
        -z, --zarr_file FILE                        ZARR datastore

"""
from timeit import default_timer as timer
from docopt import docopt
import lib.gimble

class ParameterObj(lib.gimble.RunObj):
    '''Sanitises command line arguments and stores parameters.'''

    def __init__(self, params, args):
        super().__init__(params)
        self.zstore = self._get_path(args['--zarr_file'])

def main(params):
    try:
        start_time = timer()
        args = docopt(__doc__)
        parameterObj = ParameterObj(params, args)
        store = lib.gimble.load_store(parameterObj)
        #print(store.tree())
        store.info(verbose=True)
        print("[*] Total runtime: %.3fs" % (timer() - start_time))
    except KeyboardInterrupt:
        print("\n[X] Interrupted by user after %s seconds!\n" % (timer() - start_time))
        exit(-1)