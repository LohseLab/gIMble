#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""usage: gIMble info                  -z FILE [--tree] [-h|--help]
                                            
    Options:
        -z, --zarr_file FILE                        ZARR datastore
        -t, --tree                                  Display datastructure tree

        -h --help                                   show this
"""
from timeit import default_timer as timer
from docopt import docopt
import lib.gimble

class InfoParameterObj(lib.gimble.ParameterObj):
    '''Sanitises command line arguments and stores parameters.'''

    def __init__(self, params, args):
        super().__init__(params)
        self.zstore = self._get_path(args['--zarr_file'])
        self.tree = args['--tree']

def main(params):
    try:
        start_time = timer()
        args = docopt(__doc__)
        parameterObj = InfoParameterObj(params, args)
        gimbleStore = lib.gimble.Store(path=parameterObj.zstore)
        info = gimbleStore.info(parameterObj)
        print(info)
        print("[*] Total runtime: %.3fs" % (timer() - start_time))
    except KeyboardInterrupt:
        print("\n[X] Interrupted by user after %s seconds!\n" % (timer() - start_time))
        exit(-1)