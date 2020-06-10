#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""usage: gimble windows -z <DIR> [-w <INT> -s <INT> -u <INT> -i <INT> -f -D -h]

    Options:
        -z, --zarr <DIR>                            gimble ZARR directory
        -w, --blocks <INT>                          Number of blocks in windows [default: 500]
        -s, --steps <INT>                           Number of steps (blocks) by which windows are shifted [default: 50]
        -u, --max_multiallelic <INT>                Max multiallelics per block [default: 1]
        -i, --max_missing <INT>                     Max missing per block [default: 1]
        
        -f, --force                                 Force overwrite of existing data
        -D, --debug                                 Print debug information
        -h, --help                                   show this

"""

from docopt import docopt
from timeit import default_timer as timer
from lib.gimble import RunObj
import lib.gimble

class ParameterObj(RunObj):
    '''Sanitises command line arguments and stores parameters'''
    def __init__(self, params, args):
        super().__init__(params)
        self.stage = 'windows'
        self.zstore = args['--zarr']
        self.max_multiallelic = int(args['--max_multiallelic'])
        self.max_missing = int(args['--max_missing'])
        self.window_size = int(args['--blocks'])
        self.window_step = int(args['--steps'])
        self.overwrite = True if args['--force'] else False

def main(params):
    try:
        start_time = timer()
        args = docopt(__doc__)
        #log = lib.log.get_logger(run_params)
        parameterObj = ParameterObj(params, args)
        store = lib.gimble.load_store(parameterObj)
        store.make_windows(parameterObj)
        store.dump_windows(parameterObj)
        store.add_stage(parameterObj)
        print("[*] Total runtime: %.3fs" % (timer() - start_time))
    except KeyboardInterrupt:
        print("\n[X] Interrupted by user after %s seconds!\n" % (timer() - start_time))
        exit(-1)