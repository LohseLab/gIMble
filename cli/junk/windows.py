#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""usage: gimble windows -z <DIR> [-w <INT> -s <INT> -u <INT> -i <INT> -D -h]

    Options:
        -z, --zarr <DIR>                            gimble ZARR directory
        -w, --blocks <INT>                          Number of blocks in windows [default: 500]
        -s, --steps <INT>                           Number of steps (blocks) by which windows are shifted [default: 50]
        
        -u, --max_multiallelic <INT>                Max multiallelics per block [default: 1]
        -i, --max_missing <INT>                     Max missing per block [default: 1]
        
        -D, --debug                                 Print debug information
        -h, --help                                   show this

"""

from docopt import docopt
from timeit import default_timer as timer
#from sys import stderr, exit
import lib.gimble
'''

'''

class ParameterObj(object):
    def __init__(self, args):
        print(args)
        self.zstore = args['--zarr']
        self.max_multiallelic = int(args['--max_multiallelic'])
        self.max_missing = int(args['--max_missing'])
        self.window_size = int(args['--blocks'])
        self.window_step = int(args['--steps'])

def main(run_params):
    try:
        start_time = timer()
        args = docopt(__doc__)
        #log = lib.log.get_logger(run_params)
        parameterObj = ParameterObj(args)
        store = lib.gimble.load_store(parameterObj)
        #print(store.tree())
        store.make_windows(parameterObj)
        store.dump_windows(parameterObj)
        #print(store.tree())
        #print(store.attrs())
        
        #log.info("[*] Total runtime: %.3fs" % (timer() - start_time))
    except KeyboardInterrupt:
        #log.info("\n[X] Interrupted by user after %s seconds!\n" % (timer() - start_time))
        exit(-1)

###############################################################################

if __name__ == '__main__':
    main()