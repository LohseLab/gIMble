#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""usage: gimble windows -z <DIR> [-D -h]

    Options:
        -z, --zarr <DIR>                            gimble ZARR directory
        -b, --blocks <INT>                          Number of blocks in windows [default: 500]
        -s, --steps <INT>                           Number of steps (blocks) by which windows are shifted [default: 50]
        
        -D, --debug                                 Print debug information
        -h --help                                   show this

"""

from docopt import docopt
from timeit import default_timer as timer
#from sys import stderr, exit
import lib.setup
import lib.classes
import lib.log
import lib.analysis

'''

'''

def main(run_params):
    try:
        start_time = timer()
        args = docopt(__doc__)
        log = lib.log.get_logger(run_params)
        parameterObj = lib.analysis.ParameterObj(args)
        store = lib.classes.Store(parameterObj)
        store.make_windows(parameterObj)
        print(store.tree())
        print(store.attrs())
        store.dump_windows()
        log.info("[*] Total runtime: %.3fs" % (timer() - start_time))
    except KeyboardInterrupt:
        log.info("\n[X] Interrupted by user after %s seconds!\n" % (timer() - start_time))
        exit(-1)

###############################################################################

if __name__ == '__main__':
    main()