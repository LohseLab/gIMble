#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""usage: gimble windows -z <DIR> [-w <INT> -s <INT> -u <INT> -i <INT> -f -D -h]

    Options:
        -z, --zarr <DIR>                            gimble ZARR directory
        -w, --blocks <INT>                          Number of blocks in windows [default: 500]
        -s, --steps <INT>                           Number of steps (blocks) by which windows are shifted [default: 50]
        
        -f, --force                                 Force overwrite of existing data
        -D, --debug                                 Print debug information
        -h, --help                                   show this

"""

from docopt import docopt
from timeit import default_timer as timer
import lib.gimble

class WindowsParameterObj(lib.gimble.ParameterObj):
    '''Sanitises command line arguments and stores parameters'''
    def __init__(self, params, args):
        super().__init__(params)
        self.zstore = self._get_path(args['--zarr'])
        self.window_size = self._get_int((args['--blocks']))
        self.window_step = self._get_int((args['--steps']))
        self.overwrite = True if args['--force'] else False

def main(params):
    try:
        start_time = timer()
        args = docopt(__doc__)
        parameterObj = WindowsParameterObj(params, args)
        gimbleStore = lib.gimble.Store(path=parameterObj.zstore, create=False)
        gimbleStore.windows(
            window_size=parameterObj.window_size, 
            window_step=parameterObj.window_step, 
            overwrite=parameterObj.overwrite)
        gimbleStore.log_action(module=parameterObj._MODULE, command=parameterObj._get_cmd())
        print("[*] Total runtime was %s" % (lib.gimble.format_time(timer() - start_time)))
    except KeyboardInterrupt:
        print("\n[X] Interrupted by user after %s seconds!\n" % (timer() - start_time))
        exit(-1)