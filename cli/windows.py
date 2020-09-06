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
        self.window_size = int(args['--blocks'])
        self.window_step = int(args['--steps'])
        self.overwrite = args['--force']

def main(params):
    try:
        start_time = timer()
        print("[+] Running 'gimble windows'")
        args = docopt(__doc__)
        #log = lib.log.get_logger(run_params)
        parameterObj = WindowsParameterObj(params, args)
        gimbleStore = lib.gimble.Store(path=parameterObj.zstore, create=False)
        #gimbleStore.windows(parameterObj)
        gimbleStore.info()
        #gimbleStore.tree()
        import numpy as np
        start_time = timer()
        bsfs = gimbleStore.get_block_bsfs(sample_sets='X')
        print("bsfs.shape", bsfs.shape)
        print("np.sum(bsfs)", np.sum(bsfs))
        # for sequence, bsfs in gimbleStore.yield_window_bsfs_by_seq(
        #     #):
        #     kmax_by_mutype={'m_1': 2, 'm_2': 2, 'm_3': 2, 'm_4': 2}):
        #     print("sequence", sequence, bsfs.shape)
        #     for m, c in lib.gimble.bsfs_to_counter(bsfs).items():
        #         print(m, c)
        print("[*] Total runtime: %.3fs" % (timer() - start_time))
    except KeyboardInterrupt:
        print("\n[X] Interrupted by user after %s seconds!\n" % (timer() - start_time))
        exit(-1)