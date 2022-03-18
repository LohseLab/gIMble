#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""usage: gimble gridsearch -z <FILE> -g <STR> [-d <STR>|-s <STR>] [-f] [-h|--help]
                                            
                                            
    Options:
        -h --help                                   show this
        
        -z, --zarr_file <FILE>                      Path to existing GimbleStore
        -g, --grid_label <STR>                      Label of makegrid run in GimbleStore
        -d, --tally_label <STR>                     Tally label
        -s, --sim_label <STR>                       Simulation label
        -f, --overwrite                             Overwrite results in GimbleStore

"""
from timeit import default_timer as timer
from docopt import docopt
import lib.gimble
import lib.math

class GridsearchParameterObj(lib.gimble.ParameterObj):
    '''Sanitises command line arguments and stores parameters.'''

    def __init__(self, params, args):
        super().__init__(params)
        self.zstore = self._get_path(args['--zarr_file'])
        self.tally_label = args['--tally_label']
        self.sim_label = args['--sim_label']
        self.grid_label = args['--grid_label']
        self.overwrite = args['--overwrite']
        
def main(params):
    try:
        start_time = timer()
        args = docopt(__doc__)
        parameterObj = GridsearchParameterObj(params, args)
        gimbleStore = lib.gimble.Store(path=parameterObj.zstore, create=False)
        gimbleStore.gridsearch(
            tally_label=parameterObj.tally_label,
            sim_label=parameterObj.sim_label,
            grid_label=parameterObj.grid_label,
            overwrite=parameterObj.overwrite,
            )
        print("[*] Total runtime was %s" % (lib.gimble.format_time(timer() - start_time)))
    except KeyboardInterrupt:
        print("\n[X] Interrupted by user after %s !\n" % (lib.gimble.format_time(timer() - start_time)))
        exit(-1)