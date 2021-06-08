#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""usage: gimble gridsearch -z <FILE> -g <STR> -o <STR> (-b | -w | -s <STR>) [-f] [-h|--help]
                                            
                                            
    Options:
        -h --help                                   show this
        
        -z, --zarr_file <FILE>                      Path to existing GimbleStore
        -b, --blocks                                Using blocks
        -w, --windows                               Using windows
        -g, --grid_label <STR>                          Label of makegrid run in GimbleStore
        -s, --simulation_label <STR>                    Label of simulate run in GimbleStore
        -o, --output_label <STR>                        Label under which results get stored in GimbleStore
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
        self.data_type = self._get_datatype(args)
        self.config_file = self._get_path(args['--config_file'])
        self.overwrite = args['--overwrite']
        self.config = None
        self._parse_config(self.config_file)

    def _get_datatype(self, args):
        if args['--blocks']:
            return 'blocks'
        if args['--windows']:
            return 'windows'
        if args['--simulation_label']:
            self.label = args['--simulation_label']
            return 'simulate'
        return None
        
def main(params):
    try:
        start_time = timer()
        args = docopt(__doc__)
        parameterObj = GridsearchParameterObj(params, args)
        gimbleStore = lib.gimble.Store(path=parameterObj.zstore, create=False)
        gimbleStore.gridsearch(parameterObj)
        print("[*] Total runtime was %s" % (lib.gimble.format_time(timer() - start_time)))
    except KeyboardInterrupt:
        print("\n[X] Interrupted by user after %s !\n" % (lib.gimble.format_time(timer() - start_time)))
        exit(-1)