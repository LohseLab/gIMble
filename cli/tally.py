#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""usage: gimbl tally -z <FILE> -s <STR> -d <STR> [-m <STR>] [-f] [-h|--help]
                                            
                                            
    Options:
        -h --help                                   show this
        
        -z, --zarr_file <FILE>                      Path to existing GimbleStore
        -s, --data_source <STR>                     'blocks' OR 'windows' 
        -d, --data_label <STR>                      Data label under which the tally gets saved
        -m, --maxk STR                              Max value for mutypes (values above get binned) [default: '2,2,2,2']
        -f, --overwrite                             Overwrite results in GimbleStore

"""
from timeit import default_timer as timer
from docopt import docopt
import lib.gimble
import lib.math

class TallyParameterObj(lib.gimble.ParameterObj):
    '''Sanitises command line arguments and stores parameters.'''

    def __init__(self, params, args):
        super().__init__(params)
        self.zstore = self._get_path(args['--zarr_file'])
        self.data_source = args['--data_source']
        self.data_label = args['--data_label']
        self.max_k = self._get_max_k(args['--maxk'])
        self.overwrite = args['--overwrite']
        
    def _get_max_k(self, kmax_string):
        #mutypes = ['m_1', 'm_2', 'm_3', 'm_4']
        if kmax_string == None:
            return None
        try:
            return np.array(ast.literal_eval(kmax_string))
        except ValueError:
            sys.exit("[X] Invalid k-max string (must be python list)") 
def main(params):
    try:
        start_time = timer()
        args = docopt(__doc__)
        parameterObj = GridsearchParameterObj(params, args)
        gimbleStore = lib.gimble.Store(path=parameterObj.zstore, create=False)
        gimbleStore.gridsearch(
            data_label=parameterObj.data_label,
            grid_label=parameterObj.grid_label,
            overwrite=parameterObj.overwrite,
            )
        print("[*] Total runtime was %s" % (lib.gimble.format_time(timer() - start_time)))
    except KeyboardInterrupt:
        print("\n[X] Interrupted by user after %s !\n" % (lib.gimble.format_time(timer() - start_time)))
        exit(-1)