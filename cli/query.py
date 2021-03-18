"""usage: gimbl query                    -z DIR [-b|-w|-s] [--bed|--tally|--lncls] [-m STR] [-k STR] 
                                            [-h|--help]
                                            
        -z, --zarr_file DIR              ZARR datastore
        -b, --blocks                     Query data for blocks
        -w, --windows                    Query data for windows
        -s, --windows_sum                Query data for sum of windows
        --bed                            Writes BED with variation/multiallelic/missing (windows)
        --tally                          Writes 2D bSFS for data (blocks/windows/windows_sum)
        --lncls                          Write grid and lnCls (for -k <STR>)
        -m, --maxk STR                   Max value for mutypes (values above get binned), e.g. [2, 2, 2, 2]
        -k, --key STR                    Query data based on key

        -h --help                        show this
"""

from timeit import default_timer as timer
from docopt import docopt
import lib.gimble
import numpy as np
import sys
import ast

class QueryParameterObj(lib.gimble.ParameterObj):
    '''Sanitises command line arguments and stores parameters.'''

    def __init__(self, params, args):
        super().__init__(params)
        self.zstore = self._get_path(args['--zarr_file'])
        self.data_type = self._get_data_type(args)
        self.data_format = self._get_data_format(args)
        self.max_k = self._get_max_k(args['--maxk'])
        self.key = args['--key']
        self._check_input(args)

    def _check_input(self, args):
        if self.data_type == 'blocks':
            if self.data_format == 'tally':
                return True
        if self.data_type == 'windows':
            if self.data_format == 'bed' or self.data_format == 'lncls' or self.data_format == 'tally':
                return True
        if self.data_type == 'windows_sum':
            if self.data_format == 'tally' or self.data_format == 'lncls':
                return True
        if self.data_type == None:
            sys.exit(__doc__)
        sys.exit("[X] Data type %r and data format %r are not supported" % (self.data_type, self.data_format))

    def _get_data_type(self, args):
        if args['--blocks']:
            return 'blocks'
        elif args['--windows']:
            return 'windows'
        elif args['--windows_sum']:
            return 'windows_sum'
        else:
            return None

    def _get_data_format(self, args):
        if args['--bed']:
            return 'bed'
        elif args['--tally']:
            return 'tally'
        elif args['--lncls']:
            return 'lncls'
        else:
            return None

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
        parameterObj = QueryParameterObj(params, args)
        gimbleStore = lib.gimble.Store(path=parameterObj.zstore)
        gimbleStore.query(
            parameterObj._VERSION,
            parameterObj.data_type,
            parameterObj.data_format,
            parameterObj.max_k
            )
        print("[*] Total runtime was %s" % (lib.gimble.format_time(timer() - start_time)))
    except KeyboardInterrupt:
        print("\n[X] Interrupted by user after %s !\n" % (lib.gimble.format_time(timer() - start_time)))
        exit(-1)