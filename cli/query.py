#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""usage: gimble query                  -z FILE [-b -w -s] [--bed|--simplebed] [-c FILE --lncls] [--bsfs --kmax STR] 
                                            [-h|--help]
                                            
    Options:
        -h --help                                   show this
        -z, --zarr_file FILE                        ZARR datastore
        -b, --blocks                                Query data for blocks
        -w, --windows                               Query data for windows
        -s, --windows_sum                           Query data for sum of windows
        -c, --config_file FILE                      INI Config file
        --bed                                       Writes BED with variation/multiallelic/missing (for -b, -w)
        --simplebed                                 Writes only BED coordinates (for -w)
        --bsfs                                      Writes 2D bSFS for data (for -b, -s)
        --lncls                                     Write grid and lnCls based on an INI file (requires)
        --kmax STR                                  k-max for mutypes, e.g. [2, 2, 2, 2]

"""

'''
query based on grid idx/params


output CSV
    for each windows best
'''

from timeit import default_timer as timer
from docopt import docopt
import lib.gimble
import sys
import itertools 
import ast

class QueryParameterObj(lib.gimble.ParameterObj):
    '''Sanitises command line arguments and stores parameters.'''

    def __init__(self, params, args):
        super().__init__(params)
        self.data_type = self._get_data_type(args)
        self.data_format = self._get_data_format(args)
        self._check_input()
        self.zstore = self._get_path(args['--zarr_file'])
        self.kmax_by_mutype = self._get_kmax(args['--kmax'])
        self.config_file = self._get_path(args['--config_file'])
        self.config = None
        self._parse_config(self.config_file)

    def _check_input(self):
        if self.data_type == 'blocks':
            if self.data_format == 'bsfs' or self.data_format == 'bed':
                return True
        if self.data_type == 'windows':
            if self.data_format == 'bed' or self.data_format == 'lncls' or self.data_format == 'simplebed':
                return True
        if self.data_type == 'windows_sum':
            if self.data_format == 'bsfs' or self.data_format == 'lncls':
                return True
        sys.exit("[X] Data type %r and data format %r are not supported" % (self.data_type, self.data_format))

    def _get_data_type(self, args):
        if args['--blocks']:
            return 'blocks'
        elif args['--windows']:
            return 'windows'
        elif args['--windows_sum']:
            return 'windows_sum'
        else:
            raise ValueError('unknown data_type')

    def _get_data_format(self, args):
        if args['--bed']:
            return 'bed'
        elif args['--simplebed']:
            return 'simplebed'
        elif args['--bsfs']:
            return 'bsfs'
        elif args['--lncls']:
            return 'lncls'
        else:
            raise ValueError('unknown data_format')

    def _get_kmax(self, kmax_string):
        mutypes = ['m_1', 'm_2', 'm_3', 'm_4']
        if kmax_string == None:
            return []
        try:
            return {mutype: kmax for mutype, kmax in zip(mutypes, ast.literal_eval(kmax_string))}
        except ValueError:
            sys.exit("[X] Invalid k-max string (must be python list)") 

def main(params):
    try:
        start_time = timer()
        args = docopt(__doc__)
        parameterObj = QueryParameterObj(params, args)
        gimbleStore = lib.gimble.Store(path=parameterObj.zstore)
        gimbleStore.query(parameterObj)
        print("[*] Total runtime: %.3fs" % (timer() - start_time))
    except KeyboardInterrupt:
        print("\n[X] Interrupted by user after %s seconds!\n" % (timer() - start_time))
        exit(-1)