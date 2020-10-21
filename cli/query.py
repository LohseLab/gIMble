#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""usage: gimble query                  -z FILE [-b -w -s] [--bed --bsfs --kmax STR] 
                                            [-h|--help]
                                            
    Options:
        -h --help                                   show this
        -z, --zarr_file FILE                        ZARR datastore
        -b, --blocks                                Query data for blocks
        -w, --windows                               Query data for windows
        -s, --windows_sum                           Query data for sum of windows
        --bed                                       Writes BED variation/multiallelic/missing to BED (for -b, -w)
        --bsfs                                      Writes 2D bSFS for data (for -b, -s)
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
        data_choice = [args['--blocks'], args['--windows'], args['--windows_sum']]
        format_choice = [args['--bed'], args['--bsfs']]
        self.queries = self._get_queries(data_choice, format_choice)
        self.zstore = self._get_path(args['--zarr_file'])
        self.kmax_by_mutype = self._get_kmax(args['--kmax'])

    def _get_kmax(self, kmax_string):
        mutypes = ['m_1', 'm_2', 'm_3', 'm_4']
        if kmax_string == None:
            return []
        try:
            return {mutype: kmax for mutype, kmax in zip(mutypes, ast.literal_eval(kmax_string))}
        except ValueError:
            sys.exit("[X] Invalid k-max string (must be python list)") 

    def _get_queries(self, data_choice, format_choice):
        data_choice = [choice for arg, choice in zip(data_choice, ['blocks', 'windows', 'windows_sum']) if arg]
        format_choice = [choice for arg, choice in zip(format_choice, ['bed', 'bsfs']) if arg] 
        queries = []
        for query in itertools.product(data_choice, format_choice):
            data_type, format_type = query
            if (data_type == 'windows_sum' and format_type == 'bed') or (data_type == 'windows' and format_type == 'bsfs'):
                sys.exit("[X] %r and %r are incompatible." % (data_type, format_type))
            else:
                queries.append(query)
        return queries

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