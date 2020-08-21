#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""usage: gimble query                  -z FILE [-b -w --extended_bed --popgen_metrics] 
                                            [-h|--help]
                                            
    Options:
        -h --help                                   show this
        -z, --zarr_file FILE                        ZARR datastore
        -b, --blocks                                Writes BED file of blocks
        -w, --windows                               Writes BED file of windows [TBI]
        --extended_bed                              Adds variation/multiallelic/missing to BED
        --popgen_metrics                            [TBI]

"""

'''
-u, --max_multiallelic <INT>                Max multiallelics per block [default: 2]
        -i, --max_missing <INT>                     Max missing per block [default: 2]
-W, --windows_bsfs                          Writes bSFS of windows
-B, --blocks_bsfs                           Writes bSFS of blocks
--popgen_metrics -B
causes
    - global bSFS tally as tsv
    - samble-set bSFS tally as tsv
    - popgen-metrics file
        - overall pi_A, pi_B, dxy, fst
            intra: [(0.5*m_1 + 0.5*m_2 + m_3+m_4) / blocks for sample_set in intra_sample_sets] / len(intra_sample_sets)
        - sample-set pi_A, pi_B, dxy, fst


'''
from timeit import default_timer as timer
from docopt import docopt
import lib.gimble

class QueryParameterObj(lib.gimble.ParameterObj):
    '''Sanitises command line arguments and stores parameters.'''

    def __init__(self, params, args):
        super().__init__(params)
        self.zstore = self._get_path(args['--zarr_file'])
        self.blocks = args['--blocks']
        #self.blocks_bsfs = args['--blocks_bsfs']
        self.windows = args['--windows']
        #self.windows_bsfs = args['--windows_bsfs']
        self.extended_bed = args['--extended_bed']
        self.popgen_metrics = args['--popgen_metrics']

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