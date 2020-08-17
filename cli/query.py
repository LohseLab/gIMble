#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""usage: gimble query                  -z FILE -b [-u <INT> -i <INT>] [-h|--help]
                                            
    Options:
        -h --help                                   show this
        -z, --zarr_file FILE                        ZARR datastore
        -b, --blocks                                Writes BED file of blocks
        -B, --blocks_bsfs                           Writes bSFS of blocks
        -w, --windows                               Writes BED file of windows
        -W, --windows_bsfs                          Writes bSFS of windows
        -u, --max_multiallelic <INT>                Max multiallelics per block [default: 2]
        -i, --max_missing <INT>                     Max missing per block [default: 2]
        --popgen_metrics
"""

'''
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

class ParameterObj(lib.gimble.RunObj):
    '''Sanitises command line arguments and stores parameters.'''

    def __init__(self, params, args):
        super().__init__(params)
        self.zstore = self._get_path(args['--zarr_file'])
        self.blocks = args['--blocks']
        self.block_max_multiallelic = int(args['--max_multiallelic'])
        self.block_max_missing = int(args['--max_missing'])

def main(params):
    try:
        start_time = timer()
        args = docopt(__doc__)
        parameterObj = ParameterObj(params, args)
        store = lib.gimble.load_store(parameterObj)
        #print(store.tree())
        store.write_block_bed(parameterObj)
        print("[*] Total runtime: %.3fs" % (timer() - start_time))
    except KeyboardInterrupt:
        print("\n[X] Interrupted by user after %s seconds!\n" % (timer() - start_time))
        exit(-1)