"""usage: gIMble blocks         -z <DIR> [-l <INT> -m <INT> -u <INT> -i <INT> -d -f -h]
    
    -z, --zarr <DIR>                            gIMble ZARR directory
    
    -l, --block_length <INT>                    Successively genotyped sites per block [default: 64] 
    -m, --block_span <INT>                      Maximum distance between first and last site of a block (default: '-l' * 2)
    -u, --max_multiallelic <INT>                Max multiallelics per block [default: 2]
    -i, --max_missing <INT>                     Max missing per block [default: 2]

    -f, --force                                 Force overwrite of existing data
    -d, --debug                                 Write debugging logs
    -h, --help
"""
import sys
from timeit import default_timer as timer
from docopt import docopt
#import lib.gimblelog
from lib.gimble import RunObj
import lib.gimble


class ParameterObj(RunObj):
    '''Sanitises command line arguments and stores parameters'''

    def __init__(self, params, args):
        super().__init__(params)
        self.zstore = self._get_path(args['--zarr'])
        self.block_length = self._get_int(args['--block_length'])
        self.block_span = self._get_block_span(args['--block_span'])
        self.block_gap_run = self._get_block_gap_run()
        self.block_max_multiallelic = int(args['--max_multiallelic'])
        self.block_max_missing = int(args['--max_missing'])
        self.overwrite = True if args['--force'] else False

    def _get_block_gap_run(self):
        return (self.block_span - self.block_length - 1)
    
    def _get_block_span(self, block_span):
        if block_span is None:
            return 2*self.block_length
        else:
            block_span = self._get_int(block_span)
            if block_span >= self.block_length:
                return block_span
            else:
                sys.exit("[X] Block span ('-m') must be greater than block length '-l'")

def main(params):
    try:
        start_time = timer()
        print("[+] Running 'gimble blocks'")
        args = docopt(__doc__)
        parameterObj = ParameterObj(params, args)
        gimbleStore = lib.gimble.Store(path=parameterObj.zstore, create=False)
        gimbleStore.blocks(parameterObj)
        gimbleStore.dump_blocks(parameterObj)
        gimbleStore.info()
        print("[*] Total runtime: %.3fs" % (timer() - start_time))
    except KeyboardInterrupt:
        print("\n[X] Interrupted by user after %s seconds!\n" % (timer() - start_time))
        exit(-1)