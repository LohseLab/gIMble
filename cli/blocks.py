"""usage: gIMble blocks         -z <DIR> [-l <INT> -m <INT> -r <INT> -u <INT> -i <INT> -d -f -h]
    
    -z, --zarr <DIR>                            gIMble ZARR directory
    
    -l, --block_length <INT>                    Successively genotyped sites per block [default: 64] 
    -m, --block_span <INT>                      Maximum distance between first and last site of a block (default: '-l' * 2)
    -r, --block_gap_run <INT>                   Maximum number of consecutive gaps within a block (default: '-l')
    -u, --max_multiallelic <INT>                Max multiallelics per block [default: 2]
    -i, --max_missing <INT>                     Max missing per block [default: 2]

    -f, --force                                 Force overwrite of existing data
    -d, --debug                                 Write debugging logs
    -h, --help
"""

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
        self.block_length = int(args['--block_length'])
        self.block_span = int(args['--block_span']) if not args['--block_span'] is None else 2*self.block_length
        self.block_gap_run = int(args['--block_gap_run']) if not args['--block_gap_run'] is None else self.block_length
        self.block_max_multiallelic = int(args['--max_multiallelic'])
        self.block_max_missing = int(args['--max_missing'])
        self.overwrite = True if args['--force'] else False

def main(params):
    try:
        start_time = timer()
        args = docopt(__doc__)
        #print(args)
        #log = lib.log.get_logger(run_params)
        parameterObj = ParameterObj(params, args)
        store = lib.gimble.load_store(parameterObj)
        #print(store, type(store) )
        store.make_blocks(parameterObj)
        store.dump_blocks(parameterObj)
        #print(store.tree())
        print("[*] Total runtime: %.3fs" % (timer() - start_time))
    except KeyboardInterrupt:
        print("\n[X] Interrupted by user after %s seconds!\n" % (timer() - start_time))
        exit(-1)