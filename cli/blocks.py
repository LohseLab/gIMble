"""usage: gIMble blocks         -z <DIR> [-l <INT> -m <INT> -r <INT> -u <INT> -i <INT> -h -d]
    
    -z, --zarr <DIR>                            gIMble ZARR directory
    
    -l, --block_length <INT>                    Successively genotyped sites per block [default: 64] 
    -m, --block_span <INT>                      Maximum distance between first and last site of a block [default: 80]
    -r, --block_gap_run <INT>                   Maximum number of consecutive gaps within a block [default: 1]
    
    -u, --max_multiallelic <INT>                Max multiallelics per block [default: 2]
    -i, --max_missing <INT>                     Max missing per block [default: 2]
    -d, --debug                                 Write debugging logs
    -h, --help
"""

from timeit import default_timer as timer
from docopt import docopt
#import lib.gimblelog
import lib.gimble

class ParameterObj(object):
    def __init__(self, args):
        print(args)
        self.zstore = args['--zarr']
        self.block_length = int(args['--block_length'])
        self.block_span = int(args['--block_span'])
        self.block_gap_run = int(args['--block_gap_run'])
        self.max_multiallelic = int(args['--max_multiallelic'])
        self.max_missing = int(args['--max_missing'])

def main(run_params):
    try:
        start_time = timer()
        args = docopt(__doc__)
        print(args)
        #log = lib.log.get_logger(run_params)
        parameterObj = ParameterObj(args)
        store = lib.gimble.load_store(parameterObj)
        #print(store, type(store) )
        store.make_blocks(parameterObj)
        store.dump_blocks(parameterObj)
        #print(store.tree())
        #print(store.attrs())
        print("[*] Total runtime: %.3fs" % (timer() - start_time))
    except KeyboardInterrupt:
        print("\n[X] Interrupted by user after %s seconds!\n" % (timer() - start_time))
        exit(-1)