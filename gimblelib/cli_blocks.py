"""usage: gimbl blocks         -z DIR [-l INT -m INT -u INT -i INT -d -f -h]
    
    -z, --zarr DIR                            GimbleStore (*.z) in which to generate blocks
    -l, --block_length INT                    Number of sites per block [default: 64] 
    -m, --block_span INT                      Maximum distance between first and last site in a block (default: '-l' * 2)
    -u, --max_multiallelic INT                Max number of multiallelic genotypes per block (default: round('-l' * 0.05))
    -i, --max_missing INT                     Max number of missing genotypes per block (default: round('-l' * 0.05))
    -f, --force                                 Force overwrite of existing blocks.
    -h, --help                                  Show this screen.
"""

import sys
from timeit import default_timer as timer
from docopt import docopt
import gimble.parameters
import gimble.gimble

class BlockParameterObj(gimble.parameters.ParameterObj):
    '''Sanitises command line arguments and stores parameters'''

    def __init__(self, params, args):
        super().__init__(params)
        self.zstore = self._get_path(args['--zarr'])
        self.block_length = self._get_int(args['--block_length'])
        self.block_span = self._get_block_span(args['--block_span'])
        self.block_gap_run = self._get_block_gap_run()
        self.block_max_multiallelic = self._get_max_values(args['--max_multiallelic'])
        self.block_max_missing = self._get_max_values(args['--max_missing'])
        self.overwrite = True if args['--force'] else False

    def _get_block_span(self, block_span):
        if block_span is None:
            return 2 * self.block_length
        block_span = self._get_int(block_span)
        if block_span < self.block_length:
            sys.exit("[X] Block span ('-m') must be greater or equal to block length '-l'")
        return block_span

    def _get_block_gap_run(self):
        return (self.block_span - self.block_length - 1)
    
    def _get_max_values(self, max_value):
        if max_value is None:
            return round(self.block_length * 0.05)
        return self._get_int(max_value)    
        
def main(params):
    try:
        start_time = timer()
        args = docopt(__doc__)
        print(args)
        parameterObj = BlockParameterObj(params, args)
        print("[+] Running 'gimble blocks'")
        print("[+] Parameters = [-l %s -m %s -u %s -i %s]" % (
            parameterObj.block_length, 
            parameterObj.block_span, 
            parameterObj.block_max_multiallelic, 
            parameterObj.block_max_missing))
        gimbleStore = gimble.gimble.Store(path=parameterObj.zstore)
        gimbleStore.blocks(parameterObj)
        gimbleStore.info()
        print("[*] Total runtime: %.3fs" % (timer() - start_time))
    except KeyboardInterrupt:
        print("\n[X] Interrupted by user after %s seconds!\n" % (timer() - start_time))
        exit(-1)