"""usage: gimbl tally                       -z <DIR> -d <STR> -l <STR> [-k <STR>] 
                                            [-s <STR>] [-f] [-h]
                                            

        -z, --zarr <DIR>                     Path to existing GimbleStore
        -d, --data_source <STR>              'blocks', 'windows', or 'windowsum'
        -l, --data_label <STR>               Label under which the tally gets saved
    
    [Options]
        -k, --maxk <STR>                            Max value for mutypes (values above get binned) [default: 2,2,2,2]
        -s, --sequence_ids <STR>                    Sequence IDs for which to tally blocks (comma-separated)
        -f, --overwrite                             Overwrite results in GimbleStore

        -h --help                                   show this
"""

from timeit import default_timer as timer
from docopt import docopt
import lib.gimble
import lib.math

class TallyParameterObj(lib.gimble.ParameterObj):
    '''Sanitises command line arguments and stores parameters.'''

    def __init__(self, params, args):
        super().__init__(params)
        self.zstore = self._get_path(args['--zarr'])
        self.data_source = self._get_data_source(args['--data_source'])
        self.data_label = args['--data_label']
        self.max_k = self._get_max_k(args['--maxk'])
        self.overwrite = args['--overwrite']
        self.sample_sets = 'X'
        self.sequence_ids = args['--sequence_ids']
        self.genome_file = None
    
    def _get_data_source(self, data_source):
        error = "[X] '--data_source' for tally must be 'blocks', 'windows', or 'windowsum'. Not %r" % data_source
        if not data_source in set(['blocks', 'windows', 'windowsum']):
            return error
        return data_source

def main(params):
    try:
        start_time = timer()
        args = docopt(__doc__)
        parameterObj = TallyParameterObj(params, args)
        gimbleStore = lib.gimble.Store(path=parameterObj.zstore, create=False)
        tally_key = gimbleStore.tally(
            data_source=parameterObj.data_source,
            data_label=parameterObj.data_label,
            max_k=parameterObj.max_k,
            sample_sets=parameterObj.sample_sets,
            sequence_ids=parameterObj.sequence_ids,
            genome_file=parameterObj.genome_file,
            overwrite=parameterObj.overwrite
            )
        print("[+] Tally is accessible with the key %r (use this for query)." % tally_key)
        print("[*] Total runtime was %s" % (lib.gimble.format_time(timer() - start_time)))
    except KeyboardInterrupt:
        print("\n[X] Interrupted by user after %s !\n" % (lib.gimble.format_time(timer() - start_time)))
        exit(-1)