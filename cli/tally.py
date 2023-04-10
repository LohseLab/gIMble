"""
usage: gimble tally                     -z <z> -t <t> -l <l> [-k <k>] [-s <s>] [-f] [-h]
                                           
    [Options]
        -z, --zarr=<z>                  Path to existing GimbleStore
        -t, --data_type=<t>             Type of data to tally
                                            - 'blocks': inter-population (X) blocks
                                            - 'windows': windows of inter-population (X) blocks
        -l, --data_label=<l>            Label under which the tally gets saved
        -k, --kmax=<k>                  Max count per mutation type beyond which counts 
                                            are treated as marginals. Order of mutation 
                                            types is (hetB, hetA, hetAB, fixed)
                                            [default: 2,2,2,2]
        -s, --sequence_ids=<s>          Sequence IDs for which to tally blocks (comma-separated)
        -f, --overwrite                 Overwrite results in GimbleStore
        -h, --help                      Show this
"""

from timeit import default_timer as timer
from docopt import docopt
import lib.runargs
import sys

class TallyParameterObj(lib.runargs.RunArgs):
    '''Sanitises command line arguments and stores parameters.'''

    def __init__(self, params, args):
        super().__init__(params)
        self.zstore = self._get_path(args['--zarr'])
        self.data_type = self._get_data_type(args['--data_type'])
        self.data_label = self._check_label(args['--data_label'])
        self.max_k = self._check_kmax(args['--kmax'])
        self.overwrite = args['--overwrite']
        self.sample_sets = 'X'
        self.sequence_ids = args['--sequence_ids']
        self.genome_file = None
    
    def _check_label(self, label):
        invalid_chars = set([c for c in label if not c.isalnum() and not c in set([".", "-", "_"])])
        if invalid_chars:
            sys.exit("[X] --makegrid_label contains invalid characters (%r). Should only contain alphanumericals and -_." % "".join(invalid_chars))
        return label

    def _get_data_type(self, data_type):
        error = "[X] '--data_type' for tally must be 'blocks' or 'windows'. Not %r." % data_type
        if not data_type in set(['blocks', 'windows']):
            sys.exit(error)
        return data_type

def main(params):
    try:
        start_time = timer()
        args = docopt(__doc__)
        print("[+] Running 'gimble tally' ...")
        parameterObj = TallyParameterObj(params, args)
        import lib.gimble
        gimbleStore = lib.gimble.Store(path=parameterObj.zstore, create=False)
        tally_key = gimbleStore.tally(
            data_type=parameterObj.data_type,
            data_label=parameterObj.data_label,
            max_k=parameterObj.max_k,
            sample_sets=parameterObj.sample_sets,
            sequence_ids=parameterObj.sequence_ids,
            genome_file=parameterObj.genome_file,
            overwrite=parameterObj.overwrite
            )
        print("[+] Tally is accessible with the key %r." % tally_key)
        print("[*] Total runtime was %s" % (lib.runargs.format_time(timer() - start_time)))
    except KeyboardInterrupt:
        print("\n[X] Interrupted by user after %s !\n" % (lib.runargs.format_time(timer() - start_time)))
        exit(-1)