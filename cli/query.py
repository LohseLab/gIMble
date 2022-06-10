"""usage: gimbl query                    -z <DIR> [-l <STR>] [--extended] [-h|--help] [--demes]
                                            
        -z, --zarr_f DIR                 ZARR datastore
        -l, --label <STR>                Data label
        --demes                          Pipe winning model to demes.yml (where applicable)
                                            - gridsearch
                                            - optimize
        --extended                       Write "extended" format (where applicable)
                                            - gridsearch
        -h --help                        show this
"""
from timeit import default_timer as timer
from docopt import docopt
import lib.gimble

class QueryParameterObj(lib.gimble.ParameterObj):
    '''Sanitises command line arguments and stores parameters.'''

    def __init__(self, params, args):
        super().__init__(params)
        self.zstore = self._get_path(args['--zarr_f'])
        self.data_key = args['--label']
        self.extended = args['--extended']
        self.demes = args['--demes']

def main(params):
    try:
        start_time = timer()
        args = docopt(__doc__)
        parameterObj = QueryParameterObj(params, args)
        gimbleStore = lib.gimble.Store(path=parameterObj.zstore)
        gimbleStore.query(
            parameterObj._VERSION,
            parameterObj.data_key,
            parameterObj.extended,
            parameterObj.demes
            )
        print("[*] Total runtime was %s" % (lib.gimble.format_time(timer() - start_time)))
    except KeyboardInterrupt:
        print("\n[X] Interrupted by user after %s !\n" % (lib.gimble.format_time(timer() - start_time)))
        exit(-1)