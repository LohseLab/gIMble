"""usage: gimbl query                    -z <DIR> [-l <STR>] [--extended] [--fixed_param STR] [-h|--help]
                                            
        -z, --zarr_f DIR                 ZARR datastore
        -l, --label <STR>                Data label
        
        Gridsearch results
        --extended                       Write "extended" table (all parameters!)
        --fixed_param STR                  Write "fixed_parameter" table

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
        self.fixed_param = args['--fixed_param']

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
            parameterObj.fixed_param
            )
        print("[*] Total runtime was %s" % (lib.gimble.format_time(timer() - start_time)))
    except KeyboardInterrupt:
        print("\n[X] Interrupted by user after %s !\n" % (lib.gimble.format_time(timer() - start_time)))
        exit(-1)