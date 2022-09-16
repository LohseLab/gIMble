"""usage: gimbl query                    -z <DIR> [-l <STR>] [--sliced STR] [--constrain STR] [-h|--help]
                                            
        -z, --zarr_f DIR                 ZARR datastore
        -l, --label <STR>                Data label
        
        Gridsearch results
        --sliced STR                    Write "sliced_parameter" table with best lnCL result for each slice of parameter values, e.g.:
                                            '--sliced me'
                                            '--sliced Ne_A'
                                            ... 
                                        MUST match parameters in the grid
        --constrain STR                 Write table with best lnCL result for constrained parameter values, e.g.: 
                                            '--constrained me=0.0'
                                            '--constrained me=0.0,Ne_A=100000'
                                            '--constrained me=0.0,Ne_A=100000,Ne_B=200000'
                                            ...
                                            MUST match parameters and values of the grid
        -h --help                        show this
"""
from timeit import default_timer as timer
from docopt import docopt
import lib.gimble
import sys

class QueryParameterObj(lib.gimble.ParameterObj):
    '''Sanitises command line arguments and stores parameters.'''

    def __init__(self, params, args):
        super().__init__(params)
        self.zstore = self._get_path(args['--zarr_f'])
        self.data_key = args['--label']
        self.extended = False # args['--extended']
        self.sliced_param = args['--sliced']
        self.fixed_param = args['--constrain']
        print(self.fixed_param)
        if self.fixed_param:
            try:
                self.fixed_param = {element.split("=")[0]: float(element.split("=")[1]) for element in self.fixed_param.split(",") if element}
            except (ValueError, IndexError):
                sys.exit("[X] '--fixed' is not in the right format.")

def main(params):
    try:
        start_time = timer()
        args = docopt(__doc__)
        parameterObj = QueryParameterObj(params, args)
        #print("parameterObj", parameterObj.__dict__)
        gimbleStore = lib.gimble.Store(path=parameterObj.zstore)
        gimbleStore.query(
            parameterObj._VERSION,
            parameterObj.data_key,
            parameterObj.extended,
            parameterObj.fixed_param,
            parameterObj.sliced_param
            )
        print("[*] Total runtime was %s" % (lib.gimble.format_time(timer() - start_time)))
    except KeyboardInterrupt:
        print("\n[X] Interrupted by user after %s !\n" % (lib.gimble.format_time(timer() - start_time)))
        exit(-1)