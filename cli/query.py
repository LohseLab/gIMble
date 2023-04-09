"""
usage: gimble query                     -z <z> [-d <d>] [-s <s>] [-c <c>] [-h|--help]
                                            
        -z, --zarr_file=<z>             Path to existing GimbleStore 
        -d, --data_key=<d>              Data key
        
    [Gridsearch results]
        -s, --sliced=<s>                Write "sliced_parameter" table with best lnCL result for each slice of parameter values, e.g.:
                                            '--sliced me'
                                            '--sliced Ne_A'
                                            ... 
                                            MUST match parameters in the grid
        -c, --constrained=<c>           Write table with best lnCL result for constrained parameter values, e.g.: 
                                            '--constrained me=0.0'
                                            '--constrained me=0.0,Ne_A=100000'
                                            '--constrained me=0.0,Ne_A=100000,Ne_B=200000'
                                            ...
                                            MUST match parameters and values of the grid
        -h,--help                       Show this
"""
from timeit import default_timer as timer
from docopt import docopt
import lib.gimble
import sys

class QueryParameterObj(lib.gimble.ParameterObj):
    '''Sanitises command line arguments and stores parameters.'''

    def __init__(self, params, args):
        super().__init__(params)
        self.zstore = self._get_path(args['--zarr_file'])
        self.data_key = args['--data_key']
        self.extended = False # args['--extended']
        self.sliced_param = args['--sliced']
        self.constrained = args['--constrained']
        self.diss = False
        if self.constrained:
            try:
                self.constrained = {element.split("=")[0]: float(element.split("=")[1]) for element in self.constrained.split(",") if element}
            except (ValueError, IndexError):
                sys.exit("[X] '--constrained' is not in the right format.")

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
            parameterObj.constrained,
            parameterObj.sliced_param,
            parameterObj.diss
            )
        print("[*] Total runtime was %s" % (lib.gimble.format_time(timer() - start_time)))
    except KeyboardInterrupt:
        print("\n[X] Interrupted by user after %s !\n" % (lib.gimble.format_time(timer() - start_time)))
        exit(-1)