"""
usage: gimble delete                  -z <z> -d <d> [-h]

    [Options]
        -z, --zarr=<z>                   Prefix to use for GimbleStore
        -d, --data_key=<d>               Data key to delete 
        -h, --help                       Show this
    
"""

import sys
from timeit import default_timer as timer
from docopt import docopt
import lib.runargs 

class DeleteParameterObj(lib.runargs.RunArgs):
    '''Sanitises command line arguments and stores parameters.'''

    def __init__(self, params, args):
        super().__init__(params)
        self.zstore = self._get_path(args["--zarr"])
        self.data_key = args["--data_key"]

def main(params):
    try:
        start_time = timer()
        args = docopt(__doc__)
        parameterObj = DeleteParameterObj(params, args)
        import lib.gimble
        gimbleStore = lib.gimble.Store(path=parameterObj.zstore)
        gimbleStore.delete_key(data_key=parameterObj.data_key)
        print("[*] Total runtime was %s" % (lib.runargs.format_time(timer() - start_time)))
    except KeyboardInterrupt:
        print("\n[X] Interrupted by user after %s !\n" % (lib.runargs.format_time(timer() - start_time)))
        exit(-1)