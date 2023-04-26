"""
usage: gimble list                  -z <z> [-m <m>] [-h]

    [Options]
        -z, --zarr=<z>                   Prefix to use for GimbleStore
        -m, --module=<m>                 List items in GimbleStore created by a specific module [default: all]
        -h, --help                       Show this
    
"""

'''
[ToDo]
- timestamps
- size in Gb
'''
import sys
from timeit import default_timer as timer
from docopt import docopt
import lib.runargs 

PATH_BY_MODULES = {
    'blocks': 'blocks/',       
    'windows': 'windows/',      
    'tally': 'tally',
    'simulate': 'simulate/',
    'optimize': 'optimize/',
    'makegrid': 'makegrid/',
    'gridsearch': 'gridsearch/',
}

class ListParameterObj(lib.runargs.RunArgs):
    '''Sanitises command line arguments and stores parameters.'''

    def __init__(self, params, args):
        super().__init__(params)
        self.zstore = self._get_path(args["--zarr"])
        self.module = self._check_module(args["--module"])

    def _check_module(self, module):
        if module == 'all':
            return None
        if not module in PATH_BY_MODULES:
            sys.exit("[X] Module must be one of the following: all, %s" % (", ".join(PATH_BY_MODULES.keys())))
        return module


def main(params):
    try:
        start_time = timer()
        args = docopt(__doc__)
        parameterObj = ListParameterObj(params, args)
        import lib.gimble
        gimbleStore = lib.gimble.Store(path=parameterObj.zstore)
        gimbleStore.list_keys(category=parameterObj.module)
        print("[*] Total runtime was %s" % (lib.runargs.format_time(timer() - start_time)))
    except KeyboardInterrupt:
        print("\n[X] Interrupted by user after %s !\n" % (lib.runargs.format_time(timer() - start_time)))
        exit(-1)