"""usage: gimbl makegrid                 (-z <FILE> | -o <STR>) -m <FILE> -c <FILE> [-f] [--inner_pool <INT> --outer_pool <INT>] [-h|--help]
                                            
    Options:
        -h --help                        show this

        -z, --zarr_file <FILE>           Path to existing GimbleStore
        -o, --outprefix <STR>            Prefix to use for new GimbleStore
        -m, --model_file <FILE>          Model file
        -c, --config_file <FILE>         Config file with model parameters
        -f, --overwrite                  Overwrite grid in GimbleStore
        --inner_pool INT                 Number of processes used to optimize a single data point [default: 1] 
        --outer_pool INT                 Number of data points processed in parallel [default: 1]
        
"""

'''
[To Do]
- remove --model_file
- remove --outprefix
- remove --inner/outer poll
- change to new config parser
- add --grid_id

'''
from timeit import default_timer as timer
from docopt import docopt
import lib.gimble
import lib.math

class MakeGridParameterObj(lib.gimble.ParameterObj):
    '''Sanitises command line arguments and stores parameters.'''

    def __init__(self, params, args):
        super().__init__(params)
        self.zstore = self._get_path(args["--zarr_file"])
        self.prefix = self._get_prefix(args["--outprefix"])
        self.config_file = self._get_path(args['--config_file'])
        self.model_file = self._get_path(args['--model_file'])
        #self.threads, self.gridThreads = [self._get_int(t) for t in args["--threads"].split(',')]
        self.threads = self._get_int(args['--inner_pool']) #number of workers for a single set of equations to be solved
        self.gridThreads = self._get_int(args['--outer_pool']) #number of workers for independent processes
        self.overwrite = args['--overwrite']
        self.config = None
        self.old_parse_config(self.config_file)

def main(params):
    try:
        start_time = timer()
        args = docopt(__doc__)
        parameterObj = MakeGridParameterObj(params, args)
        gimbleStore = lib.gimble.Store(
            path=parameterObj.zstore, 
            prefix=parameterObj.prefix, 
            create=(False if parameterObj.zstore else True))
        gimbleStore.makegrid(parameterObj)
        print("[*] Total runtime was %s" % (lib.gimble.format_time(timer() - start_time)))
    except KeyboardInterrupt:
        print("\n[X] Interrupted by user after %s !\n" % (lib.gimble.format_time(timer() - start_time)))
        exit(-1)