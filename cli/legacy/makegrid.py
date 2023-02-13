"""usage: gimbl makegrid                 (-z <FILE> | -o <STR>) -c <FILE> [-f] [-n <INT> ] [-h|--help]
                                            
    Options:
        -h --help                        show this

        -z, --zarr_file <FILE>           Path to existing GimbleStore
        -o, --outprefix <STR>            Prefix to use for new GimbleStore
        -c, --config_file <FILE>         Config file with model parameters
        -f, --force                      Force overwrite of grid in GimbleStore
        -n, --processes <INT>            Number of processes [default: 1] 
        
"""

'''
- Do we really want to allow entrypoint of workflow to be makegrid?
    + as opposed to only simulate and parse?
    + check for product of block_length and mutation rate
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
        self.num_cores = self._get_int(args['--processes']) #number of workers for independent processes
        self.overwrite = args['--force']
        self.config = lib.gimble.load_config(
            self.config_file, 
            self._MODULE, 
            self._CWD, 
            self._VERSION)

def main(params):
    try:
        start_time = timer()
        args = docopt(__doc__)
        parameterObj = MakeGridParameterObj(params, args)
        gimbleStore = lib.gimble.Store(
            path=parameterObj.zstore, 
            prefix=parameterObj.prefix, 
            create=(False if parameterObj.zstore else True))
        gimbleStore.makegrid(
            config=parameterObj.config,
            num_cores=parameterObj.num_cores,
            overwrite=parameterObj.overwrite,
            )
        print("[*] Total runtime was %s" % (lib.gimble.format_time(timer() - start_time)))
    except KeyboardInterrupt:
        print("\n[X] Interrupted by user after %s !\n" % (lib.gimble.format_time(timer() - start_time)))
        exit(-1)