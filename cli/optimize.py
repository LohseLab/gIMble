#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
usage: gimble optimize                  [-z FILE] -c FILE (-b [-n INT]|-w| --simID <STR>) [--track] 
                                            [-x FLOAT -i INT] [-f FLOAT] [-p] [--numc <INT>]
                                            [-h|--help]
                                            
                                            
    Options:
        -h --help                                   show this
        -z, --zarr_file FILE                        ZARR datastore
        -c, --config_file FILE                      INI config file
        -b, --blocks                                Optimize based on blocks 
        -w, --windows                               Optimize based on windows (might take very long)
        --simID STR                                 Provide name of simulation run to optimize
        --numc INT                                  Number of cores available [default: 1] 
        -n, --n_points INT                          Number of starting points [default: 1]
        -i, --iterations INT                        Number of iterations to perform when optimizing [default: 100]
        -x, --xtol_rel FLOAT                        Set relative tolerance on norm of vector of optimisation parameters, float between 0 and 1 [default: -1.0]
        -f, --ftol_rel FLOAT                        Set relative tolerance on lnCL, float between 0 and 1 [default: -1.0]
        --track                                     Track likelihood search                        
"""
from timeit import default_timer as timer
from docopt import docopt
import sys
import lib.gimble

'''
[To Do]
- change --simID to --data_id
    - blocks, windows, windowsblocks ?
- document how tols are used

- needs figuring out NLOPT
 - 1%, 0.1%, 0.01%, etc
'''
class OptimizeParameterObj(lib.gimble.ParameterObj):
    '''Sanitises command line arguments and stores parameters.'''

    def __init__(self, params, args):
        super().__init__(params)
        self.zstore = self._get_path(args['--zarr_file'])
        self.data_type = self._get_datatype(args)
        self.config_file = self._get_path(args['--config_file'])
        self.old_parse_config(self.config_file)
        #self._parse_config(self.config_file)
        self.model_file = self._get_model_f()
        self.gridThreads = self._get_int(args['--numc']) #number of workers for independent processes
        self.numPoints = self._get_int(args['--n_points'])
        self.max_eval = self._get_int(args['--iterations'])
        self.xtol_rel = self._get_float(args['--xtol_rel'])
        self.ftol_rel = self._get_float(args['--ftol_rel'])
        self.trackHistory = args['--track']

    def _get_datatype(self, args):
        choices = [args['--blocks'], args['--windows'], args['--simID']]
        if all(choices) or not any(choices):
            sys.exit("[X] Please specify either '--blocks', '--windows' or '--simID'.")
        if args['--blocks']:
            return 'blocks'
        if args['--windows']:
            return 'windows'
        if args['--simID']:
            self.label = args["--simID"]
            return 'simulate'

    def _get_model_f(self):
        return self.config['gimble']['model']

def main(params):
    try:
        start_time = timer()
        args = docopt(__doc__)
        parameterObj = OptimizeParameterObj(params, args)
        gimbleStore = lib.gimble.Store(path=parameterObj.zstore, create=False)
        gimbleStore.optimize(parameterObj)
        print("[*] Total runtime was %s" % (lib.gimble.format_time(timer() - start_time)))
    except KeyboardInterrupt:
        print("\n[X] Interrupted by user after %s !\n" % (lib.gimble.format_time(timer() - start_time)))
        exit(-1)