#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
usage: gimble optimize                  -z FILE -c FILE (--blocks [-p INT] | -s LABEL | --windows) 
                                        [-x FLOAT -f FLOAT -i INT] [-n INT]
                                        [--track] [-h|--help]

    Data
        -z, --zarr_file FILE                        ZARR datastore
        -b, --blocks                                Optimize based on blocks 
        -w, --windows                               Optimize based on windows (TBD)
        -s, --simulations LABEL                     Optimize based on simulations

    Stopping criteria of optimization
        -i, --max_iterations INT                    Maximum number of iterations to perform when optimizing [default: 100]
        -x, --xtol_rel FLOAT                        Relative tolerance on norm of vector of optimisation parameters
                                                        Float between 0 and 1, deactivate with -1 [default: -1.0]
        -f, --ftol_rel FLOAT                        Relative tolerance on lnCL 
                                                        Float between 0 and 1, deactivate with -1 [default: -1.0]
    Options
        -c, --config_file FILE                      INI config file
        -n, --num_cores INT                         Number of cores [default: 1] 
        -p, --start_points INT                      Number of starting points. (One is fixed, n-1 random) [default: 1]
        -t, --track                                 Track likelihood search (not available for simulations)   
        -h --help                                   show this                    
    
"""
from timeit import default_timer as timer
from docopt import docopt
import sys
import lib.gimble

'''
[To Do]
- --start_points has to be removed
- 

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
        self.num_cores = self._get_int(args['--num_cores'])    # number of workers for independent processes
        self.start_points = self._get_int(args['--start_points'])
        self.max_iterations = self._get_int(args['--max_iterations'])
        self.xtol_rel = self._get_float(args['--xtol_rel'])
        self.ftol_rel = self._get_float(args['--ftol_rel'])
        self.track_history = args['--track']
        self.simulations_label = None # filled in by OptimizeParameterObj._get_datatype()
        self.config = lib.gimble.load_config(self.config_file, self._MODULE, self._CWD, self._VERSION)

    def _get_datatype(self, args):
        choices = [args['--blocks'], args['--windows'], args['--simulations']]
        if all(choices) or not any(choices):
            sys.exit("[X] Please specify either '--blocks', '--windows' or '--simulations'.")
        if args['--blocks']:
            return 'blocks'
        if args['--windows']:
            return 'windows'
        if args['--simulations']:
            self.simulations_label = args["--simulations"]
            return 'simulations'

def main(params):
    try:
        start_time = timer()
        args = docopt(__doc__)
        parameterObj = OptimizeParameterObj(params, args)
        gimbleStore = lib.gimble.Store(path=parameterObj.zstore, create=False)
        #gimbleStore.optimize_before(parameterObj)
        gimbleStore.optimize(
            data_type=parameterObj.data_type,
            config=parameterObj.config,
            num_cores=parameterObj.num_cores,
            start_points=parameterObj.start_points,
            max_iterations=parameterObj.max_iterations,
            xtol_rel=parameterObj.xtol_rel,
            ftol_rel=parameterObj.ftol_rel,
            track_history=parameterObj.track_history,
            simulations_label=parameterObj.simulations_label)
        print("[*] Total runtime was %s" % (lib.gimble.format_time(timer() - start_time)))
    except KeyboardInterrupt:
        print("\n[X] Interrupted by user after %s !\n" % (lib.gimble.format_time(timer() - start_time)))
        exit(-1)