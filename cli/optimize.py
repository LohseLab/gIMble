#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
usage: gimble optimize                  -z FILE -c FILE (-t <STR>|-s <STR>) [-w]
                                        [--xtol FLOAT --ftol FLOAT -i INT] [-n INT] [-p STR] [--sage]
                                        [-f] [-h|--help]

    -z, --zarr_file FILE                            ZARR Datastore
    -c, --config_file FILE                          INI config file
    -t, --tally_label STR                           Label of tally in ZARR store ('tally/...')
    -s, --sim_label STR                             Label of simulation run in ZARR store ('simulate/...')
    -w, --windowsum                                 Apply optimize to sum of tallies in parsed/simulated windows data
                                                        [default: False]

    Stopping criteria of optimization
        -i, --max_iterations INT                    Maximum number of iterations to perform when 
                                                        optimizing, deactivate with 0 [default: 1000]
        --xtol FLOAT                                Relative tolerance on norm of vector of optimisation parameters
                                                        Float between 0 and 1, deactivate with -1 [default: -1]
        --ftol FLOAT                                Relative tolerance on lnCL 
                                                        Float between 0 and 1, deactivate with -1 [default: -1]
    Options
        -n, --num_cores INT                         Number of cores [default: 1] 
        -p, --start_point STR                       Point from which to start optimization [default: midpoint]
                                                        - 'midpoint' : midpoint between all boundary values
                                                        - 'random': based on random seed in INI file
        --sage                                      Use GeneratingFunctions via Sage [default: False]                                                        
        -f, --force                                 Force overwrite of existing analysis.
        -h --help                                   show this                    
        
"""
from timeit import default_timer as timer
from docopt import docopt
import lib.gimble
import sys

'''
'''
class OptimizeParameterObj(lib.gimble.ParameterObj):
    '''Sanitises command line arguments and stores parameters.'''

    def __init__(self, params, args):
        super().__init__(params)
        self.zstore = self._get_path(args['--zarr_file'])
        self.sim_label = args['--sim_label']
        self.windowsum = args['--windowsum']
        self.tally_label = args['--tally_label']
        self.config_file = self._get_path(args['--config_file'])

        self.max_iterations = self._get_int(args['--max_iterations'])
        self.xtol_rel = self._get_float(args['--xtol'])
        self.ftol_rel = self._get_float(args['--ftol'])

        self.num_cores = self._get_int(args['--num_cores'])    # number of workers for independent processes
        self.start_point = self._check_start_point(args['--start_point'])
        self.force = args['--force']
        self.sage = args['--sage'] 
        
        self.config = lib.gimble.load_config(self.config_file, self._MODULE, self._CWD, self._VERSION)

    def _check_start_point(self, start_point_string):
        if start_point_string in set(['random', 'midpoint']):
            return start_point_string
        sys.exit("[X] '--start_point' must be 'midpoint' or 'random'")

def main(params):
    try:
        start_time = timer()
        print("[+] Running 'gimble optimize' ...")
        args = docopt(__doc__)
        parameterObj = OptimizeParameterObj(params, args)
        gimbleStore = lib.gimble.Store(path=parameterObj.zstore, create=False)
        #gimbleStore.optimize_before(parameterObj)
        if parameterObj.sage:
            # old implementation
            gimbleStore.optimize_gf(
                config=parameterObj.config,
                sim_label=parameterObj.sim_label,
                windowsum=parameterObj.windowsum,
                tally_label=parameterObj.tally_label,
                num_cores=parameterObj.num_cores,
                start_point=parameterObj.start_point,
                max_iterations=parameterObj.max_iterations,
                xtol_rel=parameterObj.xtol_rel,
                ftol_rel=parameterObj.ftol_rel,
                overwrite=parameterObj.force)
        else:
            gimbleStore.optimize(
                config=parameterObj.config,
                sim_label=parameterObj.sim_label,
                windowsum=parameterObj.windowsum,
                tally_label=parameterObj.tally_label,
                num_cores=parameterObj.num_cores,
                start_point=parameterObj.start_point,
                max_iterations=parameterObj.max_iterations,
                xtol_rel=parameterObj.xtol_rel,
                ftol_rel=parameterObj.ftol_rel,
                overwrite=parameterObj.force)
        print("[*] Total runtime was %s" % (lib.gimble.format_time(timer() - start_time)))
    except KeyboardInterrupt:
        print("\n[X] Interrupted by user after %s !\n" % (lib.gimble.format_time(timer() - start_time)))
        exit(-1)