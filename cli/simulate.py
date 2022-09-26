#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""usage: gIMble simulate                   [-z DIR | -o DIR] -c FILE [-t INT]
                                            [(-w STR [--fixed STR])] [-f] [-h|--help]
                                            
    Options:
        -h --help                                   show this
        -z, --zarr DIR                              Path to zarr store
        -o, --outprefix DIR                         Prefix to make new zarr store
        -c, --config_file FILE                      Simulate config file (*.ini) 
        -t, --threads INT                           Threads [default: 1]
        -f, --overwrite                             Overwrite results in GimbleStore
        -w, --window_wise_bootstrap STR             Label of lncls grid to perform window-wise parametric bootstrap on                 
        --fixed STR                                 Parameter to fix to global optimum when performing window-wise parametric bootstrap
"""

from timeit import default_timer as timer
from docopt import docopt
import lib.gimble
import lib.simulate

'''
# why are these needed?
 -g, --grid                                  
 --fixed STR                                 Parameter to fix to global optimum

# simulations_labels is provided in config file under label
# when running optimize/gridsearch, there needs to be a field for simulations_label

_get_sim_grid: 
- this is about the results of a grid search, right?
- it's defined at least twice

'''                                                    

class SimulateParameterObj(lib.gimble.ParameterObj):
    """Sanitises command line arguments and stores parameters."""

    def __init__(self, params, args):
        super().__init__(params)
        self.config_file = self._get_path(args["--config_file"])
        self.zstore = self._get_path(args["--zarr"])
        self.prefix = self._get_prefix(args["--outprefix"])
        self.threads = self._get_int(args["--threads"])
        self.sim_grid = args["--window_wise_bootstrap"]
        self.overwrite = args['--overwrite']
        self.config = lib.gimble.load_config(
            self.config_file, 
            self._MODULE, 
            self._CWD, 
            self._VERSION)

def main(params):
    try:
        start_time = timer()
        args = docopt(__doc__)
        parameterObj = SimulateParameterObj(params, args)
        if parameterObj.zstore:
            gimble_store = lib.gimble.Store(path=parameterObj.zstore)
        elif parameterObj.prefix:
            gimble_store = lib.gimble.Store(prefix=parameterObj.prefix, create=True)
        else:
            sys.exit("[X] No config and no prefix specified. Should have been caught.")
        #perform recmap checks
        gimble_store.simulate(
            config=parameterObj.config,
            threads=parameterObj.threads,
            overwrite=parameterObj.overwrite
            )
        print("[*] Total runtime: %.3fs" % (timer() - start_time))
    except KeyboardInterrupt:
        print("\n[X] Interrupted by user after %s seconds!\n" % (timer() - start_time))
        exit(-1)
