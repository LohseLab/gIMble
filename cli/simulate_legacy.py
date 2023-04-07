#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""usage: gIMble simulate                       (-z DIR | -o STR) -c FILE [-t INT -f] [-h|--help]
                                            
    Options:
        -h --help                               show this
        -z, --zarr DIR                          Path to zarr store
        -o, --outprefix STR                     Prefix to make new zarr store
        -c, --config_file FILE                  Simulate config file (*.ini) 
        -t, --threads INT                       Threads [default: 1]
        -f, --overwrite                         Overwrite results in GimbleStore
        
"""

from timeit import default_timer as timer
from docopt import docopt
import lib.gimble
import lib.simulate

'''
usage: gimble simulate             (-z DIR | -o STR) [-t INT -f] [-h|--help]
                                            
    Options:
        -h --help                                   show this
        -z, --zarr DIR                              Path to zarr store
        -o, --outprefix STR                         Prefix to make new zarr store
        -c, --config_file FILE                      Simulate config file (*.ini) 
        -t, --threads INT                           Threads [default: 1]
        -f, --overwrite                             Overwrite results in GimbleStore

        [Model based]
        --model=STR
        --Ne_A=FLOAT
        --Ne_B=FLOAT
        --Ne_AB=FLOAT
        --T=FLOAT
        --me=FLOAT
        --mu=FLOAT
        --seed=INT
        
        [Gridsearch based]
        --gridsearch_key 
        --constraint
    
        [Recombination]
        --recombination_rate
        --recombination_map 

        [Data]
        --replicates
        --windows 
        --blocks
        --block_length
        --max_k
        --ploidy
        --samples_A 
        --samples_B
        --discrete_genome
'''                                                    

class SimulateParameterObj(lib.gimble.ParameterObj):
    """Sanitises command line arguments and stores parameters."""

    def __init__(self, params, args):
        super().__init__(params)
        self.config_file = self._get_path(args["--config_file"])
        self.zstore = self._get_path(args["--zarr"])
        self.prefix = self._get_prefix(args["--outprefix"])
        self.threads = self._get_int(args["--threads"])
        #self.sim_grid = args["--window_wise_bootstrap"]
        self.overwrite = args['--overwrite']
        self.config = lib.gimble.load_config(
            self.config_file, 
            self._MODULE, 
            self._CWD, 
            self._VERSION)

def main(params):
    try:
        start_time = timer()
        print("[+] Running 'gimble simulate' ...")
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
