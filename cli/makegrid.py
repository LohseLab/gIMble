#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""usage: gIMble makegrid -c <FILE> -m <FILE> [-z <FILE>] [-o <STR>] [-t <INT>] [-h|--help]
                                            
    Options:
        -h --help                                show this
        -c, --config_file <FILE>                 Config file with parameters
        -m, --model_file <FILE>                  Model file
        -z, --zarr <FILE>                        Path to zarr store
        -o, --outprefix <STR>                    Prefix to use for gimble store [default: gimble]
        -t, --threads <INT>                      Threads [default: 1]
        
"""
import pathlib
import collections
from timeit import default_timer as timer
from docopt import docopt
import sys
import lib.gimble
import lib.math

'''

Grid:
- which gridpoint gives highes Likelihood
    - overall
    - conditional, when Me is set to 0

# bootstrapping
after running simulations
- for each grid points, do simulation replicates
- somehow use rembination rate

gimble scan (combininig windows data with likelihoods)
        
grid:=
    center of grid : 
        - "global" model estimated from overal block_counts across genome in canonical workflow
        - Ne_A, Ne_B, M_e
    boundaries: 
        - all required! 
        - have to be checked for consistency! (logarithmic gridding?)
        Ne_A_min, Ne_A_max, Ne_A_steps (min=0 makes no sense, though)
        Ne_B_min, Ne_B_max, Ne_B_steps (min=0 makes no sense, though)
        M_e_min, M_e_Max, M_e_steps (def want to include min=0)
             
                          centre
step                        *
|   |   |   |   |   |   |   |   |   |   |   |   |   |   |        
min                                                     max

l = range(min, max+step, step)  # should allow for linear/log spacing                                                   
if not centre in l:  # identity check has to allow for precision-offness
    error
if M_e and not 0 in l: # identity check has to allow for precision-offness
    error   

'''

class MakeGridParameterObj(lib.gimble.ParameterObj):
    '''Sanitises command line arguments and stores parameters.'''

    def __init__(self, params, args):
        super().__init__(params)
        self.zstore = self._get_path(args["--zarr"])
        self.prefix = self._get_prefix(args["--outprefix"])
        self.config_file = self._get_path(args['--config_file'])
        self.model_file = self._get_path(args['--model_file'])
        self.threads = self._get_int(args["--threads"])
        self.config = self._parse_config(self.config_file)
        self._process_config()

def main(params):
    try:
        start_time = timer()
        args = docopt(__doc__)
        parameterObj = MakeGridParameterObj(params, args)
        if parameterObj.zstore:
            gimbleStore = lib.gimble.Store(path=parameterObj.zstore)
        elif parameterObj.prefix:
            gimbleStore = lib.gimble.Store(prefix=parameterObj.prefix, create=True)
        else:
            sys.exit("[X] No config and no prefix specified. Should have been caught.")
        
        print("[+] Generated all parameter combinations.") #in parameterObj.parameter_combinations
        
        equationSystem = lib.math.EquationSystemObj(parameterObj)
        
        #print(f'rates by variable: {equationSystem.rate_by_variable}')
        #print(f'split times:{equationSystem.split_times}')
        
        #build the equations
        equationSystem.initiate_model()
        equationSystem.ETPs = equationSystem.calculate_all_ETPs()
        
        run_count = gimbleStore._return_group_last_integer('grids')
        g = gimbleStore.data['grids'].create_dataset(f'grid_{run_count}', data=equationSystem.ETPs)
        g.attrs.put({idx:combo for idx, combo in enumerate(parameterObj.parameter_combinations)})
        print("[*] Total runtime: %.3fs" % (timer() - start_time))
    except KeyboardInterrupt:
        print("\n[X] Interrupted by user after %s seconds!\n" % (timer() - start_time))
        exit(-1)