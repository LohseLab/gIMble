#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""usage: gIMble makegrid -c <FILE> -m <FILE> -z <FILE> [-b <INT>] [-t <INT>] [-h|--help]
                                            
    Options:
        -h --help                                show this
        -c, --config_file <FILE>                 Config file with parameters
        -m, --model_file <FILE>                  Model file
        -z, --zarr <FILE>                        Path to zarr store
        -b, --blocklength <INT>                  Blocklength
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
        self.zarr_not_existing, self.zstore = self._get_path(args["--zarr"])
        self.block_length= self._get_int(args['--blocklength'], ret_none=True)
        if not (self.block_length or self.zstore):
            sys.exit('[X] Provide either blocklength or zarr store with block length in attrs.')
        _, self.config_file = self._get_path(args['--config_file'], doesNotExistError=True)
        _, self.model_file = self._get_path(args['--model_file'], doesNotExistError=True)
        #self.grid_path = self._verify_parent(args['--grid_file'])
        self.threads = self._get_int(args["--threads"])
        self.config = self._parse_config(self.config_file)
        self._process_config()

    def _get_path(self, infile, doesNotExistError=False):
        if infile is None:
            return None
        path = pathlib.Path(infile).resolve()
        if not path.exists():
            if doesNotExistError:
                sys.exit("[X] File not found: %r" % str(infile))
            else:
                parent_path=pathlib.Path(path.parent).resolve()
                if not parent_path.exists():
                    sys.exit("[X] Parent directory not found: %r" % str(parent_path))
                else:
                    print("[+] Creating new zarr store.")
                    return (True, str(path))
        return (False, str(path))

def main(params):
    try:
        start_time = timer()
        args = docopt(__doc__)
        parameterObj = MakeGridParameterObj(params, args)
        print("[+] Generated all parameter combinations.") #in parameterObj.grid
        
        equationSystem = lib.math.EquationSystemObj(parameterObj)
        #build the equations
        equationSystem.initiate_model()
        equationSystem.calculate_all_ETPs()
        
        #for zarr store: require grid
        gimbleStore = lib.gimble.Store(path=parameterObj.zstore, create=parameterObj.zarr_not_existing)
        gimbleStore.data.require_group('grids')
        if parameterObj.zarr_not_existing:
            gridcount = 0
        else:
            gridcount = gimbleStore._return_group_last_integer('grids')
        paramdict = {i:paramset for i, paramset in enumerate(parameterObj.grid)}
        g = gimbleStore.data['grids'].create_dataset(f'grid_{gridcount}', data=equationSystem.ETPs)
        g.attrs.put(paramdict)

        print("[*] Total runtime: %.3fs" % (timer() - start_time))
    except KeyboardInterrupt:
        print("\n[X] Interrupted by user after %s seconds!\n" % (timer() - start_time))
        exit(-1)