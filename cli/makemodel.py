#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""usage: gimbl makemodel                       [-o <STR>] (-1|-2|-3|-4|-5) [-h|--help]
                                            
    [Options]
        -o, --outprefix <STR>                   Prefix to use for INI config file [default: gimble]

        -1, --DIV                               Divergence model config file. 
                                                    No Migration. Third Ne for ancestral population (Ne_A_B).
        -2, --MIG_AB                            Migration model config file. 
                                                    Migration from A to B, backwards in time. B is source population.
        -3, --MIG_BA                            Migration model config file. 
                                                    Migration from B to A, backwards in time. A is source population.
        -4, --IM_AB                             IM model config file. 
                                                    Migration from A to B. Third Ne for ancestral population (Ne_A_B) 
        -5, --IM_BA                             IM model config file. 
                                                    Migration from B to A. Third Ne for ancestral population (Ne_A_B) 
        -h --help                               show this
        
"""

'''
[To do]
- write documentation (no automated filling in)

    Ne lower bounds = min(Pi_A / 4 mu, Pi_B / 4 mu) / 2
    Ne upper bounds = d_xy / (4 mu)

    T lower bound = 0
    T upper bound = dxy / (2 mu)

'''

from docopt import docopt
from timeit import default_timer as timer
import lib.gimble

MODELS = set(['--DIV', '--MIG_AB', '--MIG_BA', '--IM_AB', '--IM_BA'])

class ModelParameterObj(lib.gimble.ParameterObj):
    '''Sanitises command line arguments and stores parameters'''

    def __init__(self, params, args):
        super().__init__(params)
        self.model = self._get_model(args)
        self.outfile = self._get_outfile(args['--outprefix'])

    def _get_model(self, args):
        return [arg for arg, boolean in args.items() if boolean and arg in MODELS][0].replace('--', '')

    def _get_outfile(self, outprefix):
        return ".".join([outprefix, self.model, 'ini'])

def main(params):
    try:
        start_time = timer()
        args = docopt(__doc__)
        parameterObj = ModelParameterObj(params, args)
        lib.gimble.write_config(
            version=parameterObj._VERSION,
            model=parameterObj.model,
            outfile=parameterObj.outfile
            )
        print("[*] Total runtime was %s" % (lib.gimble.format_time(timer() - start_time)))
    except KeyboardInterrupt:
        print("\n[X] Interrupted by user after %s !\n" % (lib.gimble.format_time(timer() - start_time)))
        exit(-1)