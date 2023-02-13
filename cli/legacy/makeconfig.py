#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""usage: gimbl makeconfig          -l <STR> -t <STR> -m <INT> [-h|--help]
                                            
    [Options]
        -l, --label <STR>           User specified label for analysis.

        -t, --task <STR>            Create config INI for a specific analysis (module)
                                    'optimize': given a model and data (blocks), the module 'optimize' will 
                                                  return the parameter combinations of the model that 
                                                  maximizes the log-composite likelihood. 
                                    'makegrid': given a model, the module 'makegrid' will pre-calculate 
                                                  the probabilities of individual mutation configurations for
                                                  each of the parameter combinations specified in the grid. 
                                                  Required for running the 'gridsearch' module. 
                                    'simulate': given a model and parameters, the module 'simulate' will
                                                  generate data using msprime.
                            
        -m, --model <INT>           1 : DIV, Divergence model config file.   
                                        No Migration. Third Ne for ancestral population (Ne_A_B).
                                    2 : MIG_AB, Migration model config file. 
                                        Migration from A to B, backwards in time. B is source population.
                                    3 : MIG_BA, Migration model config file. 
                                        Migration from B to A, backwards in time. A is source population.
                                    4 : IM_AB, IM model config file. 
                                        Migration from A to B. Third Ne for ancestral population (Ne_A_B) 
                                    5 : IM_BA, IM model config file. 
                                        Migration from B to A. Third Ne for ancestral population (Ne_A_B) 
                    
        -h --help                   show this
        
"""

'''
- remove block size and get it from the tally
- remove kmax and get it from tally 
`
[--label] : defines where results of an analysis get stored!

[To do]
- write documentation (no automated filling in)

    Ne lower bounds = min(Pi_A / 4 mu, Pi_B / 4 mu) / 2
    Ne upper bounds = d_xy / (4 mu)

    T lower bound = 0
    T upper bound = dxy / (2 mu)

'''

import sys
from docopt import docopt
from timeit import default_timer as timer
import lib.gimble

MODELS_BY_INT = {
    '1': 'DIV',
    '2': 'MIG_AB',
    '3': 'MIG_BA',
    '4': 'IM_AB',
    '5': 'IM_BA'}
TASKS = set(['optimize', 'makegrid', 'simulate'])

class ModelParameterObj(lib.gimble.ParameterObj):
    '''Sanitises command line arguments and stores parameters'''

    def __init__(self, params, args):
        super().__init__(params)
        self.model = self._get_model(args)
        self.task = self._get_task(args)
        self.label = args['--label']

    def _get_task(self, args):
        if args['--task'] in TASKS:
            return args['--task']
        sys.exit('[X] Task must be one of the following: %s' % ", ".join(TASKS))

    def _get_model(self, args):
        if args['--model'] in MODELS_BY_INT:
            return MODELS_BY_INT[args['--model']]
        sys.exit('[X] Model must be one of the following: %s' % ", ".join(MODELS_BY_INT.keys()))

def main(params):
    try:
        start_time = timer()
        args = docopt(__doc__)
        parameterObj = ModelParameterObj(params, args)
        lib.gimble.write_config(
            version=parameterObj._VERSION,
            model=parameterObj.model,
            task=parameterObj.task,
            label=parameterObj.label
            )
        print("[*] Total runtime was %s" % (lib.gimble.format_time(timer() - start_time)))
    except KeyboardInterrupt:
        print("\n[X] Interrupted by user after %s !\n" % (lib.gimble.format_time(timer() - start_time)))
        exit(-1)