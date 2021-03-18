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
from sys import stderr, exit
import gimblelib.model
import lib.gimble

MODEL_BY_ARG = {
    '--DIV': 'DIV',
    '--MIG_AB': 'MIG_AB',
    '--MIG_BA': 'MIG_BA',
    '--IM_AB': 'IM_AB',
    '--IM_BA': 'IM_BA'
        }
JOIN_EVENTS_BY_MODEL = {
    'DIV': ['A_B'],
    'MIG_AB': [],
    'MIG_BA': [],
    'IM_AB': ['A_B'],
    'IM_BA': ['A_B']
}
MIGRATION_EVENTS_BY_MODEL = {
    'DIV': [],
    'MIG_AB': ['M_A_B'],
    'MIG_BA': ['M_B_A'],
    'IM_AB': ['M_A_B'],
    'IM_BA': ['M_B_A']
}

class ModelParameterObj(lib.gimble.ParameterObj):
    '''Sanitises command line arguments and stores parameters'''

    def __init__(self, params, args):
        super().__init__(params)
        self.ploidy = 2
        self.pop_ids = ['A', 'B']
        self.model = [model for arg, model in MODEL_BY_ARG.items() if args[arg]][0]
        self.samples_by_pop_id = {pop_id: 1 for pop_id in self.pop_ids}
        self.join_events = JOIN_EVENTS_BY_MODEL[self.model]
        self.migration_events = MIGRATION_EVENTS_BY_MODEL[self.model]
        self.sample_stateObj = self._get_sample_stateObj()
        self.set_of_events = self._get_set_of_events() 
        self.outfile = self._get_outfile(args['--outprefix'])

    def _get_outfile(self, outprefix):
        return ".".join([outprefix, self.model, 'ini'])
        
    def _get_ploidy(ploidy):
        return int(ploidy)
        
    def _get_set_of_events(self):
        return set(["C_%s" % pop_id for pop_id in self.pop_ids + self.join_events] + self.join_events + self.migration_events)

    def _get_sample_stateObj(self):
        pop_dict = {}
        for pop_id in self.pop_ids:
            pop_list = []
            for i in range(self.samples_by_pop_id[pop_id]):
                for j in range(self.ploidy):
                    pop_list.append("%s%s" % (pop_id.lower(), i+1))
            pop_dict[pop_id.upper()] = gimblelib.model.PopObj(pop_list)
        return gimblelib.model.StateObj(pop_dict)

def main(params):
    try:
        start_time = timer()
        args = docopt(__doc__)
        parameterObj = ModelParameterObj(params, args)
        stateGraph = gimblelib.model.graph_generator(parameterObj)
        stateGraph.write_ini(parameterObj)
        print("[*] Total runtime was %s" % (lib.gimble.format_time(timer() - start_time)))
    except KeyboardInterrupt:
        print("\n[X] Interrupted by user after %s !\n" % (lib.gimble.format_time(timer() - start_time)))
        exit(-1)