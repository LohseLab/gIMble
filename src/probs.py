#!/usr/bin/env python2
# -*- coding: utf-8 -*-

"""usage: gIMble probs -o STR [-g FILE] [-p FILE] [-C FLOAT] [-M FLOAT] [-E FLOAT] [-h|--help]

    Options:
        -h --help                       show this
        -g, --graph FILE                State graph to analyse
        -p, --paths FILE                Paths to analyse
        -o, --out_prefix STR            Prefix for output files
        -M, --migration_rate FLOAT      Migration rate [default: 0.5]
        -C, --coalescence_rate FLOAT    Coalescence rate in [pop] relative to {pop} [default: 0.75]
        
"""

from docopt import docopt
from collections import defaultdict
from timeit import default_timer as timer
from sys import stderr, exit
from pandas import read_csv
from tqdm import tqdm
from sage.all import *


'''
[Problems]
- Ipython notebook does not work on my machine, make general implementation of task
- sage needs python2 !?!
    - hacky compiling with python3 : https://wiki.sagemath.org/Python3-compatible%20code
    - official version soon: 
        - https://trac.sagemath.org/ticket/26212
        - https://trac.sagemath.org/ticket/15530
    sage 8.6 (?)
        - Add the conda-forge channel: conda config --add channels conda-forge
        - Update all packages: conda update --all
        - Create a new environment containing SageMath: conda create -n sage sage
        - Enter the new environment: source activate sage
        - Start SageMath: sage

[To do]
- Input path table
    - Is there a need for keeping E and R separate? could be one event? then we would have one 'rate' and no need for if/else's ...

- remove state graph?
- Datastructure to know which edges have a given sub-state => allows knowing paths on which mutations can be placed
'''

################################### Functions #################################

# TBD

###############################################################################

################################### Classes ###################################

class ParameterObj(object):
    def __init__(self, args):
        self.graph_file = args['--graph']
        self.path_file = args['--paths']
        self.out_prefix = args['--out_prefix']
        self.coalescence_rate_b = float(args['--coalescence_rate'])
        self.migration_rate = float(args['--migration_rate'])

    def read_paths(self):
        start_time = timer()
        paths_df = read_csv(\
            self.path_file, \
            sep="\t", \
            )
        # probably doable with defaultdict(1)
        probs = []
        T_var = var('T')
        equation_raw_by_path_idx = defaultdict(lambda: 1)
        equation_ilt_by_path_idx = {}
        R_var = var('R')
        # E_var = var('E')
        rates = {'{C}': 1, '[C]': self.coalescence_rate_b, 'M': self.migration_rate, 'E': var('E'), 'R': var('R')}
        print(rates)
        for path_idx, path, step_idx, state, event, event_count, C_mig, C_res, Ms, Es, Rs in tqdm(paths_df.values.tolist(), total=len(paths_df.index), desc="[%] ", ncols=200):
            if event == 'LCA':
                # LCA is end of path
                equation_ilt_by_path_idx[path_idx] = inverse_laplace(equation_raw_by_path_idx[path_idx] / R_var, R_var, T_var, algorithm='giac')
                # print("P(%s)=%s" % (path_idx, equation_ilt_by_path_idx[path_idx])) 
                probs.append(equation_ilt_by_path_idx[path_idx])
            else:
                counter = {'{C}': C_mig, '[C]': C_res, 'M': Ms, 'E': Es, 'R': Rs}
                equation_raw_by_path_idx[path_idx] *= (event_count * rates[event]) / ((counter['{C}']) + (counter['[C]'] * rates['[C]']) + (counter['M'] * rates['M']) + counter['E'] * rates['E'] + counter['R'] * rates['R'])
        print("[+] Parsed paths in file %s in %s seconds." % (self.path_file, timer() - start_time))
        return equation_ilt_by_path_idx

###############################################################################

################################### Main ######################################

def main():
    try:
        main_time = timer()
        args = docopt(__doc__)
        print(args)
        parameterObj = ParameterObj(args)
        equation_ilt_by_path_idx = parameterObj.read_paths()
        print(sum(equation_ilt_by_path_idx.values()).substitute(T=3.14))

    except KeyboardInterrupt:
        stderr.write("\n[X] Interrupted by user after %s seconds!\n" % (timer() - main_time))
        exit(-1)
        
###############################################################################

if __name__ == '__main__':
    main()