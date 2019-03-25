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
# from collections import Counter, defaultdict
import matplotlib as mat
mat.use("agg")
from timeit import default_timer as timer
from sys import stderr, exit
from ast import literal_eval
from pandas import DataFrame, read_csv
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
- remove state graph?

- Datastructure to know which edges have a given sub-state => allows knowing paths on which mutations can be placed
'''

################################### Functions #################################

def create_csv(out_f, header, rows, sep):
    df = DataFrame(rows, columns=header)
    print(df)
    df.to_csv(out_f, index=False, sep=sep)

def pairwise(iterable):
    it = iter(iterable)
    a = next(it, None)
    for b in it:
        yield (a, b)
        a = b

###############################################################################

################################### Classes ###################################

class ParameterObj(object):
    def __init__(self, args):
        self.graph_file = args['--graph']
        self.path_file = args['--paths']
        self.out_prefix = args['--out_prefix']
        self.coalescence_rate_b = float(args['--coalescence_rate'])
        self.migration_rate = float(args['--migration_rate'])
        
    def write_paths(self, header, paths):
        start_time = timer()
        out_f = "%s.paths.txt" % (self.get_basename())
        create_csv(out_f, header, paths, sep="\t")
        print("[+] Written paths into file %s in %s seconds." % (out_f, timer() - start_time))

    def read_paths(self):
        start_time = timer()
        paths_df = read_csv(\
            self.path_file, \
            sep="\t", \
            )
        seen_paths = set() 
        # probably doable with defaultdict(1)
        denominator = 1
        numerator = 1
        probs = []
        T_var = var('T', domain='positive')
        R_var = var('R')
        E_var = var('E')
        for path_idx, path, step_idx, state, event, event_count, C_mig, C_res, Ms, Es, Rs in tqdm(paths_df.values.tolist(), total=len(paths_df.index), desc="[%] ", ncols=200):
            # setup of symbolic variables
            if event == 'LCA':
                continue
            if path_idx not in seen_paths:
                if len(seen_paths) > 0:
                    print((numerator/denominator)/R_var)
                    print("P(%s)=%s" % (path_idx, inverse_laplace((numerator/denominator)/R_var, R_var, T_var, algorithm='giac'))) # make it figure out which var (E/R) to declare/use
                    probs.append(inverse_laplace((numerator/denominator)/R_var, R_var, T_var, algorithm='giac'))
                seen_paths.add(path_idx)
                denominator = 1
                numerator = 1
            #print(path_idx, path, step_idx, state, event, event_count, C_mig, C_res, Ms, Es, Rs)
            counter = {'{C}': C_mig, '[C]': C_res, 'M': Ms, 'E': Es, 'R': Rs}
            if event == '{C}':
                numerator *= event_count 
            elif event == '[C]':
                numerator *= event_count * self.coalescence_rate_b
            elif event == 'M':
                numerator *= event_count * self.migration_rate
            elif event == 'E':
                numerator *= event_count * E_var
            elif event == 'R':
                numerator *= event_count * R_var
            else:
                exit("[X] Should never happen ...")
            denominator *= (counter['{C}']) + (counter['[C]'] * self.coalescence_rate_b) + (counter['M'] * self.migration_rate) + counter['E'] * E_var + counter['R'] * R_var
        print((numerator/denominator)/R_var)
        print("P(%s)=%s" % (path_idx, inverse_laplace((numerator/denominator)/R_var, R_var, T_var, algorithm='giac'))) # make it figure out which var (E/R) to declare/use
        probs.append(inverse_laplace((numerator/denominator)/R_var, R_var, T_var, algorithm='giac'))
        print("[+] Parsed paths in file %s in %s seconds." % (self.path_file, timer() - start_time))
        print(sum(probs).substitute(T=3.14))
    
    def write_graph(self, state_graph):
        start_time = timer()
        out_f = "%s.gml" % (self.get_basename())
        # states have to be stringified ...
        for idx, node in state_graph.nodes(data=True):
            # print(idx, node)
            node['state'] = str(node.get('state', None))
        nx.write_gml(state_graph, out_f)
        print("[+] Written graph into file %s in %s seconds." % (out_f, timer() - start_time))
        return out_f

    def read_graph(self):
        start_time = timer()
        state_graph = nx.read_gml(self.graph_file)
        for idx, node in state_graph.nodes(data=True):
            node['state'] = literal_eval(node['state'])
            if 'pos' in node:
                node['pos'] = tuple(node['pos'])
        print("[+] Read graph from file %s in %s seconds." % (self.graph_file, timer() - start_time))
        return state_graph

###############################################################################

################################### Main ######################################

def main():
    try:
        main_time = timer()
        args = docopt(__doc__)
        print(args)
        parameterObj = ParameterObj(args)
        #state_graph = parameterObj.read_graph()
        path = parameterObj.read_paths()

    except KeyboardInterrupt:
        stderr.write("\n[X] Interrupted by user after %s seconds!\n" % (timer() - main_time))
        exit(-1)
        
###############################################################################

if __name__ == '__main__':
    main()