#!/usr/bin/env python2
# -*- coding: utf-8 -*-

"""usage: gIMble probs -o STR [-g FILE] [-p FILE] [-C FLOAT] [-M FLOAT] [-t INT] [-m FLOAT] [-h|--help]

    Options:
        -h --help                       show this
        -g, --graph FILE                State graph to analyse
        -p, --paths FILE                Paths to analyse
        -o, --out_prefix STR            Prefix for output files
        -M, --migration_rate FLOAT      Migration rate [default: 0.5]
        -C, --coalescence_rate FLOAT    Coalescence rate in [pop] relative to {pop} [default: 0.75]
        -t, --threads INT               Threads [default: 1]
        -m, --mutation_rate FLOAT       Mutation rate [default: 1.0]
"""

from __future__ import division
from docopt import docopt
from collections import defaultdict, Counter
from timeit import default_timer as timer
from sys import stderr, exit
from pandas import read_csv
from tqdm import tqdm
from sage.all import *
from multiprocessing import Pool
from contextlib import contextmanager
from itertools import combinations_with_replacement as combi_with_replacement 
from itertools import combinations as combi
from itertools import product as prod

'''
-m, --mutation_rate FLOAT       Mutation rate [default: 0.1]

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
- filter paths that can't accomodate mutations
- Marginals
    - keep mutation rates symbolic for path/mutation computations so that marginals can be calculated
    - iteratively turn off combinations of rates, None, A, B, AB, fixed, A+B, A+AB, A+fixed, ...A
    - save results in datastructure/pickle (4D)
- ILT
    - progressbar for SUM-ILT
    - implement switch for doing ILT-SUM vs SUM-ILT
    - implement interface for ILT-calculations for mpath (http://mpmath.org/) and pygiac (http://www-fourier.ujf-grenoble.fr/~parisse/giac_fr.html#python)
    - compare sage-giac/sage-sympy/sage-maxima/mpath/pygiac ILTs
'''

################################### Functions #################################

def multinomial(lst):
    res, i = 1, 1
    for a in lst:
        for j in range(1,a+1):
            res *= i
            res //= j
            i += 1
    return res

def weak_compositions(n, k):
    '''
    https://stackoverflow.com/questions/4647120/next-composition-of-n-into-k-parts-does-anyone-have-a-working-algorithm
    '''
    return ([t for t in combinations(range(n + k - 1), k - 1)])

def choose(n, k):
    """
    A fast way to calculate binomial coefficients by Andrew Dalke.
    """
    if 0 <= k <= n:
        ntok = 1
        ktok = 1
        for t in range(1, min(k, n - k) + 1):
            ntok *= n
            ktok *= t
            n -= 1
        return ntok // ktok
    else:
        return 0


@contextmanager
def poolcontext(*args, **kwargs):
    pool = Pool(*args, **kwargs)
    yield pool
    pool.terminate()

def calculate_equation(params):
    data_idx, data_point, pathObj, max_by_mutype, rates = params
    mutation_equation = pathObj.infer_mutation_equation(data_point, rates)
    equation = pathObj.multiply_with_path_equation(mutation_equation)
    equation_ilt = inverse_laplace(equation / rates['E'], var('E'), var('T'), algorithm='giac').substitute(T=1.4)
    return (data_idx, pathObj.path_id, equation_ilt)
###############################################################################

################################### Classes ###################################

class NodeObj(object):
    def __init__(self, node_id, mutype_counter, event_counter):
        self.node_id = node_id
        self.mutype_counter = mutype_counter
        self.event_counter = event_counter

class PathObj(object):
    def __init__(self, path_id):
        self.path_id = path_id
        self.nodeObjs = []
        self.event_ids = []
        self.event_counts = []
        self.path_probability = 1
        self.mutable_nodes_by_mutype = defaultdict(list)
        self.path_numerators = []
        self.path_denominators = []

    def add_step(self, nodeObj, event_id, event_count):
        self.nodeObjs.append(nodeObj)
        self.event_ids.append(event_id)
        self.event_counts.append(event_count)

    def infer_path_probability(self, rates):
        for nodeObj, event_id, event_count in self.yield_step():
            self.path_numerators.append((event_count * rates[event_id]))
            self.path_denominators.append(
                        (nodeObj.event_counter['C_ancestor'] * rates['C_ancestor']) + 
                        (nodeObj.event_counter['C_derived'] * rates['C_derived']) + 
                        (nodeObj.event_counter['M'] * rates['M']) + 
                        (nodeObj.event_counter['E'] * rates['E']) +
                        (nodeObj.mutype_counter['hetA'] * rates['theta']['hetA']) +
                        (nodeObj.mutype_counter['hetB'] * rates['theta']['hetB']) +
                        (nodeObj.mutype_counter['hetAB'] * rates['theta']['hetAB']) +
                        (nodeObj.mutype_counter['fixed'] * rates['theta']['fixed'])
                )

    def multiply_with_path_equation(self, mutation_equation):
        initial = 1
        #print("[+] numerators:", self.path_numerators)
        #print("[+] denominators:", self.path_denominators)
        for numerator, denominator in zip(self.path_numerators, self.path_denominators):
            initial *= (numerator / denominator)
        mutation_equation *= initial
        return mutation_equation

    def infer_mutation_equation(self, data_point, rates):
        #print("#data:", data_point, "at path", self.path_id)
        MUTYPES = ['hetA', 'fixed', 'hetB', 'hetAB']
        # theta = {'hetA': var('theta_A'), 'fixed': var('theta_fixed'), 'hetB': var('theta_B'), 'hetAB': var('theta_AB')}
        mutation_equation = 0
        for slots in prod(*[list(combi_with_replacement(self.mutable_nodes_by_mutype[mutype], data_point[mutype])) for mutype in MUTYPES]):
            # each slots is a way of placing data_point on given path
            #print("[+] Slots", slots)
            mutypes_by_idx = defaultdict(list)
            # mutypes_by_idx is number of ways a mutype can be placed at node
            for mutype, slot in zip(MUTYPES, slots):
                if slot:
                    for idx in list(slot):
                        mutypes_by_idx[idx].append(mutype)
            mutation_part = 1
            for node_idx, mutypes in mutypes_by_idx.items():
                #print("[+] Node %s (%s)" % (node_idx, mutypes))
                mutype_counter = Counter(mutypes)
                #print("[+] Counts: %s" % mutype_counter)
                node_multinomial = multinomial(mutype_counter.values())
                mutation_part *= node_multinomial
                #print("[+] Multinomial equation for node %s (based on %s) = %s" % (node_idx, mutype_counter.values(), node_multinomial))
                nodeObj = self.nodeObjs[node_idx]
                denominator = self.path_denominators[node_idx]
                #print("[+] Path equation for node %s" % (node_idx), self.path_denominators[node_idx])
                #for mutype in mutypes_by_idx[node_idx]:
                #    denominator += (mutype_counter[mutype] * rates['theta'][mutype])
                #    print("[+] Denominator of equation for node %s" % (node_idx), denominator)
                for mutype in mutypes:
                    mutation_part *= ((nodeObj.mutype_counter[mutype] * rates['theta'][mutype]) / denominator)
                    #print("[+] Mutation equation for %s" % (mutype))
            #print("[+] mutation equation of slot: ", mutation_part)
            mutation_equation += mutation_part
        return mutation_equation

    def infer_mutable_nodes(self):
        for idx, nodeObj in enumerate(self.nodeObjs):
            for mutype, count in nodeObj.mutype_counter.items():
                if count:
                    self.mutable_nodes_by_mutype[mutype].append(idx)

    def yield_step(self):
        for nodeObj, event_id, event_count in zip(self.nodeObjs, self.event_ids, self.event_counts):
            yield (nodeObj, event_id, event_count)

    def __str__(self):
        return "# Path %s:\n%s" % (self.path_id, 
            "\n".join(["%s: %s %s" % (nodeObj.node_id, event_id, str(nodeObj.mutype_counter)) for nodeObj, event_id in zip(self.nodeObjs, self.event_ids)])
            )

class ParameterObj(object):
    def __init__(self, args):
        self.graph_file = args['--graph']
        self.path_file = args['--paths']
        self.out_prefix = args['--out_prefix']
        self.threads = int(args['--threads'])
        self.max_by_mutype = {'hetA': 2, 'hetB': 2, 'hetAB': 2, 'fixed': 2}
        self.rates = {
            'C_ancestor': 1, 
            'C_derived': 1.3, 
            'M': 2.34, 
            'E': var('E'),
            'theta' : {
                'hetA': 0.6,
                'hetB': 0.6,
                'hetAB': 0.6,
                'fixed': 0.6
                }
            }
        self.symbolic_rates = {
            'C_ancestor': var('C_anc'), 
            'C_derived': var('C_der'), 
            'M': var('M'),
            'E': var('E'),
            'theta' : {
                'hetA': var('theta_A'),
                'hetB': var('theta_B'),
                'hetAB': var('theta_BA'),
                'fixed': var('theta_fix'),
                }
            }

    def read_paths(self):
        start_time = timer()
        print("[+] Reading paths in file %s..." % self.path_file)
        paths_df = read_csv(\
            self.path_file, \
            sep="\t", \
            )
        pathObj_by_path_id = {}
        nodeObj_by_node_id = {}
        for path_id, node_id, event_id, event_count, C_ancestor, C_derived, Ms, Es, hetA, fixed, hetB, hetAB in tqdm(paths_df.values.tolist(), total=len(paths_df.index), desc="[%] ", ncols=200):
            if not event_id == 'LCA':
                if not node_id in nodeObj_by_node_id:
                    mutype_counter = Counter({'hetA': hetA, 'fixed': fixed, 'hetB': hetB, 'hetAB': hetAB})
                    event_counter = Counter({'C_ancestor': C_ancestor, 'C_derived': C_derived, 'M': Ms, 'E': Es})
                    nodeObj_by_node_id[node_id] = NodeObj(node_id, mutype_counter, event_counter)
                # paths
                if not path_id in pathObj_by_path_id:
                    pathObj_by_path_id[path_id] = PathObj(path_id)
                pathObj_by_path_id[path_id].add_step(nodeObj_by_node_id[node_id], event_id, event_count)
        print("[+] Parsed paths in file %s in %s seconds." % (self.path_file, timer() - start_time))
        return (pathObj_by_path_id, nodeObj_by_node_id) 

    def analyse_paths(self, data, pathObj_by_path_id, nodeObj_by_node_id):
        start_time = timer()
        equation_ilt_by_path_id_by_data_idx = defaultdict(dict)
        #print("[+] Analysing %s paths and %s data points (%s * %s = %s) with %s threads..." % (len(pathObj_by_idx), len(data), len(pathObj_by_idx), len(data), len(pathObj_by_idx) * len(data), self.threads))
        params = []

        rates = self.rates
        #rates = self.symbolic_rates

        # infer path_probabilities for each pathObj (only has to be done once) ... 
        for path_id, pathObj in pathObj_by_path_id.items():
            pathObj.infer_path_probability(rates)
            pathObj.infer_mutable_nodes()
            #print(pathObj)
        
        # deal with mutations
        for data_idx, data_point in enumerate(data):
            #print(data_idx, data_point)
            for path_id, pathObj in pathObj_by_path_id.items():
                # can mutations be placed on path
                #mutypes_in_data = [mutype for mutype, count in data_point.items() if count > 0]
                # Mutations can be placed on given path
                params.append((data_idx, data_point, pathObj, self.max_by_mutype, rates))
                #    else:
                #        print("[-] Path %s (%s) can't accomodate mutations: %s" % (pathObj.path_id, pathObj.mutable_nodes_by_mutype, data_point))
        #print(len(params))
        if self.threads == 1:
            for param in tqdm(params, total=len(params), desc="[%] ", ncols=200):
                data_idx, path_id, equation_ilt = calculate_equation(param)
                equation_ilt_by_path_id_by_data_idx[data_idx][path_id] = equation_ilt
        else:
            with poolcontext(processes=self.threads) as pool:
                with tqdm(params, total=len(params), desc="[%] ", ncols=200) as pbar:
                    for data_idx, path_id, equation_ilt in pool.imap_unordered(calculate_equation, params):
                        equation_ilt_by_path_id_by_data_idx[data_idx][path_id] = equation_ilt
                        pbar.update()
        print("[+] Analysed paths in %s seconds." % (timer() - start_time))
        return equation_ilt_by_path_id_by_data_idx

###############################################################################

################################### Main ######################################

def main():
    try:
        main_time = timer()
        args = docopt(__doc__)
        print(args)
        parameterObj = ParameterObj(args)
        pathObj_by_path_id, nodeObj_by_node_id = parameterObj.read_paths()
        data = [
            Counter({'hetA':2, 'hetB':1, 'fixed':0, 'hetAB':0}),
            Counter({'hetA':0, 'hetB':0, 'fixed':0, 'hetAB':0}),
            Counter({'hetA':0, 'hetB':1, 'fixed':0, 'hetAB':0}),
            Counter({'hetA':0, 'hetB':2, 'fixed':0, 'hetAB':0}),
            Counter({'hetA':0, 'hetB':0, 'fixed':0, 'hetAB':1}),
            Counter({'hetA':0, 'hetB':1, 'fixed':0, 'hetAB':1}),
            Counter({'hetA':0, 'hetB':2, 'fixed':0, 'hetAB':1}),
            Counter({'hetA':0, 'hetB':0, 'fixed':0, 'hetAB':2}),
            Counter({'hetA':0, 'hetB':1, 'fixed':0, 'hetAB':2}),
            Counter({'hetA':0, 'hetB':2, 'fixed':0, 'hetAB':2}),
            Counter({'hetA':0, 'hetB':0, 'fixed':1, 'hetAB':0}),
            Counter({'hetA':0, 'hetB':1, 'fixed':1, 'hetAB':0}),
            Counter({'hetA':0, 'hetB':2, 'fixed':1, 'hetAB':0}),
            Counter({'hetA':0, 'hetB':0, 'fixed':2, 'hetAB':0}),
            Counter({'hetA':0, 'hetB':1, 'fixed':2, 'hetAB':0}),
            Counter({'hetA':0, 'hetB':2, 'fixed':2, 'hetAB':0}),
            Counter({'hetA':1, 'hetB':1, 'fixed':0, 'hetAB':0}),
            Counter({'hetA':1, 'hetB':2, 'fixed':0, 'hetAB':0}),
            Counter({'hetA':1, 'hetB':0, 'fixed':0, 'hetAB':1}),
            Counter({'hetA':1, 'hetB':1, 'fixed':0, 'hetAB':1}),
            Counter({'hetA':1, 'hetB':2, 'fixed':0, 'hetAB':1}),
            Counter({'hetA':1, 'hetB':0, 'fixed':0, 'hetAB':2}),
            Counter({'hetA':1, 'hetB':1, 'fixed':0, 'hetAB':2}),
            Counter({'hetA':1, 'hetB':2, 'fixed':0, 'hetAB':2}),
            Counter({'hetA':1, 'hetB':0, 'fixed':1, 'hetAB':0}),
            Counter({'hetA':1, 'hetB':1, 'fixed':1, 'hetAB':0}),
            Counter({'hetA':1, 'hetB':2, 'fixed':1, 'hetAB':0}),
            Counter({'hetA':1, 'hetB':0, 'fixed':2, 'hetAB':0}),
            Counter({'hetA':1, 'hetB':1, 'fixed':2, 'hetAB':0}),
            Counter({'hetA':1, 'hetB':2, 'fixed':2, 'hetAB':0}),
            Counter({'hetA':2, 'hetB':0, 'fixed':0, 'hetAB':0}),
            Counter({'hetA':2, 'hetB':2, 'fixed':0, 'hetAB':0}),
            Counter({'hetA':2, 'hetB':0, 'fixed':0, 'hetAB':1}),
            Counter({'hetA':2, 'hetB':1, 'fixed':0, 'hetAB':1}),
            Counter({'hetA':2, 'hetB':2, 'fixed':0, 'hetAB':1}),
            Counter({'hetA':2, 'hetB':0, 'fixed':0, 'hetAB':2}),
            Counter({'hetA':2, 'hetB':1, 'fixed':0, 'hetAB':2}),
            Counter({'hetA':2, 'hetB':2, 'fixed':0, 'hetAB':2}),
            Counter({'hetA':2, 'hetB':0, 'fixed':1, 'hetAB':0}),
            Counter({'hetA':2, 'hetB':1, 'fixed':1, 'hetAB':0}),
            Counter({'hetA':2, 'hetB':2, 'fixed':1, 'hetAB':0}),
            Counter({'hetA':2, 'hetB':0, 'fixed':2, 'hetAB':0}),
            Counter({'hetA':2, 'hetB':1, 'fixed':2, 'hetAB':0}),
            Counter({'hetA':2, 'hetB':2, 'fixed':2, 'hetAB':0})
        ]

        equation_ilt_by_path_id_by_data_idx = parameterObj.analyse_paths(data, pathObj_by_path_id, nodeObj_by_node_id)
        ilvs = []
        for data_idx in equation_ilt_by_path_id_by_data_idx.keys():
            ilv = sum(equation_ilt_by_path_id_by_data_idx[data_idx].values())
            # sum_equation = sum(equation_ilt_by_path_id_by_data_idx[data_idx].values())
            # ilv = inverse_laplace(sum_equation / var('E'), var('E'), var('T'), algorithm='giac').substitute(T=1.4)
            print(data[data_idx], "P(all)=%s" % ilv)
            ilvs.append(ilv)
        print("Total:", sum(ilvs))
    except KeyboardInterrupt:
        stderr.write("\n[X] Interrupted by user after %s seconds!\n" % (timer() - main_time))
        exit(-1)
        
###############################################################################

if __name__ == '__main__':
    main()