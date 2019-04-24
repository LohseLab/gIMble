#!/usr/bin/env python2
# -*- coding: utf-8 -*-

"""usage: gIMble probs         -p FILE [-A FLOAT] [-D FLOAT] [-M FLOAT] [-m FLOAT] [-T FLOAT] [-t INT] [-h|--help]

    Options:
        -h --help                                   show this
        -p, --paths FILE                            Paths to analyse
        -A, --ancestor_coalescence_rate FLOAT   Coalescence rate in [ancestor] pop [default: 1.0]
        -D, --derived_coalescence_rate FLOAT    Coalescence rate in {derived} pop [default: 1.0]
        -M, --migration_rate FLOAT                  Migration rate per generation [default: 2.35]
        -m, --mutation_rate FLOAT                   Mutation rate/lineage [default: 0.6]
        -T, --T_var FLOAT                           T [default: 1.4]
        -t, --threads INT                           Threads [default: 1]
        
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
from itertools import product as prod
from pandas import DataFrame

'''
[To do]
- filter paths that can't accomodate mutations
- Marginals
    - keep mutation rates symbolic for path/mutation computations so that marginals can be calculated
    - iteratively turn off combinations of rates, None, A, B, AB, fixed, A+B, A+AB, A+fixed, ...A
    - save results in datastructure/pickle (4D)
- ILT
    - implement interface for ILT-calculations for mpath (http://mpmath.org/) and pygiac (http://www-fourier.ujf-grenoble.fr/~parisse/giac_fr.html#python)
    - compare sage-giac/sage-sympy/sage-maxima/mpath/pygiac ILTs
'''

################################### CONSTANTS #################################

MUTYPES = ['hetA', 'fixed', 'hetB', 'hetAB']

###############################################################################

################################### Functions #################################

def create_csv(out_f, header, rows, sep):
    df = DataFrame(rows, columns=header)
    print(df)
    df.to_csv(out_f, index=False, sep=sep)

def multinomial(lst):
    '''
    function to calculate multinomial coefficient from list, e.g.: [2,1]
    https://stackoverflow.com/questions/46374185/does-python-have-a-function-which-computes-multinomial-coefficients
    '''
    res, i = 1, 1
    for a in lst:
        for j in range(1,a+1):
            res *= i
            res //= j
            i += 1
    return res

@contextmanager
def poolcontext(*args, **kwargs):
    pool = Pool(*args, **kwargs)
    yield pool
    pool.terminate()

def calculate_ilt(params):
    data_idx, data_point, pathObj, max_by_mutype, rates = params
    mutation_equation = pathObj.infer_mutation_equation(data_point, rates)
    equation = pathObj.multiply_with_path_equation(mutation_equation)
    equation_ilt = inverse_laplace(equation / rates['E'], var('E'), var('T'), algorithm='giac').substitute(T=rates['T'])
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
        for numerator, denominator in zip(self.path_numerators, self.path_denominators):
            initial *= (numerator / denominator)
        mutation_equation *= initial
        return mutation_equation

    def infer_mutation_equation(self, data_point, rates):
        mutation_equation = 0
        for slots in prod(*[list(combi_with_replacement(self.mutable_nodes_by_mutype[mutype], data_point[mutype])) for mutype in MUTYPES]):
            # each slots is a way of placing data_point on given path
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
        self.path_file = args['--paths']
        self.out_prefix = ".".join(args['--paths'].split(".")[0:-1])
        self.threads = int(args['--threads'])
        self.max_by_mutype = {'hetA': 2, 'hetB': 2, 'hetAB': 2, 'fixed': 2}
        self.rates = {
            'M': float(args['--migration_rate']),
            'C_ancestor': float(args['--ancestor_coalescence_rate']),
            'C_derived': float(args['--derived_coalescence_rate']),
            'T': float(args['--T_var']),
            'E': var('E'),
            'm': float(args['--mutation_rate']),
            'theta' : {
                'hetA': float(args['--mutation_rate']),
                'hetB': float(args['--mutation_rate']),
                'hetAB': float(args['--mutation_rate']),
                'fixed': float(args['--mutation_rate'])
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

    def write_probs(self, prob_by_data_idx, data):
        header = MUTYPES + ['probability']
        rows = []
        prob_sum = 0
        for data_idx, prob in prob_by_data_idx.items():
            rows.append([data[data_idx][mutype] for mutype in MUTYPES] + [prob])
            prob_sum += prob
        rows.append(['*', '*', '*', '*',  prob_sum])
        out_f = "%s.C_anc=%s.C_der=%s.M=%s.theta=%s.T=%s.t=%s.tsv" % (self.out_prefix, self.rates['C_ancestor'], self.rates['C_derived'], self.rates['M'], self.rates['m'], self.rates['T'], self.threads)
        create_csv(out_f, header, rows, "\t")
            

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
        
        print("[+] Analysing %s paths and %s data points (%s * %s = %s) with %s threads..." % (len(pathObj_by_path_id), len(data), len(pathObj_by_path_id), len(data), len(pathObj_by_path_id) * len(data), self.threads))
        rates = self.rates
        #rates = self.symbolic_rates

        # infer path_probabilities for each pathObj (only has to be done once) ... 
        for path_id, pathObj in pathObj_by_path_id.items():
            pathObj.infer_path_probability(rates)
            pathObj.infer_mutable_nodes()
        
        # deal with mutations
        params = []
        for data_idx, data_point in enumerate(data):
            for path_id, pathObj in pathObj_by_path_id.items():
                # can mutations be placed on path
                #mutypes_in_data = [mutype for mutype, count in data_point.items() if count > 0]
                # Mutations can be placed on given path
                params.append((data_idx, data_point, pathObj, self.max_by_mutype, rates))
                #    else:
                #        print("[-] Path %s (%s) can't accomodate mutations: %s" % (pathObj.path_id, pathObj.mutable_nodes_by_mutype, data_point))
        #print(len(params))
        ilt_by_path_id_by_data_idx = defaultdict(dict)
        if self.threads == 1:
            for param in tqdm(params, total=len(params), desc="[%] ", ncols=200):
                data_idx, path_id, ilt = calculate_ilt(param)
                ilt_by_path_id_by_data_idx[data_idx][path_id] = ilt
        else:
            with poolcontext(processes=self.threads) as pool:
                with tqdm(params, total=len(params), desc="[%] ", ncols=200) as pbar:
                    for data_idx, path_id, ilt in pool.imap_unordered(calculate_ilt, params):
                        ilt_by_path_id_by_data_idx[data_idx][path_id] = ilt
                        pbar.update()
        prob_by_data_idx = {}
        for data_idx in ilt_by_path_id_by_data_idx.keys():
            prob = sum(ilt_by_path_id_by_data_idx[data_idx].values())
            prob_by_data_idx[data_idx] = prob
        print("[+] Analysed paths in %s seconds." % (timer() - start_time))
        return prob_by_data_idx

###############################################################################

################################### Main ######################################

def main():
    try:
        main_time = timer()
        args = docopt(__doc__)
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
        prob_by_data_idx = parameterObj.analyse_paths(data, pathObj_by_path_id, nodeObj_by_node_id)
        parameterObj.write_probs(prob_by_data_idx, data)
        print("[+] Total runtime: %s seconds" % (timer() - main_time))
    except KeyboardInterrupt:
        stderr.write("\n[X] Interrupted by user after %s seconds!\n" % (timer() - main_time))
        exit(-1)
        
###############################################################################

if __name__ == '__main__':
    main()