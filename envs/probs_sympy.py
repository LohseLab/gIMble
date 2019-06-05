#!/usr/bin/env python2
# -*- coding: utf-8 -*-

"""usage: gIMble probs_sympy                      -p FILE [-h|--help]
                                            [-A FLOAT] [-D FLOAT] [-M FLOAT] [-m FLOAT] [-T FLOAT] 
                                            [-t INT] [-P INT] 

    Options:
        -h --help                                   show this
        -p, --paths FILE                            Paths to analyse
        -A, --ancestor_coalescence_rate FLOAT       Coalescence rate in [ancestor] pop [default: 1.0]
        -D, --derived_coalescence_rate FLOAT        Coalescence rate in {derived} pop [default: 1.0]
        -M, --migration_rate FLOAT                  Migration rate per generation [default: 2.35]
        -m, --mutation_rate FLOAT                   Mutation rate/lineage [default: 0.6]
        -T, --T_var FLOAT                           T [default: 1.4]
        -P, --precision INT                         Floating point precision of probabilities [default: 30]
        -t, --threads INT                           Threads [default: 1]
        
"""

from __future__ import division
from docopt import docopt
from collections import defaultdict, Counter, OrderedDict
from timeit import default_timer as timer
from sys import stderr, exit
from pandas import read_csv
from tqdm import tqdm
from sage.all import *
from scipy.special import xlogy
from scipy.optimize import minimize
import numpy as np
from multiprocessing import Pool
from contextlib import contextmanager
from itertools import combinations_with_replacement as combi_with_replacement
from itertools import combinations as combi
from itertools import chain as iter_chain 
from itertools import product as prod
from pandas import DataFrame

'''
[To do]
- filter paths that can't accomodate mutations √
- Marginals
    - keep mutation rates symbolic for path/mutation computations so that marginals can be calculated √
    - save results in datastructure/pickle (4D)
- ILT
    - implement interface for ILT-calculations for mpath (http://mpmath.org/) and pygiac (http://www-fourier.ujf-grenoble.fr/~parisse/giac_fr.html#python)
    - compare sage-giac/sage-sympy/sage-maxima/mpath/pygiac ILTs
'''

################################### CONSTANTS #################################

MUTYPES = ['hetA', 'fixed', 'hetB', 'hetAB']

###############################################################################

################################### Functions #################################

def calculate_likelihood(x0, args):
        start_time = timer()
        print("[+] Calculating probabilities (datapoints and marginals)...")
        # unpack parameters
        data, symbolic_equations_by_data_idx, parameterObj.base_rate, params = args
        base_rate, mutation_rate, split_time = params
        split_time = x0[0]
        base_rate['M'] = x0[1]

        # symbolic_equations_by_data_idx are completely symbolic 
        equation_by_data_idx = {data_idx: sum(equation).subs(base_rate) for data_idx, equation in symbolic_equations_by_data_idx.items()}


        # Setting up datastructure
        shape = tuple(self.max_by_mutype[mutype] + 2 for mutype in MUTYPES)
        A = np.zeros(shape, np.float64)
        prob_by_data_string = OrderedDict()
        data_queries_seen = set()
        theta_var_by_mutype = {
            'hetA' : var('theta_A'),
            'fixed' : var('theta_fix'),
            'hetB' : var('theta_B'),
            'hetAB' : var('theta_AB')
        }
        # Generating sets of mutation_rate combinations for marginals
        sets_of_mutation_rates = list(powerset(MUTYPES))
        for mutypes in tqdm(reversed(sets_of_mutation_rates), total=len(sets_of_mutation_rates), desc="[%] ", ncols=200):
            mutypes_masked = [mutype for mutype in MUTYPES if not mutype in mutypes]
            for data_idx, equation in equation_by_data_idx.items():                
                data_point = data[data_idx]
                rate_by_mutype = {}
                _data_vector = []
                _data_string = []
                for mutype in MUTYPES:
                    if mutype in mutypes:
                        rate_by_mutype[theta_var_by_mutype[mutype]] = mutation_rate
                        _data_vector.append(data_point[mutype]) 
                        _data_string.append(data_point[mutype])
                    else:
                        rate_by_mutype[theta_var_by_mutype[mutype]] = 0
                        _data_vector.append(self.max_by_mutype[mutype] + 1)
                        _data_string.append("%s+" % (self.max_by_mutype[mutype] + 1))
                    data_vector = tuple(_data_vector)
                    data_string = tuple(_data_string)
                if not mutypes_masked: 
                    if not data_string in data_queries_seen:
                        probability_ilt = float(inverse_laplace(equation.substitute(rate_by_mutype) / var('E'), var('E'), var('T'), algorithm='giac').substitute(T=split_time))
                        prob_by_data_string[data_string] = probability_ilt
                        A[data_vector] = probability_ilt
                        data_queries_seen.add(data_string)
                else:
                    data_query = tuple(
                        [(data_point[mutype] if mutype in mutypes else slice(0, self.max_by_mutype[mutype] + 2)) 
                            for mutype in MUTYPES])
                    if not (data_query[1] > 0 and data_query[3] > 0):
                        if not data_string in data_queries_seen:
                            probability_ilt = float(inverse_laplace(equation.substitute(rate_by_mutype) / var('E'), var('E'), var('T'), algorithm='giac').substitute(T=split_time))
                            probability = probability_ilt - np.sum(A[data_query].flatten())  
                            prob_by_data_string[data_string] = probability
                            A[data_vector] = probability
                            data_queries_seen.add(data_string)
        composite_likelihood = np.sum((xlogy(np.sign(A), A) * A))
        print("[=] Calculated probabilities in %s seconds..." % (timer() - start_time))
        # return prob_by_data_string, composite_likelihood
        return composite_likelihood

def powerset(iterable):
    "powerset([1,2,3]) --> () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)"
    s = list(iterable)
    return iter_chain.from_iterable(combi(s, r) for r in range(len(s)+1))

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

def calculate_equation(params):
    data_idx, data_point, pathObj, max_by_mutype, base_rate, symbolic_mutation_rate = params
    mutation_equation = pathObj.infer_mutation_equation(data_point, symbolic_mutation_rate)
    equation = pathObj.multiply_with_path_equation(mutation_equation)
    return (data_idx, pathObj.path_id, equation)

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
        self.mutable_nodes_by_mutype = defaultdict(list)
        self.path_numerators = []
        self.path_denominators = []
        self.path_equation = 1

    def add_step(self, nodeObj, event_id, event_count):
        self.nodeObjs.append(nodeObj)
        self.event_ids.append(event_id)
        self.event_counts.append(event_count)

    def infer_path_quation(self, symbolic_rate, symbolic_mutation_rate):
        self.path_numerators = []
        self.path_denominators = []
        for nodeObj, event_id, event_count in self.yield_step():
            self.path_numerators.append((event_count * symbolic_rate[event_id]))
            self.path_denominators.append(
                        (nodeObj.event_counter['C_ancestor'] * symbolic_rate['C_ancestor']) + 
                        (nodeObj.event_counter['C_derived'] * symbolic_rate['C_derived']) + 
                        (nodeObj.event_counter['M'] * symbolic_rate['M']) + 
                        (nodeObj.event_counter['E'] * symbolic_rate['E']) +
                        (nodeObj.mutype_counter['hetA'] * symbolic_mutation_rate['hetA']) +
                        (nodeObj.mutype_counter['hetB'] * symbolic_mutation_rate['hetB']) +
                        (nodeObj.mutype_counter['hetAB'] * symbolic_mutation_rate['hetAB']) +
                        (nodeObj.mutype_counter['fixed'] * symbolic_mutation_rate['fixed'])
                )
        for numerator, denominator in zip(self.path_numerators, self.path_denominators):
            self.path_equation *= (numerator / denominator)

    def multiply_with_path_equation(self, mutation_equation):
        initial = 1
        for numerator, denominator in zip(self.path_numerators, self.path_denominators):
            initial *= (numerator / denominator)
        mutation_equation *= initial
        return mutation_equation

    def infer_mutation_equation(self, data_point, symbolic_mutation_rate):
        mutation_equation = 0
        for slots in prod(*[list(combi_with_replacement(self.mutable_nodes_by_mutype[mutype], data_point[mutype])) for mutype in MUTYPES]):
            mutypes_by_idx = defaultdict(list)
            for mutype, slot in zip(MUTYPES, slots):
                if slot:
                    for idx in list(slot):
                        mutypes_by_idx[idx].append(mutype)
            mutation_part = 1
            for node_idx, mutypes in mutypes_by_idx.items():
                mutype_counter = Counter(mutypes)
                node_multinomial = multinomial(mutype_counter.values())
                mutation_part *= node_multinomial
                nodeObj = self.nodeObjs[node_idx]
                denominator = self.path_denominators[node_idx]
                for mutype in mutypes:
                    mutation_part *= ((nodeObj.mutype_counter[mutype] * symbolic_mutation_rate[mutype]) / denominator)
            mutation_equation += mutation_part
        return mutation_equation

    def infer_mutable_nodes(self):
        for idx, nodeObj in enumerate(self.nodeObjs):
            for mutype, count in nodeObj.mutype_counter.items():
                if count:
                    self.mutable_nodes_by_mutype[mutype].append(idx)

    def is_compatible_with_data_point(self, data_point):
        if sum(data_point.values()) == 0:
            return True
        else:
            mutypes_in_data_point = set([mutype for mutype, count in data_point.items() if count > 0])
            if all([(True if mutype in self.mutable_nodes_by_mutype else False) for mutype in mutypes_in_data_point]):
                return True
        return False

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
        self.precision = RealField(int(args['--precision']) * log(10,2)) 
        self.user_rate = {
            'M' : Rational(float(args['--migration_rate'])),
            'C_ancestor' : Rational(1.0),
            'C_derived' : Rational(float(args['--derived_coalescence_rate'])),
            'T' : Rational(float(args['--T_var'])),
            'theta' : Rational(float(args['--mutation_rate']))
        }
        self.base_rate = {
            var('M') : Rational(float(args['--migration_rate']) / 2),
            var('C_ancestor') : Rational(1.0),
            var('C_derived') : Rational(float(args['--derived_coalescence_rate'])),
            var('E') : var('E', domain='positive')
        }
        self.split_time = Rational(float(args['--T_var']))
        self.exodus_rate = var('E', domain='positive')
        self.mutation_rate = Rational(float(args['--mutation_rate']) / 2)
        self.symbolic_mutation_rate = {
            'hetA' : var('theta_A', domain='positive'),
            'hetB' : var('theta_B', domain='positive'),
            'hetAB' : var('theta_AB', domain='positive'),
            'fixed' : var('theta_fix', domain='positive')
        }
        self.symbolic_rate = {
            'C_ancestor' : var('C_ancestor', domain='positive'), 
            'C_derived' : var('C_derived', domain='positive'), 
            'M' : var('M', domain='positive'),
            'E' : var('E', domain='positive')
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
        print("[=] Parsed %s paths from file %s in %s seconds." % (len(pathObj_by_path_id), self.path_file, timer() - start_time))
        return (pathObj_by_path_id, nodeObj_by_node_id) 

    def prepare_paths(self, pathObj_by_path_id):
        # infer path_probabilities for each pathObj (only has to be done once) ... 
        start_time = timer()
        print("[+] Preparing %s paths (path equations and mutable nodes) ..." % (len(pathObj_by_path_id)))
        for path_id, pathObj in tqdm(pathObj_by_path_id.items(), total=len(pathObj_by_path_id), desc="[%] ", ncols=200):
            pathObj.infer_mutable_nodes()
            pathObj.infer_path_quation(self.symbolic_rate, self.symbolic_mutation_rate)
        print("[=] Prepared paths in %s seconds." % (timer() - start_time))
        return pathObj_by_path_id

    def generate_equations(self, data, pathObj_by_path_id):
        start_time = timer()
        params = []
        equations_by_data_idx = defaultdict(list)
        for data_idx, data_point in enumerate(data):
            for path_id, pathObj in pathObj_by_path_id.items():
                if pathObj.is_compatible_with_data_point(data_point):
                    params.append((data_idx, data_point, pathObj, self.max_by_mutype, self.base_rate, self.symbolic_mutation_rate))
                else:
                    #equationn is , rates0 if data point can't be placed on path ...
                    equations_by_data_idx[data_idx].append(0)                            
        print("[+] Analysing %s combinations of paths and data points with %s threads..." % (len(params), self.threads))
        if self.threads == 1:
            for param in tqdm(params, total=len(params), desc="[%] ", ncols=200):
                data_idx, path_id, equation = calculate_equation(param)
                equations_by_data_idx[data_idx].append(equation)
        else:
            with poolcontext(processes=self.threads) as pool:
                with tqdm(params, total=len(params), desc="[%] ", ncols=200) as pbar:
                    for data_idx, path_id, equation in pool.imap_unordered(calculate_equation, params):
                        equations_by_data_idx[data_idx].append(equation)
                        pbar.update()
        print("[=] Analysed paths in %s seconds." % (timer() - start_time))
        return equations_by_data_idx

    def calculate_initial_probabilities(self, data, symbolic_equations_by_data_idx):
        start_time = timer()
        print("[+] Calculating PODs ...")
        
        base_rate = self.base_rate
        mutation_rate = self.mutation_rate
        split_time = self.split_time

        # symbolic_equations_by_data_idx are completely symbolic 
        equation_by_data_idx = {data_idx: sum(equation).subs(base_rate) for data_idx, equation in symbolic_equations_by_data_idx.items()}

        # Setting up datastructure
        shape = tuple(self.max_by_mutype[mutype] + 2 for mutype in MUTYPES)
        A = np.zeros(shape, np.float64)
        prob_by_data_string = OrderedDict()
        data_queries_seen = set()
        theta_var_by_mutype = {
            'hetA' : var('theta_A'),
            'fixed' : var('theta_fix'),
            'hetB' : var('theta_B'),
            'hetAB' : var('theta_AB')
        }
        # Generating sets of mutation_rate combinations for marginals
        sets_of_mutation_rates = list(powerset(MUTYPES))
        #for mutypes in tqdm(reversed(sets_of_mutation_rates), total=len(sets_of_mutation_rates), desc="[%] ", ncols=200):
        for mutypes in reversed(sets_of_mutation_rates):
            mutypes_masked = [mutype for mutype in MUTYPES if not mutype in mutypes]
            for data_idx, equation in equation_by_data_idx.items():                
                data_point = data[data_idx]
                rate_by_mutype = {}
                _data_vector = []
                _data_string = []
                for mutype in MUTYPES:
                    if mutype in mutypes:
                        rate_by_mutype[theta_var_by_mutype[mutype]] = mutation_rate
                        _data_vector.append(data_point[mutype]) 
                        _data_string.append(data_point[mutype])
                    else:
                        rate_by_mutype[theta_var_by_mutype[mutype]] = 0
                        _data_vector.append(self.max_by_mutype[mutype] + 1)
                        _data_string.append("%s+" % (self.max_by_mutype[mutype] + 1))
                    data_vector = tuple(_data_vector)
                    data_string = tuple(_data_string)
                if not mutypes_masked: 
                    if not data_string in data_queries_seen:
                        probability_ilt = inverse_laplace(equation.substitute(rate_by_mutype) / var('E'), var('E'), var('T'), algorithm='giac').substitute(T=split_time)
                        prob_by_data_string[data_string] = probability_ilt
                        A[data_vector] = probability_ilt
                        data_queries_seen.add(data_string)
                else:
                    
                    data_query = tuple(
                        [(data_point[mutype] if mutype in mutypes else slice(0, self.max_by_mutype[mutype] + 2)) 
                            for mutype in MUTYPES])
                    if not (data_query[1] > 0 and data_query[3] > 0):
                        if not data_string in data_queries_seen:
                            print("marginal", data_idx, data_point, data_query)
                            probability_ilt = inverse_laplace(equation.substitute(rate_by_mutype) / var('E'), var('E'), var('T'), algorithm='giac').substitute(T=split_time)
                            probability = probability_ilt - sum(A[data_query].flatten())  
                            prob_by_data_string[data_string] = probability
                            A[data_vector] = probability
                            print(data_vector, rate_by_mutype, A[data_vector])
                            data_queries_seen.add(data_string)
        
        composite_likelihood = self.precision(-np.sum((xlogy(np.sign(A), A) * A)))
        print("[=] Calculated PODs (L=%s) in %s seconds..." % (composite_likelihood, timer() - start_time))
        exit()
        return composite_likelihood, A

    def calculate_probabilities(self, x0, *args):
        start_time = timer()
        #print("[+] Calculating probabilities (datapoints and marginals)...")
        # unpack parameters
        PoD, data, symbolic_equations_by_data_idx = args
        split_time = Rational(x0[0])
        mutation_rate = Rational(x0[1])
        #if any([(True if x < 0.1 else False) for x in x0]):
        #    return 100

        # symbolic_equations_by_data_idx are completely symbolic 
        equation_by_data_idx = {data_idx: sum(equation).subs(self.base_rate) for data_idx, equation in symbolic_equations_by_data_idx.items()}

        # Setting up datastructure
        shape = tuple(self.max_by_mutype[mutype] + 2 for mutype in MUTYPES)
        A = np.zeros(shape, np.float64)
        prob_by_data_string = OrderedDict()
        data_queries_seen = set()
        theta_var_by_mutype = {
            'hetA' : var('theta_A'),
            'fixed' : var('theta_fix'),
            'hetB' : var('theta_B'),
            'hetAB' : var('theta_AB')
        }
        # Generating sets of mutation_rate combinations for marginals
        sets_of_mutation_rates = list(powerset(MUTYPES))
        #for mutypes in tqdm(reversed(sets_of_mutation_rates), total=len(sets_of_mutation_rates), desc="[%] ", ncols=200):
        for mutypes in reversed(sets_of_mutation_rates):
            mutypes_masked = [mutype for mutype in MUTYPES if not mutype in mutypes]
            for data_idx, equation in equation_by_data_idx.items():                
                data_point = data[data_idx]
                rate_by_mutype = {}
                _data_vector = []
                for mutype in MUTYPES:
                    if mutype in mutypes:
                        rate_by_mutype[theta_var_by_mutype[mutype]] = mutation_rate
                        _data_vector.append(data_point[mutype])
                    else:
                        rate_by_mutype[theta_var_by_mutype[mutype]] = 0
                        _data_vector.append(self.max_by_mutype[mutype] + 1)
                    data_vector = tuple(_data_vector)

                if not mutypes_masked: 
                    if not data_vector in data_queries_seen:
                        probability_ilt = inverse_laplace(equation.substitute(rate_by_mutype) / var('E'), var('E'), var('T'), algorithm='giac').substitute(T=split_time)
    
                        #print(data_string, probability_ilt)
                        A[data_vector] = probability_ilt
                        data_queries_seen.add(data_vector)
                else:
                    data_query = tuple(
                        [(data_point[mutype] if mutype in mutypes else slice(0, self.max_by_mutype[mutype] + 2)) 
                            for mutype in MUTYPES])
                    if not (data_query[1] > 0 and data_query[3] > 0):
                        if not data_vector in data_queries_seen:
                            probability_ilt = inverse_laplace(equation.substitute(rate_by_mutype) / var('E'), var('E'), var('T'), algorithm='giac').substitute(T=split_time)
                            probability = probability_ilt - sum(A[data_query].flatten())
        
                            #print(data_string, probability)
                            A[data_vector] = probability
                            data_queries_seen.add(data_vector)
        
        #print("[=] Calculated probabilities in %s seconds..." % (timer() - start_time))
        composite_likelihood = self.precision(-np.sum((xlogy(np.sign(A), A) * PoD)))
        print('L=-%s\tT=%s\ttheta=%s\t[M=%s\tC_ancestor=%s,\tC_derived=%s]' % (composite_likelihood, self.precision(split_time), self.precision(mutation_rate), self.precision(self.base_rate[var('M')]), self.base_rate[var('C_ancestor')], self.base_rate[var('C_derived')]))
        return composite_likelihood

    def write_probs(self, prob_by_data_string, data):
        header = MUTYPES + ['probability']
        rows = []
        prob_sum = 0
        for data_string, prob in prob_by_data_string.items():
            rows.append([count for count, mutype in zip(data_string, MUTYPES)] + [prob])
            try: 
                prob_sum += float(prob)
            except TypeError:
                print("[!] TypeError when summing: %s" % prob)
        rows.append(['*', '*', '*', '*',  prob_sum])
        out_f = "%s.C_anc=%s.C_der=%s.M=%s.theta=%s.T=%s.threads=%s.tsv" % (self.out_prefix, self.user_rate['C_ancestor'], self.user_rate['C_derived'], self.user_rate['M'], self.user_rate['theta'], self.split_time, self.threads)
        create_csv(out_f, header, rows, "\t")

    def write_path_equations(self, pathObj_by_path_id):
        out_f = "%s.C_anc=%s.C_der=%s.M=%s.theta=%s.T=%s.threads=%s.path_equations.txt" % (self.out_prefix, self.user_rate['C_ancestor'], self.user_rate['C_derived'], self.user_rate['M'], self.user_rate['theta'], self.split_time, self.threads)
        out_lines = []
        for path_id, pathObj in pathObj_by_path_id.items():
            out_lines.append(str(pathObj.path_equation))
        with open(out_f, 'w') as out_fh:
            out_fh.write(str(out_lines))
            out_fh.write("\n")

    def setup_data_space(self):
        data = []
        start_time = timer()
        print("[+] Generating base data points (based on MAX_MUTYPES: %s) ..." % (", ".join(["%s=%s" % (mutype, self.max_by_mutype[mutype]) for mutype in MUTYPES])))
        # works only for equal max_mutypes ...
        for i, data_point in enumerate(prod(*[range(0, self.max_by_mutype[mutype] + 1) for mutype in MUTYPES])):
            counter = Counter()
            for mutype, count in zip(MUTYPES, data_point):
                counter[mutype] = count
            # Four-gamete-test : exclude data points w/ counter['fixed'] > 0 AND counter['hetAB'] > 0
            if not (counter['fixed'] > 0 and counter['hetAB'] > 0):
                data.append(counter)
        print("[=] Generated %s base data points in %s seconds." % (len(data), timer() - start_time))
        return data

def generate_initial_simplex(boundaries, seed=0):
    # boundaries = {
    # 'a': (0, 10), 
    # 'b': (20, 30),
    # 'c': (30, 40)
    # }
    np.random.seed(seed)
    simplex = [] 
    for i in range(len(boundaries) + 1):
        vertex = []
        for parameter, (_min, _max) in boundaries.items():
            value = numpy.random.uniform([_min, _max])
            vertex.append(value)
        simplex.append(tuple(vertex))
    return simplex


            # simplex = [(a1, b2),
                       # (a2, b1),
                       # (a1, b1)]

            # random.shuffle(data)
         #       for elem in data: process(elem) 
# 
        


###############################################################################

################################### Main ######################################

def main():
    try:
        main_time = timer()
        args = docopt(__doc__)
        print(args)
        parameterObj = ParameterObj(args)
        print(parameterObj.base_rate)
        pathObj_by_path_id, nodeObj_by_node_id = parameterObj.read_paths()
        pathObj_by_path_id = parameterObj.prepare_paths(pathObj_by_path_id)
        data = parameterObj.setup_data_space()
        symbolic_equations_by_data_idx = parameterObj.generate_equations(data, pathObj_by_path_id)
        # prob_by_data_string, composite_likelihood = parameterObj.calculate_probabilities(data, symbolic_equations_by_data_idx, params)
        
        pod_composite_likelihood, pod = parameterObj.calculate_initial_probabilities(data, symbolic_equations_by_data_idx)
        # TBE_PARAMETERS_BY_MODEL = {
        #    'M': ['T', 'theta', 'M'], 
        #    'E': ['T', 'theta'], 
        #    'R': ['T', 'theta'], 
        #    'M+E': ['T', 'theta', 'M'], 
        #    'M+R': ['T', 'theta']
        #    } 
        #estimate_parameters()
        #initial_simplex = generate_initial_simplex()
        initial_simplex = [
                            (Rational(0.1), Rational(2.0)), 
                            (Rational(2.0), Rational(0.1)), 
                            #(Rational(0.1), Rational(0.1))
                            (Rational(2.0), Rational(2.0))
                            ]
        x0 = (Rational(2.0), Rational(0.1)) # T, theta
        res = minimize(
            parameterObj.calculate_probabilities, 
            x0, 
            args=(pod, data, symbolic_equations_by_data_idx), 
            method="Nelder-Mead", 
            options={
                'initial_simplex': initial_simplex, 
                'disp': True, 
                'xatol': 1e-1, 
                'fatol': 1e-5, 
                'adaptive': False}) 
        
        print(res)
        print("[+] Total runtime: %s seconds" % (timer() - main_time))
        exit()
        parameterObj.write_probs(prob_by_data_string, data)
        parameterObj.write_path_equations(pathObj_by_path_id)
        print("[+] Total runtime: %s seconds" % (timer() - main_time))
    except KeyboardInterrupt:
        stderr.write("\n[X] Interrupted by user after %s seconds!\n" % (timer() - main_time))
        exit(-1)
        
###############################################################################

if __name__ == '__main__':
    main()