#!/usr/bin/env python3
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

# stdlib
from __future__ import division, print_function
from collections import defaultdict, Counter, OrderedDict
from timeit import default_timer as timer
from sys import stderr, exit, stdout
from multiprocessing import Pool
from contextlib import contextmanager
from itertools import combinations_with_replacement as combi_with_replacement
from itertools import combinations as combi
from itertools import chain as iter_chain 
from itertools import product as prod
import subprocess
# external
from docopt import docopt
from pandas import read_csv
from tqdm import tqdm
import numpy as np
import sympy
from scipy.special import xlogy
from scipy.optimize import minimize
from pandas import DataFrame

'''
[To do]
- explicit imports !
- change linewidth pep8

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

#THETA_VAR_BY_MUTYPE = {
#            'hetA' : sympy.Symbol('theta_A'),
#            'fixed' : sympy.Symbol('theta_fix'),
#            'hetB' : sympy.Symbol('theta_B'),
#            'hetAB' : sympy.Symbol('theta_AB')
#        }

TBE_PARAMETERS_BY_MODEL = {
           'Mig': ['Time', 'theta', 'Mig'], 
           'bigL': ['Time', 'theta'], 
           'R': ['Time', 'theta'], 
           'M+E': ['Time', 'theta', 'Mig'], 
           'M+R': ['Time', 'theta']
           }

FOUR_GAMETE_VIOLATION = set(['hetAB', 'fixed'])
###############################################################################

################################### Functions #################################

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

def inverse_laplace_transform(params):
    vector, marginal_query, equation, rate_by_mutype, split_time = params
    equation_giac = str(equation).replace("**", "^")
    mutation_substitution = ",[%s]" % ",".join(["%s=%s" % (mutype, rate) for mutype, rate in rate_by_mutype.items()])
    assumptions = "assume(bigL >= 0)"
    #print()
    '''
    [ratnormal] rewrites an expression using its irreducible representation. The expression is viewed as a multivariate rational fraction with coefficients in Q (or
        Q[i]). The variables are generalized identifiers which are assumed to be algebraically independent. Unlike with normal, an algebraic extension is considered
        as a generalized identifier. Therefore ratnormal is faster but might miss some
        simplifications if the expression contains radicals or algebraically dependent transcendental functions.
    '''
    invlaplace_string = "%s; invlaplace(ratnormal(subst(%s%s)/bigL), bigL, T).subst(T, %s)" % (assumptions, equation_giac, mutation_substitution, float(split_time))
    print(invlaplace_string)
    #start_time = timer()
    try:
        process = subprocess.run(["giac", invlaplace_string], stderr=subprocess.PIPE, stdout=subprocess.PIPE, check=True, encoding='utf-8')
        #print("elapsed: %s" % (timer() - start_time))
        print(process.stderr)
        print(process.stdout)
        probability = process.stdout.rstrip("\n").split(",")[1]

    except subprocess.CalledProcessError:
        exit("[X] giac could not run.")
    #print(probability)
    return sympy.Rational(str(probability)), vector, marginal_query

@contextmanager
def poolcontext(*args, **kwargs):
    pool = Pool(*args, **kwargs)
    yield pool
    pool.terminate()

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

    def infer_path_quation(self, parameterObj):
        self.path_numerators = []
        self.path_denominators = []
        for nodeObj, event_id, event_count in self.yield_step():
            self.path_numerators.append((event_count * parameterObj.symbolic_rate[event_id]))
            self.path_denominators.append(
                        (nodeObj.event_counter['C_ancestor'] * parameterObj.symbolic_rate['C_ancestor']) + 
                        (nodeObj.event_counter['C_derived'] * parameterObj.symbolic_rate['C_derived']) + 
                        (nodeObj.event_counter['Mig'] * parameterObj.symbolic_rate['Mig']) + 
                        (nodeObj.event_counter['bigL'] * parameterObj.symbolic_rate['bigL']) +
                        (nodeObj.mutype_counter['hetA'] * parameterObj.symbolic_mutation_rate['hetA']) +
                        (nodeObj.mutype_counter['hetB'] * parameterObj.symbolic_mutation_rate['hetB']) +
                        (nodeObj.mutype_counter['hetAB'] * parameterObj.symbolic_mutation_rate['hetAB']) +
                        (nodeObj.mutype_counter['fixed'] * parameterObj.symbolic_mutation_rate['fixed'])
                )
        for numerator, denominator in zip(self.path_numerators, self.path_denominators):
            self.path_equation *= numerator / denominator

    def multiply_with_path_equation(self, mutation_equation):
        initial = 1
        for numerator, denominator in zip(self.path_numerators, self.path_denominators):
            initial *= numerator / denominator
        mutation_equation *= initial
        return mutation_equation

    def infer_mutation_equation(self, data_tuple, parameterObj):
        mutation_equation = 0
        for slots in prod(*[list(combi_with_replacement(self.mutable_nodes_by_mutype[mutype], count)) for count, mutype in zip(data_tuple, MUTYPES)]):
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
                    mutation_part *= (nodeObj.mutype_counter[mutype] * parameterObj.symbolic_mutation_rate[mutype]) / denominator
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
        #self.precision = RealField(int(args['--precision']) * log(10,2)) 
        self.user_rate = {
            'Mig' : sympy.Rational(float(args['--migration_rate'])),
            'C_ancestor' : sympy.Rational(1.0),
            'C_derived' : sympy.Rational(float(args['--derived_coalescence_rate'])),
            'Time' : sympy.Rational(float(args['--T_var'])),
            'theta' : sympy.Rational(float(args['--mutation_rate']))
        }
        self.C_ancestor = sympy.Symbol('C_ancestor', positive=True)
        self.C_derived = sympy.Symbol('C_derived', positive=True)
        self.M = sympy.Symbol('Mig', positive=True)
        self.bigL = sympy.Symbol('bigL', positive=True)
        self.T = sympy.Symbol('Time', positive=True)
        self.theta_A = sympy.Symbol('theta_A', positive=True)
        self.theta_B = sympy.Symbol('theta_B', positive=True)
        self.theta_fix = sympy.Symbol('theta_fix', positive=True)
        self.theta_AB = sympy.Symbol('theta_AB', positive=True)
        self.base_rate = {
            self.M : sympy.Rational(str(float(args['--migration_rate']) / 2)),
            self.C_ancestor : sympy.Rational(str(1.0)),
            self.C_derived : sympy.Rational(str(float(args['--derived_coalescence_rate'])))
        }
        self.split_time = sympy.Rational(str(float(args['--T_var'])))
        self.mutation_rate = sympy.Rational(str(float(args['--mutation_rate'])/ 2))
        self.symbolic_mutation_rate = {
            'hetA' : self.theta_A,
            'hetB' : self.theta_B,
            'hetAB' : self.theta_AB,
            'fixed' : self.theta_fix
        }

        self.symbolic_rate = {
            'C_ancestor' : self.C_ancestor,
            'C_derived' : self.C_derived, 
            'Mig' : self.M,
            'bigL' : self.bigL
        }
        self.data = []

    def generate_data_space(self):
        start_time = timer()
        print("[+] Generating base data points (based on MAX_MUTYPES: %s) ..." % (", ".join(["%s=%s" % (mutype, self.max_by_mutype[mutype]) for mutype in MUTYPES])))
        # works only for equal max_mutypes ...
        for i, data_point in enumerate(prod(*[range(0, self.max_by_mutype[mutype] + 1) for mutype in MUTYPES])):
            counter = Counter()
            for mutype, count in zip(MUTYPES, data_point):
                counter[mutype] = count
            # Four-gamete-test : exclude data points w/ counter['fixed'] > 0 AND counter['hetAB'] > 0
            if not (counter['fixed'] > 0 and counter['hetAB'] > 0):
                self.data.append(counter)
        print("[=] Generated %s base data points in %s seconds." % (len(self.data), timer() - start_time))

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
        out_f = "%s.C_anc=%s.C_der=%s.M=%s.theta=%s.T=%s.threads=%s.tsv" % (self.out_prefix, self.user_rate['C_ancestor'], self.user_rate['C_derived'], self.user_rate['Mig'], self.user_rate['theta'], self.split_time, self.threads)
        create_csv(out_f, header, rows, "\t")

    def write_path_equations(self, pathObj_by_path_id):
        out_f = "%s.C_anc=%s.C_der=%s.M=%s.theta=%s.T=%s.threads=%s.path_equations.txt" % (self.out_prefix, self.user_rate['C_ancestor'], self.user_rate['C_derived'], self.user_rate['Mig'], self.user_rate['theta'], self.split_time, self.threads)
        out_lines = []
        for path_id, pathObj in pathObj_by_path_id.items():
            out_lines.append(str(pathObj.path_equation))
        with open(out_f, 'w') as out_fh:
            out_fh.write(str(out_lines))
            out_fh.write("\n")

def infer_composite_likelihood(x0, *args):
    # unpack parameters
    start_time = timer()
    symbolic_equations_by_data_tuple, pod, parameterObj = args
    split_time, mutation_rate = parameterObj.split_time, parameterObj.mutation_rate
    if not x0 is None:
        split_time = sympy.Rational(str(x0[0]))
        mutation_rate = sympy.Rational(str(x0[1]))
    # symbolic_equations_by_data_idx are completely symbolic 
    # print([(type(rate), type(number)) for rate, number in parameterObj.base_rate.items()])
    # print(symbolic_equations_by_data_tuple[(0,0,0,0)])
    # for data_tuple, equation_l in symbolic_equations_by_data_tuple.items():
    #     sumeq = 0
    #     for eq in equation_l:
    #         print([(rate, number) for rate, number in parameterObj.base_rate.items()])
    #         print("before", eq, type(eq))
    #         try:
    #             subeq = eq.subs(parameterObj.base_rate)
    #             print("SUB", subeq)
    #             sumeq += subeq 
    #         except:
    #             pass
            
    #     print("SUM", sumeq)
    # exit()
    # start = timer()
    # equation_by_data_tuple = {
    #     data_tuple: sum(equation).subs(parameterObj.base_rate, simultaneous=True)
    #     for data_tuple, equation in symbolic_equations_by_data_tuple.items()
    #     }
    # print("1: %s" % (timer() - start))
    equation_by_data_tuple = {
        data_tuple: sum(equation).xreplace(parameterObj.base_rate)
        for data_tuple, equation in symbolic_equations_by_data_tuple.items()
        }
    # start = timer()
    # equation_by_data_tuple = {
    #     data_tuple: sum(equation).evalf(subs=parameterObj.base_rate)
    #     for data_tuple, equation in symbolic_equations_by_data_tuple.items()
    #     }
    # print("3: %s" % (timer() - start))
    # Setting up datastructure
    shape = tuple(parameterObj.max_by_mutype[mutype] + 2 for mutype in MUTYPES)
    probability_matrix = np.zeros(shape, np.float64)
    # Generating sets of mutation_rate combinations for marginals
    #marginal_queries = []
    #rates_by_mutype_by_vector = {}
    vector_seen = set([])
    #equations_by_vector = defaultdict(list)
    for mutypes in reversed(list(powerset(MUTYPES))):
        mutypes_masked = set([mutype for mutype in MUTYPES if not mutype in mutypes])
        if not FOUR_GAMETE_VIOLATION.issubset(mutypes_masked):
            equations_by_vector = defaultdict(list)
            marginal_queries = []
            vectors = []
            rates_by_mutype_by_vector = {}
            for data_tuple in sorted(equation_by_data_tuple.keys()):               
                vector = tuple([count if not mutype in mutypes_masked else parameterObj.max_by_mutype[mutype] + 1 for count, mutype in zip(data_tuple, MUTYPES)])
                if not vector in vector_seen:
                    vector_seen.add(vector)
                    if not len([mutype for mutype, count in zip(MUTYPES, vector) if mutype in FOUR_GAMETE_VIOLATION and count > 0]) > 1: 
                        vectors.append(vector)
                        rates_by_mutype_by_vector[vector] = {parameterObj.symbolic_mutation_rate[mutype]: (mutation_rate if mutype in mutypes else 0) for mutype in MUTYPES}
                        if not mutypes_masked:
                            marginal_query = None
                        else:
                            marginal_query = tuple(
                                [(data_tuple[idx] if not mutype in mutypes_masked else slice(0, parameterObj.max_by_mutype[mutype] + 2)) 
                                    for idx, mutype in enumerate(MUTYPES)])
                        marginal_queries.append(marginal_query)
                        equations_by_vector[vector] = equation_by_data_tuple[data_tuple]
            # can be parallelised 
            params = [(vector, marginal_query, equations_by_vector[vector], rates_by_mutype_by_vector[vector], split_time) for vector, marginal_query in zip(vectors, marginal_queries)]
            #print(params)
            if parameterObj.threads < 2:
                for param in params:
                    probability, vector, marginal_query = inverse_laplace_transform(param)
                    if marginal_query:
                        probability -= sum(probability_matrix[marginal_query].flatten())
                    probability_matrix[vector] = probability
            else:
                with poolcontext(processes=parameterObj.threads) as pool:
                    for probability, vector, marginal_query in pool.imap_unordered(inverse_laplace_transform, params):
                        if marginal_query:
                            probability -= sum(probability_matrix[marginal_query].flatten())
                        probability_matrix[vector] = probability                            

    ##print("setup", round(timer() - start_time, 2)) # setup just takes : 0.5s
    #start_time = timer()
    #for vector, marginal_query in zip(vectors, marginal_queries):
    #    # eq = sympy.apart(equation.xreplace(rate_by_mutype) / parameterObj.bigL)
    #    probability = inverse_laplace_transform(equations_by_vector[vector], rates_by_mutype_by_vector[vector], split_time)
    #    if marginal_query:
    #        probability -= sum(probability_matrix[marginal_query].flatten())
    #    #print(vector, sympy.N(probability, 12), "%ss" % round(timer() - timeilt, 2))
    #    probability_matrix[vector] = probability
    if pod is None:
        composite_likelihood = -np.sum((xlogy(np.sign(probability_matrix), probability_matrix) * probability_matrix))
        print('[+] L=-%s\tT=%s\ttheta=%s\t' % (composite_likelihood, split_time, mutation_rate))
        return (probability_matrix, composite_likelihood)
    composite_likelihood = -np.sum((xlogy(np.sign(probability_matrix), probability_matrix) * pod))
    print(" " * 100, end='\r'),
    print("[O] L=-%s \t T=%s \t theta=%s \t iteration=%ss" % (
        composite_likelihood, 
        round(split_time, 3), 
        round(mutation_rate, 3),
        round(timer() - start_time, 2)
        ), end='\r')
    stdout.flush()
    return composite_likelihood

def calculate_pods(symbolic_equations_by_data_tuple, parameterObj):
    start_time = timer()
    print("[+] Calculating PODs ...")
    pod, composite_likelihood = infer_composite_likelihood(None, symbolic_equations_by_data_tuple, None, parameterObj)
    print("[=] Calculated PODs (L=%s) in %s seconds..." % (composite_likelihood, timer() - start_time))
    return pod

def prepare_paths(parameterObj):
    # infer path_probabilities for each pathObj (only has to be done once) ... 
    start_time = timer()
    pathObj_by_path_id = read_paths(parameterObj)
    print("[+] Preparing %s paths (path equations and mutable nodes) ..." % (len(pathObj_by_path_id)))
    for path_id, pathObj in tqdm(pathObj_by_path_id.items(), total=len(pathObj_by_path_id), desc="[%] ", ncols=200):
        pathObj.infer_mutable_nodes()
        pathObj.infer_path_quation(parameterObj)
    print("[=] Prepared paths in %s seconds." % (timer() - start_time))
    return pathObj_by_path_id

def read_paths(parameterObj):
    start_time = timer()
    print("[+] Reading paths in file %s..." % parameterObj.path_file)
    paths_df = read_csv(\
        parameterObj.path_file, \
        sep="\t", \
        )
    pathObj_by_path_id = {}
    nodeObj_by_node_id = {}
    for path_id, node_id, event_id, event_count, C_ancestor, C_derived, Ms, Es, hetA, fixed, hetB, hetAB in tqdm(paths_df.values.tolist(), total=len(paths_df.index), desc="[%] ", ncols=200):
        if not event_id == 'LCA':
            if event_id == 'E':
                event_id = 'bigL'
            if not node_id in nodeObj_by_node_id:
                mutype_counter = Counter({'hetA': hetA, 'fixed': fixed, 'hetB': hetB, 'hetAB': hetAB})
                event_counter = Counter({'C_ancestor': C_ancestor, 'C_derived': C_derived, 'Mig': Ms, 'bigL': Es})
                nodeObj_by_node_id[node_id] = NodeObj(node_id, mutype_counter, event_counter)
            # paths
            if not path_id in pathObj_by_path_id:
                pathObj_by_path_id[path_id] = PathObj(path_id)
            pathObj_by_path_id[path_id].add_step(nodeObj_by_node_id[node_id], event_id, event_count)
    print("[=] Parsed %s paths from file %s in %s seconds." % (len(pathObj_by_path_id), parameterObj.path_file, timer() - start_time))
    return pathObj_by_path_id         

def generate_equations(pathObj_by_path_id, parameterObj):
    start_time = timer()
    args = []
    equations_by_data_tuple = defaultdict(list)
    for data_point in parameterObj.data:
        data_tuple = tuple([data_point[mutype] for mutype in MUTYPES])
        for path_id, pathObj in pathObj_by_path_id.items():
            if pathObj.is_compatible_with_data_point(data_point):
                args.append((  
                    data_tuple, 
                    pathObj, 
                    parameterObj))
            else:
                #equation is 0 if data point can't be placed on path ...
                equations_by_data_tuple[data_tuple].append(0)                            
    print("[+] Analysing %s combinations of paths and data points with %s threads..." % (len(args), parameterObj.threads))
    if parameterObj.threads == 1:
        for param in tqdm(args, total=len(args), desc="[%] ", ncols=200):
            data_tuple, equation = build_equation(param)
            equations_by_data_tuple[data_tuple].append(equation)
    else:
        with poolcontext(processes=parameterObj.threads) as pool:
            with tqdm(args, total=len(args), desc="[%] ", ncols=200) as pbar:
                for data_tuple, equation in pool.imap_unordered(build_equation, args):
                    equations_by_data_tuple[data_tuple].append(equation)
                    pbar.update()
    print("[=] Analysed paths in %s seconds." % (timer() - start_time))
    return equations_by_data_tuple

def build_equation(args):
    data_tuple, pathObj, parameterObj = args
    mutation_equation = pathObj.infer_mutation_equation(data_tuple, parameterObj)
    equation = pathObj.multiply_with_path_equation(mutation_equation)
    return data_tuple, equation

def generate_initial_simplex(boundaries, seed):
    np.random.seed(seed)
    simplex = [] 
    for i in range(len(boundaries) + 1):
        vertex = []
        for parameter, (minval, maxval) in boundaries.items():
            value = np.random.uniform(minval, maxval, 1)
            vertex.append(sympy.Rational(str(value[0])))
        simplex.append(tuple(vertex))
    return simplex

def estimate_parameters(symbolic_equations_by_data_tuple, boundaries, pod, parameterObj, seed):
    print("[+] Optimising parameters: %s ..." % (", ".join(boundaries.keys())))    
    start_time = timer()
    initial_simplex = generate_initial_simplex(boundaries, seed)
    res = minimize(
        infer_composite_likelihood, 
        (0, 0), 
        args=(symbolic_equations_by_data_tuple, pod, parameterObj), 
        method="Nelder-Mead", 
        options={
            'initial_simplex': initial_simplex, 
            'maxfev' : 200,
            'disp': False, 
            'xatol': 1e-2, 
            'fatol': 1e-2, # needs to be scaled by number of blocks
            'adaptive': False})
    print()
    if res.success:
        estimated_parameters = OrderedDict({key: value for (key, _), value in zip(boundaries.items(), res.x)})
        estimated_parameters_string = ", ".join(["%s=%s" % (key, round(value, 1)) for key, value in estimated_parameters.items()])
        print("[+] Parameters estimated in %ss using %s iterations (Composite Likelihood = -%s): %s" % (timer() - start_time, res.nit, res.fun, estimated_parameters_string))
    else:
        print("[-] No covergence reached after %s iterations (%ss elapsed)" % (res.nit, timer() - start_time))
    return estimated_parameters
###############################################################################

################################### Main ######################################

def main():
    try:
        main_time = timer()
        args = docopt(__doc__)
        print(args)
        parameterObj = ParameterObj(args)
        print(parameterObj.base_rate)
        parameterObj.generate_data_space()
        pathObj_by_path_id = prepare_paths(parameterObj)
        symbolic_equations_by_data_tuple = generate_equations(pathObj_by_path_id, parameterObj)
        pod = calculate_pods(symbolic_equations_by_data_tuple, parameterObj)
        boundaries = OrderedDict({
            'Time' : (0.0, 4.0),
            'theta' : (0.0, 1.0)
        })
        estimate_parameters(symbolic_equations_by_data_tuple, boundaries, pod, parameterObj, 12345)
        #parameterObj.write_probs(prob_by_data_string, data)
        #parameterObj.write_path_equations(pathObj_by_path_id)
        print("[+] Total runtime: %s seconds" % (timer() - main_time))
    except KeyboardInterrupt:
        stderr.write("\n[X] Interrupted by user after %s seconds!\n" % (timer() - main_time))
        exit(-1)
        
###############################################################################

if __name__ == '__main__':
    main()