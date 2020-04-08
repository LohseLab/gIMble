#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""usage: gIMble inference                  -p FILE [-h|--help]
                                            [-A FLOAT] [-D FLOAT] [-M FLOAT] [-m FLOAT] [-T FLOAT] 
                                            [--A&B_coalescence_rate]
                                            [-t INT] [-P INT] 

    Options:
        -h --help                                   show this
        -p, --paths FILE                            Paths to analyse
        -A, --A_coalescence_rate FLOAT              Coalescence rate in [ancestor] pop [default: 1.0]
        -D, --B_coalescence_rate FLOAT              Coalescence rate in {derived} pop [default: 1.0]
        --A&B_coalescence_rate                       Coalescence rate in {join=ed} pop [default: 1.0]
        -M, --migration_rate FLOAT                  Migration rate per generation [default: 2.35]
        -m, --mutation_rate FLOAT                   Mutation rate/lineage [default: 0.6]
        -T, --T_var FLOAT                           T [default: 1.4]
        -P, --precision INT                         Floating point precision of probabilities [default: 30]
        -t, --threads INT                           Threads [default: 1]
        
"""

from timeit import default_timer as timer
from docopt import docopt
from tqdm import tqdm
import numpy as np
import pandas as pd
import collections
import lib.gimble
import lib.functions
import itertools
import sage.all
import scipy
'''
[Truths]
- All things concerning inference are done here for now. Partition code into input/gimble.py later on ...

[To Do]
    - Events
        C_A
        C_B
        C_AB
        J_AB
        M_AB
        M_BA
'''
#mutype_by_lineage = {
#    ('aa') : 'fixed',
#    ('bb') : 'fixed',
#    ('a') : 'hetA',
#    ('abb') : 'hetA',
#    ('aab') : 'hetB',
#    ('b') : 'hetB',
#    ('ab') : 'hetAB'
#    }

MUTYPES = ['hetA', 'fixed', 'hetB', 'hetAB']

THETA_VAR_BY_MUTYPE = {
            'hetA' : sage.all.var('theta_A'),
            'fixed' : sage.all.var('theta_fix'),
            'hetB' : sage.all.var('theta_B'),
            'hetAB' : sage.all.var('theta_AB')
        }

TBE_PARAMETERS_BY_MODEL = {
           'M': ['T', 'theta', 'M'], 
           'E': ['T', 'theta'], 
           'R': ['T', 'theta'], 
           'M+E': ['T', 'theta', 'M'], 
           'M+R': ['T', 'theta']
           }

FOUR_GAMETE_VIOLATION = set(['hetAB', 'fixed'])

def powerset(iterable):
    "powerset([1,2,3]) --> () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)"
    s = list(iterable)
    return itertools.chain.from_iterable(itertools.combinations(s, r) for r in range(len(s)+1))

def infer_composite_likelihood(x0, *args):
    # unpack parameters
    start_time = timer()
    symbolic_equations_by_data_tuple, pod, parameterObj = args
    split_time, mutation_rate = parameterObj.split_time, parameterObj.mutation_rate
    if not x0 is None:
        split_time = sage.all.Rational(x0[0])
        mutation_rate = sage.all.Rational(x0[1])
    equation_by_data_tuple = {data_tuple: sum(equation).subs(parameterObj.base_rate) for data_tuple, equation in symbolic_equations_by_data_tuple.items()}
 
    # Setting up datastructure
    shape = tuple(parameterObj.max_by_mutype[mutype] + 2 for mutype in MUTYPES)
    probability_matrix = np.zeros(shape, np.float64)
    # Generating sets of mutation_rate combinations for marginals
    marginal_queries = []
    rates_by_mutype_by_vector = {}
    vectors = []
    vector_seen = set([])
    equations_by_vector = collections.defaultdict(list)
    for mutypes in reversed(list(powerset(MUTYPES))):
        mutypes_masked = set([mutype for mutype in MUTYPES if not mutype in mutypes])
        if not mutypes_masked == FOUR_GAMETE_VIOLATION:
            for data_tuple in sorted(equation_by_data_tuple.keys()):                
                vector = tuple([count if not mutype in mutypes_masked else parameterObj.max_by_mutype[mutype] + 1 for count, mutype in zip(data_tuple, MUTYPES)])
                if not vector in vector_seen:
                    vector_seen.add(vector)
                    if not len([mutype for mutype, count in zip(MUTYPES, vector) if mutype in FOUR_GAMETE_VIOLATION and count > 0]) > 1: 
                        vectors.append(vector)
                        rates_by_mutype_by_vector[vector] = {THETA_VAR_BY_MUTYPE[mutype]: (mutation_rate if mutype in mutypes else 0) for mutype in MUTYPES}
                        if not mutypes_masked:
                            marginal_query = None
                        else:
                            marginal_query = tuple(
                                [(data_tuple[idx] if not mutype in mutypes_masked else slice(0, parameterObj.max_by_mutype[mutype] + 2)) 
                                    for idx, mutype in enumerate(MUTYPES)])
                        marginal_queries.append(marginal_query)
                        equations_by_vector[vector] = equation_by_data_tuple[data_tuple]
    for vector, marginal_query in zip(vectors, marginal_queries):
        rate_by_mutype = rates_by_mutype_by_vector[vector]
        equation = equations_by_vector[vector]
        if marginal_query:
            probability = sage.all.inverse_laplace(equation.substitute(rate_by_mutype) / sage.all.var('E'), sage.all.var('E'), sage.all.var('T'), algorithm='giac').substitute(T=split_time) - sum(probability_matrix[marginal_query].flatten())
        else:
            probability = sage.all.inverse_laplace(equation.substitute(rate_by_mutype) / sage.all.var('E'), sage.all.var('E'), sage.all.var('T'), algorithm='giac').substitute(T=split_time)
        probability_matrix[vector] = probability
    if pod is None:
        composite_likelihood = -np.sum((scipy.special.xlogy(np.sign(probability_matrix), probability_matrix) * probability_matrix))
        print('[+] L=-%s\tT=%s\ttheta=%s\t' % (composite_likelihood, parameterObj.precision(split_time), parameterObj.precision(mutation_rate)))
        return (probability_matrix, composite_likelihood)
    composite_likelihood = -np.sum((scipy.special.xlogy(np.sign(probability_matrix), probability_matrix) * pod))
    #print(" " * 100, end='\r')
    print("[O] L=-%s \t T=%s \t theta=%s \t iteration=%ss" % (
        composite_likelihood, 
        round(split_time, 3), 
        round(mutation_rate, 3),
        round(timer() - start_time, 2)
        ))
    #stdout.flush()
    return composite_likelihood

def calculate_pods(symbolic_equations_by_data_tuple, parameterObj):
    start_time = timer()
    print("[+] Calculating PODs ...")
    pod, composite_likelihood = infer_composite_likelihood(None, symbolic_equations_by_data_tuple, None, parameterObj)
    print("[=] Calculated PODs (L=%s) in %s seconds..." % (composite_likelihood, timer() - start_time))
    return pod

def generate_equations(pathObj_by_path_id, parameterObj):
    start_time = timer()
    args = []
    equations_by_data_tuple = collections.defaultdict(list)
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
    #else:
#    #    with poolcontext(processes=parameterObj.threads) as pool:
    #        with tqdm(args, total=len(args), desc="[%] ", ncols=200) as pbar:
    #            for data_tuple, equation in pool.imap_unordered(build_equation, args):
    #                equations_by_data_tuple[data_tuple].append(equation)
    #                pbar.update()
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
            vertex.append(sage.all.Rational(float(value)))
        simplex.append(tuple(vertex))
    return simplex

def estimate_parameters(symbolic_equations_by_data_tuple, boundaries, pod, parameterObj, seed):
    print("[+] Optimising parameters: %s ..." % (", ".join(boundaries.keys())))    
    start_time = timer()
    initial_simplex = generate_initial_simplex(boundaries, seed)
    res = scipy.optimize.minimize(
        infer_composite_likelihood, 
        (0, 0), 
        args=(symbolic_equations_by_data_tuple, pod, parameterObj), 
        method="Nelder-Mead", 
        options={
            'initial_simplex': initial_simplex, 
            'maxfev' : 200,
            'disp': False, 
            'xatol': 1e-2, 
            'fatol': 1e-6, 
            'adaptive': False})
    print()
    if res.success:
        estimated_parameters = collections.OrderedDict({key: value for (key, _), value in zip(boundaries.items(), res.x)})
        estimated_parameters_string = ", ".join(["=".join([key, str(sage.all.numerical_approx(value, digits=1))]) for key, value in estimated_parameters.items()])
        print("[+] Parameters estimated in %ss using %s iterations (Composite Likelihood = -%s): %s" % (timer() - start_time, res.nit, res.fun, estimated_parameters_string))
    else:
        print("[-] No covergence reached after %s iterations (%ss elapsed)" % (res.nit, timer() - start_time))
    return estimated_parameters

def multinomial(lst):
    '''
    https://stackoverflow.com/questions/46374185/does-python-have-a-function-which-computes-multinomial-coefficients
    '''
    res, i = 1, 1
    for a in lst:
        for j in range(1, a + 1):
            res *= i
            res //= j
            i += 1
    return res

def prepare_paths(parameterObj):
    # infer path_probabilities for each pathObj (only has to be done once) ... 
    start_time = timer()
    pathObj_by_path_id = read_paths(parameterObj)
    print("[+] Preparing %s paths (path equations and mutable nodes) ..." % (len(pathObj_by_path_id)))
    for path_id, pathObj in tqdm(pathObj_by_path_id.items(), total=len(pathObj_by_path_id), desc="[%] ", ncols=200):
        print(path_id)
        print(pathObj)
        pathObj.infer_mutable_nodes()
        pathObj.infer_path_quation(parameterObj)
    print("[=] Prepared paths in %s seconds." % (timer() - start_time))
    return pathObj_by_path_id

def read_paths(parameterObj):
    '''
    when reading paths 
    - mutation labels get summarised (for now by removing integers)
        
        a1        : 'hetA'          
        a1a1      : 'fixed'             
        a1a1b1    : 'hetB'              
        a1a1b1b1  : None                  
        a1b1      : 'hetAB'              
        a1b1b1    : 'hetA'              
        b1        : 'hetB'          
        b1b1      : 'fixed'              
            
    '''
    mutype_by_lineage = {
        ('a1a1') : 'fixed',
        ('b1b1') : 'fixed',
        ('a1') : 'hetA',
        ('a1b1b1') : 'hetA',
        ('a1a1b1') : 'hetB',
        ('b1') : 'hetB',
        ('a1b1') : 'hetAB'
    }
    start_time = timer()
    print("[+] Reading paths in file %s..." % parameterObj.path_file)
    paths_df = pd.read_csv(\
        parameterObj.path_file, \
        sep="\t", \
        )
    pathObj_by_path_id = {}
    nodeObj_by_node_id = {}
    events = [event for event in paths_df.columns[7:-8]]
    print("events", events) 
    for line in tqdm(paths_df.values.tolist(), total=len(paths_df.index), desc="[%] ", ncols=200):
        value_by_key = {key: value for key, value in zip(paths_df.columns, line)}
        if not value_by_key['event'] == 'LCA':
            node_id = value_by_key['node_id']
            path_id = value_by_key['path_idx']
            if not node_id in nodeObj_by_node_id:
                mutype_counter = collections.Counter()
                event_counter = collections.Counter()
                for lineage, mutype in mutype_by_lineage.items():
                    mutype_counter[mutype] += value_by_key[lineage]
                for event in events:
                    event_counter[event] += value_by_key[event]
                #mutype_counter = collections.Counter({'hetA': hetA, 'fixed': fixed, 'hetB': hetB, 'hetAB': hetAB})
                #event_counter = collections.Counter({'C_ancestor': C_ancestor, 'C_derived': C_derived, 'M': Ms, 'E': Es})
                nodeObj_by_node_id[node_id] = NodeObj(node_id, mutype_counter, event_counter)
            # paths
            if not path_id in pathObj_by_path_id:
                pathObj_by_path_id[path_id] = PathObj(path_id)
            pathObj_by_path_id[path_id].add_step(nodeObj_by_node_id[node_id], value_by_key['event'], value_by_key['count'])
    print("[=] Parsed %s paths from file %s in %s seconds." % (len(pathObj_by_path_id), parameterObj.path_file, timer() - start_time))
    return pathObj_by_path_id   

class ParameterObj(object):
    def __init__(self, args):
        self.path_file = args['--paths']
        self.out_prefix = ".".join(args['--paths'].split(".")[0:-1])
        self.threads = int(args['--threads'])
        self.max_by_mutype = {'hetA': 2, 'hetB': 2, 'hetAB': 2, 'fixed': 2}
        self.precision = sage.all.RealField(int(args['--precision']) * sage.all.log(10,2)) 
        self.user_rate = {
            'M_AB' : sage.all.Rational(float(args['--migration_rate'])),
            'M_BA' : sage.all.Rational(float(args['--migration_rate'])),
            'C_A' : sage.all.Rational(1.0),
            'C_B' : sage.all.Rational(float(args['--B_coalescence_rate'])),
            'C_AB' : sage.all.Rational(float(args['--A&B_coalescence_rate'])),
            'T' : sage.all.Rational(float(args['--T_var'])),
            'theta' : sage.all.Rational(float(args['--mutation_rate']))
        }
        self.base_rate = {
            sage.all.var('M') : sage.all.Rational(float(args['--migration_rate']) / 2),
            sage.all.var('M_AB') : sage.all.Rational(float(args['--migration_rate']) / 2),
            sage.all.var('M_BA') : sage.all.Rational(float(args['--migration_rate']) / 2),
            sage.all.var('C_A') : sage.all.Rational(1.0),
            sage.all.var('C_B') : sage.all.Rational(float(args['--B_coalescence_rate'])),
            sage.all.var('C_AB') : sage.all.Rational(float(args['--A&B_coalescence_rate'])),
            sage.all.var('E') : sage.all.var('E') #domain='positive')
        }
        self.split_time = sage.all.Rational(float(args['--T_var']))
        self.exodus_rate = sage.all.var('E'), #domain='positive')
        self.mutation_rate = sage.all.Rational(float(args['--mutation_rate']) / 2)
        self.symbolic_mutation_rate = {
            'hetA' : sage.all.var('theta_A'), #domain='positive'),
            'hetB' : sage.all.var('theta_B'), #domain='positive'),
            'hetAB' : sage.all.var('theta_AB'), #domain='positive'),
            'fixed' : sage.all.var('theta_fix') #domain='positive')
        }
        self.symbolic_rate = {
            'C_A' : sage.all.var('C_A'), #domain='positive'), 
            'C_B' : sage.all.var('C_B'), #domain='positive'), 
            'C_AB' : sage.all.var('C_AB'), #domain='positive'), 
            'M_AB' : sage.all.var('M_AB'),
            'M_BA' : sage.all.var('M_BA'),
            'M' : sage.all.var('M'), #domain='positive'),
            'A&B' : sage.all.var('AB') #domain='positive')
        }
        self.data = []

    def generate_data_space(self):
        start_time = timer()
        print("[+] Generating base data points (based on MAX_MUTYPES: %s) ..." % (", ".join(["%s=%s" % (mutype, self.max_by_mutype[mutype]) for mutype in MUTYPES])))
        # works only for equal max_mutypes ...
        for i, data_point in enumerate(itertools.product(*[range(0, self.max_by_mutype[mutype] + 1) for mutype in MUTYPES])):
            counter = collections.Counter()
            for mutype, count in zip(MUTYPES, data_point):
                counter[mutype] = count
            # Four-gamete-test : exclude data points w/ counter['fixed'] > 0 AND counter['hetAB'] > 0
            if not (counter['fixed'] > 0 and counter['hetAB'] > 0):
                self.data.append(counter)
        print("[=] Generated %s base data points in %s seconds." % (len(self.data), timer() - start_time))
#
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
        lib.functions.create_csv(out_f, header, rows, "\t")
#
    def write_path_equations(self, pathObj_by_path_id):
        out_f = "%s.C_anc=%s.C_der=%s.M=%s.theta=%s.T=%s.threads=%s.path_equations.txt" % (self.out_prefix, self.user_rate['C_ancestor'], self.user_rate['C_derived'], self.user_rate['M'], self.user_rate['theta'], self.split_time, self.threads)
        out_lines = []
        for path_id, pathObj in pathObj_by_path_id.items():
            out_lines.append(str(pathObj.path_equation))
        with open(out_f, 'w') as out_fh:
            out_fh.write(str(out_lines))
            out_fh.write("\n")       

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
        self.mutable_nodes_by_mutype = collections.defaultdict(list)
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
                        (nodeObj.event_counter['C_A'] * parameterObj.symbolic_rate['C_A']) + 
                        (nodeObj.event_counter['C_B'] * parameterObj.symbolic_rate['C_B']) + 
                        (nodeObj.event_counter['C_AB'] * parameterObj.symbolic_rate['C_AB']) + 
                        (nodeObj.event_counter['M_AB'] * parameterObj.symbolic_rate['M_AB']) + 
                        (nodeObj.event_counter['M_BA'] * parameterObj.symbolic_rate['M_BA']) + 
                        (nodeObj.event_counter['A&B'] * parameterObj.symbolic_rate['A&B']) +
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
        for slots in itertools.product(*[list(itertools.combinations_with_replacement(self.mutable_nodes_by_mutype[mutype], count)) for count, mutype in zip(data_tuple, MUTYPES)]):
            mutypes_by_idx = collections.defaultdict(list)
            for mutype, slot in zip(MUTYPES, slots):
                if slot:
                    for idx in list(slot):
                        mutypes_by_idx[idx].append(mutype)
            mutation_part = 1
            for node_idx, mutypes in mutypes_by_idx.items():
                mutype_counter = collections.Counter(mutypes)
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

#class ParameterObj(object):
#    def __init__(self, args):
#        print(args)
#        self.zarr = args['--zarr']
#        self.population_sizes_by_pop_id = self.get_population_sizes(args['--population_sizes'])
#        self.population_names_by_pop_id = self.get_population_names(args['--population_names'])
#        self.model = args['--model']
#        
#        self.migration = float(args['--migration'])
#        self.theta = float(args['--theta'])
#        self.time = float(args['--time'])
#        self.kmax = int(args['--kmax'])

def main(run_params):
    try:
        start_time = timer()    
        args = docopt(__doc__)
        print(args)
        parameterObj = ParameterObj(args)
        print(parameterObj.base_rate)
        parameterObj.generate_data_space()
        pathObj_by_path_id = prepare_paths(parameterObj)
        #Start
        symbolic_equations_by_data_tuple = generate_equations(pathObj_by_path_id, parameterObj)
        pod = calculate_pods(symbolic_equations_by_data_tuple, parameterObj)
        print(pod)
        boundaries = collections.OrderedDict({
            'T' : (0.0, 4.0),
            'theta' : (0.0, 1.0)
        })
        estimated_parameters = estimate_parameters(symbolic_equations_by_data_tuple, boundaries, pod, parameterObj, 12345)
        print(estimated_parameters)
        #parameterObj.write_probs(prob_by_data_string, data)
        #parameterObj.write_path_equations(pathObj_by_path_id)
        #print("[+] Total runtime: %s seconds" % (timer() - main_time))
        print("[*] Total runtime: %.3fs" % (timer() - start_time))
    except KeyboardInterrupt:
        print("\n[X] Interrupted by user after %s seconds!\n" % (timer() - start_time))
        exit(-1)