# stdlib
from __future__ import division, print_function
import sys
import collections
from timeit import default_timer as timer
from multiprocessing import Pool
from contextlib import contextmanager
import itertools
import subprocess
# external
import pandas as pd
from tqdm import tqdm
import numpy as np
import sympy
from scipy.special import xlogy
from scipy.optimize import minimize
from lib.classes import EntityCollection
from lib.functions import format_percentage
from lib.functions import check_path, check_file, check_prefix, format_count
import pathlib


################################### CONSTANTS #################################

MUTYPES = ['hetA', 'fixed', 'hetB', 'hetAB']
FOUR_GAMETE_VIOLATION_IDX = [idx for idx, mutype in enumerate(MUTYPES) if mutype in set(['fixed', 'hetAB'])]
FOUR_GAMETE_VIOLATION = set(['fixed', 'hetAB'])
#THETA_VAR_BY_MUTYPE = {
#            'hetA' : sympy.Symbol('theta_A'),
#            'fixed' : sympy.Symbol('theta_fix'),
#            'hetB' : sympy.Symbol('theta_B'),
#            'hetAB' : sympy.Symbol('theta_AB')
#        }

PARAMETERS_BY_MODEL_NAME = {
           #'model.M.txt': ['theta', 'Mig'], 
           'model.divergence.txt': ['Time', 'theta', 'C_derived'], 
           'model.divergence.2.txt': ['Time', 'theta', 'C_derived'], 
           'model.IM.M_A2D.MM_D2A.txt': ['Time', 'theta', 'Migration', 'C_derived'], 
           'model.IM.M_D2A.MM_D2A.txt': ['Time', 'theta', 'Migration', 'C_derived']
           }



###############################################################################

################################### Functions #################################

class ParameterObj(object):
    def __init__(self, args):
        # initial user parameters
        self.mode = None
        self.model_file = check_file(args.get('--model', None))
        self.model_name = pathlib.Path(args['--model']).name
        self.time = float(args['--time']) if not args['--time'] is None else None
        self.time_low = float(args['--time_low']) if not args['--time_low'] is None else None
        self.time_high = float(args['--time_high']) if not args['--time_high'] is None else None
        self.theta = float(args['--theta']) if not args['--theta'] is None else None
        self.theta_low = float(args['--theta_low']) if not args['--theta_low'] is None else None
        self.theta_high = float(args['--theta_high']) if not args['--theta_high'] is None else None
        self.migration = float(args['--migration']) if not args['--migration'] is None else None
        self.migration_low = float(args['--migration_low']) if not args['--migration_low'] is None else None
        self.migration_high = float(args['--migration_high']) if not args['--migration_high'] is None else None
        self.derived_coalescence = float(args['--derived_Ne']) if not args['--derived_Ne'] is None else None
        self.derived_coalescence_low = float(args['--derived_low']) if not args['--derived_low'] is None else None
        self.derived_coalescence_high = float(args['--derived_high']) if not args['--derived_high'] is None else None
        self.ancestor_population_id = args['--ancestor_population']
        self.kmax = int(args['--kmax']) 
        self.sample_file = check_file(args.get('--sample_file', None))
        self.genome_file = check_file(args.get('--genome_file', None))
        self.windows_file = check_file(args['--windows_hd5']) if not args['--windows_hd5'] is None else None
        self.variants_file = check_file(args['--variants_hd5']) if not args['--variants_hd5'] is None else None
        self.threads = int(args['--threads'])
        self.path = check_path(args.get('--prefix', None))
        self.prefix = check_prefix(args.get('--prefix', None))
        self.dataset = self.prefix if self.path is None else (self.path / self.prefix) 
        self.parameters_numeric = {}
        self.parameters_symbolic = {}
        self.boundaries = {}
        self.check_parameters()
        ##########################################################
        self.max_by_mutype = {mutype: self.kmax for mutype in MUTYPES}
        self.mutuple_space = []

#        self.user_rate = {
#            'Mig' : sympy.Rational(float(args['--migration_rate'])),
#            'C_ancestor' : sympy.Rational(1.0),
#            'C_derived' : sympy.Rational(float(args['--derived_coalescence_rate'])),
#            'Time' : sympy.Rational(float(args['--T_var'])),
#            'theta' : sympy.Rational(float(args['--mutation_rate']))
#        }
#        self.C_ancestor = sympy.Symbol('C_ancestor', positive=True)
#        self.C_derived = sympy.Symbol('C_derived', positive=True)
#        self.M = sympy.Symbol('Mig', positive=True)
#        self.bigL = sympy.Symbol('bigL', positive=True)
#        self.T = sympy.Symbol('Time', positive=True)
#        self.theta_A = sympy.Symbol('theta_A', positive=True)
#        self.theta_B = sympy.Symbol('theta_B', positive=True)
#        self.theta_fix = sympy.Symbol('theta_fix', positive=True)
#        self.theta_AB = sympy.Symbol('theta_AB', positive=True)
#        self.base_rate = {
#            self.M : sympy.Rational(str(float(args['--migration_rate']) / 2)),
#            self.C_ancestor : sympy.Rational(str(1.0)),
#            self.C_derived : sympy.Rational(str(float(args['--derived_coalescence_rate'])))
#        }
#        self.split_time = sympy.Rational(str(float(args['--T_var'])))
#        self.mutation_rate = sympy.Rational(str(float(args['--mutation_rate'])/ 2))
#        self.symbolic_mutation_rate = {
#            'hetA' : self.theta_A,
#            'hetB' : self.theta_B,
#            'hetAB' : self.theta_AB,
#            'fixed' : self.theta_fix
#        }
#
#        self.symbolic_rate = {
#            'C_ancestor' : self.C_ancestor,
#            'C_derived' : self.C_derived, 
#            'Mig' : self.M,
#            'bigL' : self.bigL
#        }
        #self.data = []

    def get_mutuple_counters(self, entityCollection):
        if self.mode == 'variants':
            infile = self.variants_file
        elif self.mode == 'windows':
            infile = self.windows_file
        else:
            pass
        mutype_hdf5_store = pd.HDFStore(infile)
        mutype_df = pd.read_hdf(mutype_hdf5_store, key='mutypes')
        mutype_hdf5_store.close()
        shape = tuple(self.max_by_mutype[mutype] + 2 for mutype in MUTYPES)
        mutuple_count_matrix = np.zeros(shape, np.float64)
        #print(self.ancestor_population_id, (entityCollection.populationObjs[0].id, entityCollection.populationObjs[1].id))
        #print("before")
        #print(mutype_df)
        if self.ancestor_population_id == entityCollection.populationObjs[0].id:
            # mutuples do not have to be flipped
            print("[+] Ancestor is %s ..." % self.ancestor_population_id)
        elif self.ancestor_population_id == entityCollection.populationObjs[1].id:
            mutype_df.rename(columns={'hetA': 'hetB', 'hetB': 'hetA'}, inplace=True)
            print("[+] Ancestor is %s (hetA and hetB will be flipped)... " % self.ancestor_population_id)
        #print("before")
        #print(mutype_df)
        # this has to be changed if order of mutypes changes
        FGV_count = 0
        kmax_binned_count = 0
        total_count = mutype_df['count'].sum()
        for count, hetA, fixed, hetB, hetAB in tqdm(mutype_df[['count'] + MUTYPES].values, total=len(mutype_df.index), desc="[%] ", ncols=100):
            #print(hetA, fixed, hetB, hetAB)
            mutuple = (hetA, fixed, hetB, hetAB)
            if mutuple[1] > 0 and mutuple[3] > 0:
                FGV_count += count  
            if any([count > self.kmax for count in mutuple]):
                kmax_binned_count += count
            mutuple_vector = tuple([count if not count > self.max_by_mutype[mutype] else self.max_by_mutype[mutype] + 1 for count, mutype in zip(mutuple, MUTYPES)])
            
            mutuple_count_matrix[mutuple_vector] += count
            #print(count, hetA, fixed, hetB, hetAB, mutuple_vector, mutuple_count_matrix)
        print("[=] Total mutuple count = %s" % (format_count(total_count)))
        print("[=] Counts excluded due to four-gamete-violations = %s (%s)" % (format_count(FGV_count), format_percentage(FGV_count / total_count)))
        print("[=] Counts binned due to kmax = %s (%s)" % (format_count(kmax_binned_count), format_percentage(kmax_binned_count / total_count)))
        return mutuple_count_matrix

    def check_parameters(self):
        required_parameters = PARAMETERS_BY_MODEL_NAME[self.model_name]
        missing_parameters = []
        if 'Migration' in required_parameters:
            # migration rate per lineage (do the bounds also have to be divided?) 
            self.parameters_symbolic['Migration'] = sympy.Symbol('Migration', positive=True)
            self.parameters_numeric['Migration'] = sympy.Rational(0)
            if not self.migration is None:
                self.parameters_numeric['Migration'] = sympy.Rational(str(self.migration / 2))
            elif (self.migration_low and self.migration_high):
                if self.migration_low < self.migration_high:
                    self.boundaries['Migration'] = (sympy.Rational(str(self.migration_low / 2)), sympy.Rational(str(self.migration_high / 2)))
                else:
                    sys.exit("[X] Lower bound must be smaller than upper bound: %s" % 'Migration')
            else:
                missing_parameters.append('Migration')
                self.parameters_numeric['Migration'] = sympy.Rational(str(0))
        else:
            # set to 0 if no --migration
            self.parameters_symbolic['Migration'] = sympy.Symbol('Migration', positive=True)
            self.parameters_numeric['Migration'] = sympy.Rational(0)
        if 'C_derived' in required_parameters:
            # migration rate per lineage (do the bounds also have to be divided?) 
            self.parameters_symbolic['C_derived'] = sympy.Symbol('C_derived', positive=True)
            if not self.derived_coalescence is None:
                self.parameters_numeric['C_derived'] = sympy.Rational(str(self.derived_coalescence))
            elif (self.derived_coalescence_low and self.derived_coalescence_high):
                if self.derived_coalescence_low < self.derived_coalescence_high:
                    self.boundaries['C_derived'] = (sympy.Rational(str(self.derived_coalescence_low)), sympy.Rational(str(self.derived_coalescence_high)))
                else:
                    sys.exit("[X] Lower bound must be smaller than upper bound: %s" % 'C_derived')
            else:
                missing_parameters.append('C_derived')
        else:
            # set to 0 if no --migration
           pass
        if 'theta' in required_parameters:
            # mutation rate per lineage (do the bounds also have to be divided?) 
            self.parameters_symbolic['theta'] = sympy.Symbol('theta', positive=True)
            self.parameters_symbolic['hetA'] = sympy.Symbol('hetA', positive=True)
            self.parameters_symbolic['fixed'] = sympy.Symbol('fixed', positive=True)
            self.parameters_symbolic['hetB'] = sympy.Symbol('hetB', positive=True)
            self.parameters_symbolic['hetAB'] = sympy.Symbol('hetAB', positive=True)
            if not self.theta is None:
                self.parameters_numeric['theta'] = sympy.Rational(str(self.theta / 2))
            elif (self.theta_low and self.theta_high):
                if self.theta_low < self.theta_high:
                    self.boundaries['theta'] = (sympy.Rational(str(self.theta_low / 2)), sympy.Rational(str(self.theta_high / 2)))
                else:
                    sys.exit("[X] Lower bound must be smaller than upper bound: %s" % 'theta')   
            else:
                missing_parameters.append('theta')
        if 'Time' in required_parameters:
            self.parameters_symbolic['Time'] = sympy.Symbol('Time', positive=True)
            if not self.time is None:
                self.parameters_numeric['Time'] = sympy.Rational(str(self.time))
            elif (self.time_low and self.time_high):
                if self.time_low < self.time_high:
                    self.boundaries['Time'] = (self.time_low, self.time_high)
                else:
                    sys.exit("[X] Lower bound must be smaller than upper bound: %s" % 'Time')      
            else:
                missing_parameters.append('Time')
        if missing_parameters:
            sys.exit("[X] Please specify values or lower/upper bounds for parameter(s): %s" % ",".join(missing_parameters))
        if self.ancestor_population_id is None:
            sys.exit("[X] Please specify a population id using '--ancestor_population'.")
        else:
            self.parameters_numeric['C_ancestor'] = sympy.Rational(1.0)
            self.parameters_symbolic['C_ancestor'] = sympy.Symbol('C_ancestor', positive=True)
        self.parameters_symbolic['BigL'] = sympy.Symbol('BigL', positive=True)
        self.parameters_numeric['BigL'] = sympy.Symbol('BigL', positive=True)
        if self.windows_file:
            self.mode = 'windows'
        if self.variants_file:
            self.mode = 'variants'
        if self.mode is None:
            sys.exit("[X] Please provide input files using '--windows_hd5' or '--variants_hd5'")
        #print("self.parameters_numeric", self.parameters_numeric)
        #print("self.parameters_symbolic", self.parameters_symbolic)

    def generate_mutuple_space(self):
        print("[+] Generating all mutuples for kmax = %s ..." % self.kmax)
        # works only for equal max_mutypes ...
        for mutuple in itertools.product(*[range(0, self.max_by_mutype[mutype] + 1) for mutype in MUTYPES]):
            if not (mutuple[FOUR_GAMETE_VIOLATION_IDX[0]] > 0 and mutuple[FOUR_GAMETE_VIOLATION_IDX[1]] > 0): # FGV
                self.mutuple_space.append(mutuple)
        print("[=] Generated %s mutuples (with kmax = %s)" % (len(self.mutuple_space), self.kmax))

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

def powerset(iterable):
    "powerset([1,2,3]) --> () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)"
    s = list(iterable)
    return itertools.chain.from_iterable(itertools.combinations(s, r) for r in range(len(s)+1))

def create_csv(out_f, header, rows, sep):
    df = pd.DataFrame(rows, columns=header)
    #print(df)
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
    assumptions = "assume(BigL >= 0)"
    #print()
    '''
    [ratnormal] rewrites an expression using its irreducible representation. The expression is viewed as a multivariate rational fraction with coefficients in Q (or
        Q[i]). The variables are generalized identifiers which are assumed to be algebraically independent. Unlike with normal, an algebraic extension is considered
        as a generalized identifier. Therefore ratnormal is faster but might miss some
        simplifications if the expression contains radicals or algebraically dependent transcendental functions.
    '''
    invlaplace_string = "%s; invlaplace(ratnormal(subst(%s%s) / BigL), BigL, T).subst(T, %s)" % (assumptions, equation_giac, mutation_substitution, float(split_time))
    try:
        process = subprocess.run(["giac", invlaplace_string], stderr=subprocess.PIPE, stdout=subprocess.PIPE, check=True, encoding='utf-8')
        probability = process.stdout.rstrip("\n").split(",")[1]
    except subprocess.CalledProcessError:
        exit("[X] giac could not run.")
    #if probability == 'undef':
    #    print(invlaplace_string)
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
            self.path_numerators.append((event_count * parameterObj.parameters_symbolic[event_id]))
            self.path_denominators.append(
                        (nodeObj.event_counter['C_ancestor'] * parameterObj.parameters_symbolic['C_ancestor']) + 
                        (nodeObj.event_counter['C_derived'] * parameterObj.parameters_symbolic['C_derived']) + 
                        (nodeObj.event_counter['Migration'] * parameterObj.parameters_symbolic['Migration']) + 
                        (nodeObj.event_counter['BigL'] * parameterObj.parameters_symbolic['BigL']) +
                        (nodeObj.mutype_counter['hetA'] * parameterObj.parameters_symbolic['hetA']) +
                        (nodeObj.mutype_counter['hetB'] * parameterObj.parameters_symbolic['hetB']) +
                        (nodeObj.mutype_counter['hetAB'] * parameterObj.parameters_symbolic['hetAB']) +
                        (nodeObj.mutype_counter['fixed'] * parameterObj.parameters_symbolic['fixed'])
                )
        for numerator, denominator in zip(self.path_numerators, self.path_denominators):
            self.path_equation *= numerator / denominator

    def multiply_with_path_equation(self, mutation_equation):
        initial = 1
        for numerator, denominator in zip(self.path_numerators, self.path_denominators):
            initial *= numerator / denominator
        mutation_equation *= initial
        return mutation_equation

    def infer_mutation_equation(self, mutuple, parameterObj):
        mutation_equation = 0
        for slots in itertools.product(*[list(itertools.combinations_with_replacement(self.mutable_nodes_by_mutype[mutype], count)) for count, mutype in zip(mutuple, MUTYPES)]):
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
                    mutation_part *= (nodeObj.mutype_counter[mutype] * parameterObj.parameters_symbolic[mutype]) / denominator
            mutation_equation += mutation_part
        return mutation_equation

    def infer_mutable_nodes(self):
        for idx, nodeObj in enumerate(self.nodeObjs):
            for mutype, count in nodeObj.mutype_counter.items():
                if count:
                    self.mutable_nodes_by_mutype[mutype].append(idx)

    def is_compatible_with_data_point(self, mutuple):
        if sum(mutuple) == 0:
            return True
        else:
            mutypes_in_mutuple = set([mutype for mutype, count in zip(MUTYPES, mutuple) if count > 0])
            if all([(True if mutype in self.mutable_nodes_by_mutype else False) for mutype in mutypes_in_mutuple]):
                return True
        return False

    def yield_step(self):
        for nodeObj, event_id, event_count in zip(self.nodeObjs, self.event_ids, self.event_counts):
            yield (nodeObj, event_id, event_count)

    def __str__(self):
        return "# Path %s:\n%s" % (self.path_id, 
            "\n".join(["%s: %s %s" % (nodeObj.node_id, event_id, str(nodeObj.mutype_counter)) for nodeObj, event_id in zip(self.nodeObjs, self.event_ids)])
            )

def infer_composite_likelihood(x0, *args):
    simplex_parameters, symbolic_equations_by_mutuple, mutuple_count_matrix, parameterObj = args
    numeric_value_by_symbol_for_substitution = {}
    numeric_value_by_parameter = {}
    theta_by_symbol = {}
    x0_by_parameter = {}
    #print(x0[0], type(x0[0]))
    if simplex_parameters:
        # round here? if float too long then causes undef in ILT
        x0_by_parameter = {parameter: sympy.Rational(str(round(x, 2))) for parameter, x in zip(simplex_parameters, list(x0))}
        #print(x0_by_parameter)
    for parameter, symbol in parameterObj.parameters_symbolic.items():
        if parameter in x0_by_parameter:
            numeric_value_by_symbol_for_substitution[symbol] = x0_by_parameter[parameter]
            numeric_value_by_parameter[parameter] = sympy.Rational(str(x0_by_parameter[parameter]))
            if parameter == 'theta':
                theta_by_symbol[symbol] = sympy.Rational(str(x0_by_parameter[parameter]))
        else:
            if parameter in set(MUTYPES):
                theta_by_symbol[symbol] = x0_by_parameter['theta'] if 'theta' in x0_by_parameter else parameterObj.parameters_numeric['theta']
            else:
                numeric_value_by_symbol_for_substitution[symbol] = parameterObj.parameters_numeric[parameter]
                numeric_value_by_parameter[parameter] = parameterObj.parameters_numeric[parameter]             
    #print(numeric_value_by_parameter)
    #print(numeric_value_by_symbol_for_substitution)      
    equation_by_mutuple = {
        mutuple: sum(equation).xreplace(numeric_value_by_symbol_for_substitution)
        for mutuple, equation in symbolic_equations_by_mutuple.items()
        }
    # Setting up datastructure
    shape = tuple(parameterObj.max_by_mutype[mutype] + 2 for mutype in MUTYPES)
    probability_matrix = np.zeros(shape, np.float64)
    # Generating sets of mutation_rate combinations for marginals
    vector_seen = set([])
    for mutypes in reversed(list(powerset(MUTYPES))):
        mutypes_masked = set([mutype for mutype in MUTYPES if not mutype in mutypes])
        if not FOUR_GAMETE_VIOLATION.issubset(mutypes_masked):
            equations_by_vector = collections.defaultdict(list)
            marginal_queries = []
            vectors = []
            rates_by_mutype_by_vector = {}
            for data_tuple in sorted(equation_by_mutuple.keys()):               
                vector = tuple([count if not mutype in mutypes_masked else parameterObj.max_by_mutype[mutype] + 1 for count, mutype in zip(data_tuple, MUTYPES)])
                if not vector in vector_seen:
                    vector_seen.add(vector)
                    if not len([mutype for mutype, count in zip(MUTYPES, vector) if mutype in FOUR_GAMETE_VIOLATION and count > 0]) > 1: 
                        vectors.append(vector)
                        rates_by_mutype_by_vector[vector] = {parameterObj.parameters_symbolic[mutype]: (numeric_value_by_parameter['theta'] if mutype in mutypes else 0) for mutype in MUTYPES}
                        if not mutypes_masked:
                            marginal_query = None
                        else:
                            marginal_query = tuple(
                                [(data_tuple[idx] if not mutype in mutypes_masked else slice(0, parameterObj.max_by_mutype[mutype] + 2)) 
                                    for idx, mutype in enumerate(MUTYPES)])
                        marginal_queries.append(marginal_query)
                        equations_by_vector[vector] = equation_by_mutuple[data_tuple]
            # can be parallelised 
            params = [(vector, marginal_query, equations_by_vector[vector], rates_by_mutype_by_vector[vector], numeric_value_by_parameter['Time']) for vector, marginal_query in zip(vectors, marginal_queries)]
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
    if simplex_parameters is None:
        composite_likelihood = -np.sum((xlogy(np.sign(probability_matrix), probability_matrix) * mutuple_count_matrix))
        print('[+] L=-%s\t%s' % (composite_likelihood, ", ".join(["%s=%s" % (param, round(value, 4)) for param, value in numeric_value_by_parameter.items() if not param == 'BigL'])))
        return composite_likelihood
    composite_likelihood = -np.sum((xlogy(np.sign(probability_matrix), probability_matrix) * mutuple_count_matrix))
    print(" " * 100, end='\r'),
    print("[O] L=-%s\t%s" % (composite_likelihood, ", ".join(["%s=%s" % (param, round(value, 4)) for param, value in numeric_value_by_parameter.items() if not param == 'BigL'])))
    return composite_likelihood

def calculate_likelihood(symbolic_equations_by_mutuple, mutuple_count_matrix, parameterObj):
    start_time = timer()
    print("[+] Calculating composite Likelihood ...")
    x0, simplex_parameters = None, None
    composite_likelihood = infer_composite_likelihood(x0, simplex_parameters, symbolic_equations_by_mutuple, mutuple_count_matrix, parameterObj)
    print("[=] Calculated composite Likelihood (L=%s) in %s seconds..." % (composite_likelihood, timer() - start_time))

def prepare_paths(parameterObj):
    # infer path_probabilities for each pathObj (only has to be done once) ... 
    start_time = timer()
    pathObj_by_path_id = read_paths(parameterObj)
    print("[+] Preparing %s paths (path equations and mutable nodes) ..." % (len(pathObj_by_path_id)))
    for path_id, pathObj in tqdm(pathObj_by_path_id.items(), total=len(pathObj_by_path_id), desc="[%] ", ncols=100):
        pathObj.infer_mutable_nodes()
        pathObj.infer_path_quation(parameterObj)
    print("[=] Prepared paths in %s seconds." % (timer() - start_time))
    return pathObj_by_path_id

def read_paths(parameterObj):
    start_time = timer()
    print("[+] Reading paths in file %s..." % parameterObj.model_file)
    paths_df = pd.read_csv(parameterObj.model_file, sep="\t")
    pathObj_by_path_id = {}
    nodeObj_by_node_id = {}
    for path_id, node_id, event_id, event_count, C_ancestor, C_derived, Migration, BigL, hetA, fixed, hetB, hetAB in tqdm(paths_df.values.tolist(), total=len(paths_df.index), desc="[%] ", ncols=100):
        if not event_id == 'LCA':
            if not node_id in nodeObj_by_node_id:
                mutype_counter = collections.Counter({'hetA': hetA, 'fixed': fixed, 'hetB': hetB, 'hetAB': hetAB})
                event_counter = collections.Counter({'C_ancestor': C_ancestor, 'C_derived': C_derived, 'Migration': Migration, 'BigL': BigL})
                nodeObj_by_node_id[node_id] = NodeObj(node_id, mutype_counter, event_counter)
            if not path_id in pathObj_by_path_id:
                pathObj_by_path_id[path_id] = PathObj(path_id)
            pathObj_by_path_id[path_id].add_step(nodeObj_by_node_id[node_id], event_id, event_count)
    print("[=] Parsed %s paths from file %s in %s seconds." % (len(pathObj_by_path_id), parameterObj.model_file, timer() - start_time))
    return pathObj_by_path_id         

def generate_equations(pathObj_by_path_id, parameterObj):
    start_time = timer()
    args = []
    equations_by_data_tuple = collections.defaultdict(list)
    for mutuple in parameterObj.mutuple_space:
        for path_id, pathObj in pathObj_by_path_id.items():
            if pathObj.is_compatible_with_data_point(mutuple):
                args.append((  
                    mutuple, 
                    pathObj, 
                    parameterObj))
            else:
                #equation is 0 if data point can't be placed on path ...
                equations_by_data_tuple[mutuple].append(0)                            
    print("[+] Analysing %s combinations of paths and data points with %s threads..." % (len(args), parameterObj.threads))
    if parameterObj.threads == 1:
        for param in tqdm(args, total=len(args), desc="[%] ", ncols=100):
            data_tuple, equation = build_equation(param)
            equations_by_data_tuple[data_tuple].append(equation)
    else:
        with poolcontext(processes=parameterObj.threads) as pool:
            with tqdm(args, total=len(args), desc="[%] ", ncols=100) as pbar:
                for data_tuple, equation in pool.imap_unordered(build_equation, args):
                    equations_by_data_tuple[data_tuple].append(equation)
                    pbar.update()
    print("[=] Analysed paths in %s seconds." % (timer() - start_time))
    return equations_by_data_tuple

def build_equation(args):
    mutuple, pathObj, parameterObj = args
    mutation_equation = pathObj.infer_mutation_equation(mutuple, parameterObj)
    equation = pathObj.multiply_with_path_equation(mutation_equation)
    return mutuple, equation

def generate_initial_simplex(boundaries, seed):
    np.random.seed(seed)
    simplex = [] 
    simplex_order = list(boundaries.keys())
    for i in range(len(boundaries) + 1):
        vertex = []
        for parameter, (minval, maxval) in boundaries.items():
            value = sympy.Rational(str(round(np.random.uniform(minval, maxval, 1)[0], 2)))
            vertex.append(value)
        simplex.append(tuple(vertex))
    return simplex, simplex_order

def task_generate_entityCollection(parameterObj):
    start = timer()
    print("[#] Building entities based on samples and sequences...")
    entityCollection = EntityCollection()
    entityCollection.parse_sample_file(parameterObj)
    print("[+] Read %s samples from %s populations and generated %s pairs in %.3fs." % (\
        entityCollection.count('samples'), \
        entityCollection.count('populations'), \
        entityCollection.count('pairs'), \
        timer() - start))
    entityCollection.parse_genome_file(parameterObj)
    print("[+] Read %s sequences with total length of %s b in %.3fs" % (\
        entityCollection.count('sequences'), \
        entityCollection.count('bases'), \
        timer() - start))
    return entityCollection

def estimate_parameters(symbolic_equations_by_mutuple, mutuple_count_matrix, parameterObj, seed):
    print("[+] Optimising parameters: %s ..." % (", ".join(parameterObj.boundaries.keys())))    
    start_time = timer()
    simplex_values, simplex_parameters = generate_initial_simplex(parameterObj.boundaries, seed)
    x0 = tuple([0] * len(parameterObj.boundaries.keys()))
    #block_count = mutuple_count_matrix.flatten().sum()
    res = minimize(
        infer_composite_likelihood, 
        x0, 
        args=(simplex_parameters, symbolic_equations_by_mutuple, mutuple_count_matrix, parameterObj), 
        method="Nelder-Mead", 
        options={
            'initial_simplex': simplex_values, 
            'maxfev' : 200,
            'maxiter': 200,
            'disp': False, 
            'xatol': 1e-3, 
            'fatol': 1e-3, # * block_count, # needs to be scaled by number of blocks
            'adaptive': True})
    print()
    if res.success:
        estimated_parameters = collections.OrderedDict({key: value for (key, _), value in zip(parameterObj.boundaries.items(), res.x)})
        estimated_parameters_string = ", ".join(["%s=%s" % (key, round(value, 4)) for key, value in estimated_parameters.items()])
        print("[+] Parameters estimated in %ss using %s iterations (Composite Likelihood = -%s): %s" % (timer() - start_time, res.nit, res.fun, estimated_parameters_string))
    else:
        print("[-] No covergence reached after %s iterations (%ss elapsed)" % (res.nit, timer() - start_time))
    return estimated_parameters