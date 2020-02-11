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
import pprint
pp = pprint.PrettyPrinter(indent=4)



'''
variable definition
equation building
    sage
    giac:   
        - write to file


'''

################################### CONSTANTS #################################

MUTYPES = ['hetA', 'fixed', 'hetB', 'hetAB']
MUTYPE_SET = set(MUTYPES)
EVENTS = ['C_ancestor', 'C_derived', 'Migration', 'BigL']
EVENT_SPACE = EVENTS + MUTYPES

FOUR_GAMETE_VIOLATION_IDX = [idx for idx, mutype in enumerate(MUTYPES) if mutype in set(['fixed', 'hetAB'])]
FOUR_GAMETE_VIOLATION = set(['fixed', 'hetAB'])

PARAMETERS_BY_MODEL_NAME = {
           #'model.M.txt': ['theta', 'Mig'], 
           'model.divergence.txt': ['Time', 'theta', 'C_derived'], 
           'model.divergence.1.txt': ['Time', 'theta', 'C_derived'], 
           'model.IM.M_A2D.MM_D2A.txt': ['Time', 'theta', 'Migration', 'C_derived'], 
           'model.IM.M_D2A.MM_D2A.txt': ['Time', 'theta', 'Migration', 'C_derived']
           }



###############################################################################

################################### Functions #################################

'''

1. Parse model
    - based on model check requirements from header
        [0:4]: path_idx, origin, idx, label
        - Coalescence rates: one for each C=* in header
        - Migration rates : one for each M=* in header
        - Join rates : one for each * in header

    => read paths 
Goals
1. Model to equations: needs
    - paths
    - values
2. equations to probabilities
    - data 
'''


class ParameterObj(object):
    def __init__(self, args):
        self.model_file = check_file(args.get('--model', None))
        self.zarr_store = None

        self.time = float(args['--time']) if not args['--time'] is None else None
        self.theta = float(args['--theta']) if not args['--theta'] is None else None
        self.migration = float(args['--migration']) if not args['--migration'] is None else None
        '''

        self.pop_ids = 'A,B,A&B'
        self.pop_Ns = '1.2,1.2,1.4'
        slef
        '''

        self.derived_coalescence = float(args['--derived_Ne']) if not args['--derived_Ne'] is None else None


        self.kmax = int(args['--kmax'])
        self.seed = int(args['--seed']) 

        self.sample_file = check_file(args.get('--sample_file', None))
        self.genome_file = check_file(args.get('--genome_file', None))
        #self.windows_file = check_file(args['--windows_hd5']) if not args['--windows_hd5'] is None else None
        self.windows_file = None
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
            self.parameters_numeric['hetA'] = sympy.Rational(str(self.theta / 2)) # test
            self.parameters_numeric['fixed'] = sympy.Rational(str(self.theta / 2)) # test
            self.parameters_numeric['hetB'] = sympy.Rational(str(self.theta / 2)) # test
            self.parameters_numeric['hetAB'] = sympy.Rational(str(self.theta / 2)) # test
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
        print("[=] Generated %s mutuples" % (len(self.mutuple_space)))

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

def giacify_equations(equations, rates=None, split_time=None):
    '''
    return equation as GIAC-compatible string

    if rates and split_time 
        -> symbols are substituted by values
    else:
        -> symbolic equation
    '''
    giac_equations = []
    for equation in equations:
        giac_string = []        
        giac_string.append("normal(invlaplace(normal(")
        if not rates is None and split_time is None:
            giac_string.append(str(equation.xreplace(rates)).replace("**", "^"))
            giac_string.append(") / BigL, BigL, Time)).subst(Time, %s)" % str(float(split_time)))
        else:
            giac_string.append(str(equation).replace("**", "^"))
            giac_string.append(") / BigL, BigL, Time))")
        giac_equations.append("".join(giac_string))
    return " + ".join(giac_equations)

def inverse_laplace_transform(params):
    # with tempfile.NamedTemporaryFile(mode='w') as temp_fh:
    #     temp_fh.write("%s\n" % "\n".join(giac_strings))
    #     try:
    #         process = subprocess.run(["giac", temp_fh.name], stderr=subprocess.PIPE, stdout=subprocess.PIPE, check=True, encoding='utf-8')
    #     except subprocess.CalledProcessError:
    #         exit("[X] giac could not run.")
    #     probabilities = [np.float64(sympy.Rational(str(probability))) for probability in process.stdout.split(",\n")]
    vector, giac_string = params
    try:
        process = subprocess.run(["giac", giac_string], stderr=subprocess.PIPE, stdout=subprocess.PIPE, check=True, encoding='utf-8')
        results = str(process.stdout)
    except subprocess.CalledProcessError:
        exit("[X] giac could not run.")
    if results == 'undef':
        print(giac_string, results)
    return vector, sum([np.float64(sympy.Rational(result)) for result in results.split(",")])

@contextmanager
def poolcontext(*args, **kwargs):
    pool = Pool(*args, **kwargs)
    yield pool
    pool.terminate()

###############################################################################

################################### Classes ###################################

def read_model(parameterObj):
    start_time = timer()
    print("[+] Reading paths in file %s..." % parameterObj.model_file)
    paths_df = pd.read_csv(parameterObj.model_file, sep="\t")
    pathObj_by_path_id = {}
    for path_id, node_id, event, count, C_ancestor, C_derived, Migration, BigL, hetA, fixed, hetB, hetAB in tqdm(
            paths_df.values.tolist(), total=len(paths_df.index), desc="[%] ", ncols=100):
        if not event == 'LCA':
            event_space = collections.Counter({
                event_type: event_count for event_type, event_count in zip(
                    EVENT_SPACE, [C_ancestor, C_derived, Migration, BigL, hetA, fixed, hetB, hetAB]) 
                        if event_count
                })
            event_count = collections.Counter({event: count})
            if not path_id in pathObj_by_path_id:
                pathObj_by_path_id[path_id] = PathObj(path_id)
            pathObj_by_path_id[path_id].add_event(event_count, event_space)
    print("[=] Parsed %s paths from file %s in %s seconds." % (len(pathObj_by_path_id), parameterObj.model_file, timer() - start_time))
    #for path_id, pathObj in pathObj_by_path_id.items():
    #    pp.pprint(pathObj.__dict__)
    return pathObj_by_path_id 

def generate_event_equations(parameterObj, pathObj_by_path_id):
    event_equations_by_mutuple = collections.defaultdict(list)
    for mutuple in parameterObj.mutuple_space:
        print("## mutuple", mutuple)
        for path_id, pathObj in pathObj_by_path_id.items():
            if pathObj.is_compatible_with_data_point(mutuple):
                eq_numerators, eq_denominators = [], []
                # path equation
                for event_counter, event_space in zip(pathObj.event_counters, pathObj.event_spaces):
                    event_type = list(event_counter.keys())[0]
                    event_count = list(event_counter.values())[0]
                    eq_numerators.append(event_count * parameterObj.parameters_symbolic[event_type])
                    denominator = []
                    for event_space_type, event_space_count in event_space.items():
                        denominator.append(event_space_count * parameterObj.parameters_symbolic[event_space_type])
                    eq_denominators.append(sum(denominator))
                event_equation = 1
                for eq_numerator, eq_denominator in zip(eq_numerators, eq_denominators):
                    event_equation *= eq_numerator / eq_denominator
                mutation_equation = 0
                for slots in itertools.product(*[           
                    list(itertools.combinations_with_replacement(pathObj.slots_by_mutype[mutype], count)) 
                           for count, mutype in zip(mutuple, MUTYPES)]):
                    mutypes_by_idx = collections.defaultdict(list)
                    for _mutype, slot in zip(MUTYPES, slots):
                        if slot: 
                            for idx in list(slot):
                                mutypes_by_idx[idx].append(_mutype) 
                    mutation_part = 1
                    for idx, mutypes in mutypes_by_idx.items():
                        mutype_counter = collections.Counter(mutypes)
                        mutation_part *= multinomial(mutype_counter.values())
                        
                        for mutype in mutypes:
                            mutation_part *= (pathObj.event_spaces[idx][mutype] * parameterObj.parameters_symbolic[mutype]) / eq_denominators[idx]
                    mutation_equation += mutation_part
                mutation_equation *= event_equation 
                event_equations_by_mutuple[mutuple].append(mutation_equation)
                
    return event_equations_by_mutuple

def infer_composite_likelihood(event_equations_by_mutuple, raw_rates, mutuple_count_matrix, parameterObj):
    # Setting up datastructure
    shape = tuple(parameterObj.max_by_mutype[mutype] + 2 for mutype in MUTYPES)
    probability_matrix = np.zeros(shape, np.float64)
    marginal_query_by_vector = {}
    params = []
    vector_space = sorted(list(itertools.product(range(parameterObj.kmax+2), repeat=len(MUTYPES))))
    for vector in vector_space:
        # check if not FGV
        if not all([vector[idx] > 0 for idx in FOUR_GAMETE_VIOLATION_IDX]): 
            # set mutype rates, marginals ...
            rates = {k: v for k, v in raw_rates.items()}
            if not vector in event_equations_by_mutuple:
                marginal_query = []
                equation_mutype = []
                for idx, mutype in enumerate(MUTYPES):
                    if vector[idx] > parameterObj.kmax:
                        symbol = parameterObj.parameters_symbolic[mutype]
                        rates[symbol] = 0
                        marginal_query.append(slice(0, parameterObj.kmax + 2))
                        equation_mutype.append(0)
                    else:
                        marginal_query.append(vector[idx])
                        equation_mutype.append(vector[idx])
                marginal_query_by_vector[vector] = tuple(marginal_query)
                equations = event_equations_by_mutuple[tuple(equation_mutype)]
            else:
                marginal_query_by_vector[vector] = None
                equations = event_equations_by_mutuple[vector]        
            giac_string = giacify_equations(equations, rates, parameterObj.parameters_numeric['Time'])
            params.append([vector, giac_string])
    probability_by_vector = {}
    if parameterObj.threads < 2:
        for param in tqdm(params, total=len(params), desc="[%] ", ncols=100):
            vector, probability = inverse_laplace_transform(param)
            probability_by_vector[vector] = probability
    else:
        with poolcontext(processes=parameterObj.threads) as pool:
            with tqdm(params, total=len(params), desc="[%] ", ncols=100) as pbar:
                for vector, probability in pool.imap_unordered(inverse_laplace_transform, params):
                    probability_by_vector[vector] = probability
                    pbar.update()
    for vector, probability in probability_by_vector.items():
        if not marginal_query_by_vector[vector] is None:
            probability_matrix[vector] = probability - sum(probability_matrix[marginal_query_by_vector[vector]].flatten())
        else:
            probability_matrix[vector] = probability
    for vector, _ in params:
        print(vector, "\t", probability_matrix[vector])
    print(len(params), sum(probability_matrix.flatten()))

    print_params = { (param):(value if not param == 'theta' else value * 2) for param, value in parameterObj.numeric_value_by_parameter.items()}
#    #if simplex_parameters is None:
    #    composite_likelihood = -np.sum((xlogy(np.sign(probability_matrix), probability_matrix) * mutuple_count_matrix))
    #    print('[+] L=-%s\t%s' % (composite_likelihood, ", ".join(["%s=%s" % (param, round(value, 4)) for param, value in print_params.items() if not param == 'BigL'])))
    #    return composite_likelihood
    composite_likelihood = -np.sum((xlogy(np.sign(probability_matrix), probability_matrix) * mutuple_count_matrix))
    print(" " * 100, end='\r'),
    print("[O] L=-%s\t%s" % (composite_likelihood, ", ".join(["%s=%s" % (param, round(value, 4)) for param, value in print_params.items() if not param == 'BigL'])))
    return composite_likelihood

def calculate_likelihood(event_equations_by_mutuple, mutuple_count_matrix, parameterObj):
    start_time = timer()
    print("[+] Calculating composite Likelihood ...")
    rates = {
        parameterObj.parameters_symbolic['hetA']: parameterObj.parameters_numeric['hetA'],
        parameterObj.parameters_symbolic['fixed']: parameterObj.parameters_numeric['fixed'],
        parameterObj.parameters_symbolic['hetB']: parameterObj.parameters_numeric['hetB'],
        parameterObj.parameters_symbolic['hetAB']: parameterObj.parameters_numeric['hetAB'],
        parameterObj.parameters_symbolic['C_ancestor']: parameterObj.parameters_numeric['C_ancestor'],
        parameterObj.parameters_symbolic['C_derived']: parameterObj.parameters_numeric['C_derived'],
        parameterObj.parameters_symbolic['Migration']: parameterObj.parameters_numeric['Migration']
        }
    composite_likelihood = infer_composite_likelihood(event_equations_by_mutuple, rates, mutuple_count_matrix, parameterObj)
    print("[=] Calculated composite Likelihood (L=%s) in %s seconds..." % (composite_likelihood, timer() - start_time))

class PathObj(object):
    def __init__(self, path_id):
        self.path_id = path_id
        self.steps = 0
        self.event_counters = []    # list of event:count (what happens)
        self.event_spaces = []      # list of events:counts (all things that can happen)
        self.slots_by_mutype = collections.defaultdict(list)

        self.path_numerators = []
        self.path_denominators = []
        self.path_equation = 1

    def add_event(self, event_count, event_space):
        self.event_counters.append(event_count)
        self.event_spaces.append(event_space)
        for mutype in MUTYPE_SET.intersection(event_space):
            self.slots_by_mutype[mutype].append(self.steps)
        self.steps += 1

    def is_compatible_with_data_point(self, mutuple):
        if sum(mutuple) == 0:
            return True
        else:
            mutypes_in_mutuple = set([mutype for mutype, count in zip(MUTYPES, mutuple) if count > 0])
            if all([(True if mutype in self.slots_by_mutype else False) for mutype in mutypes_in_mutuple]):
                return True
        return False

    def yield_step(self):
        for nodeObj, event_id, event_count in zip(self.nodeObjs, self.event_ids, self.event_counts):
            yield (nodeObj, event_id, event_count)

    def __str__(self):
        return "# Path %s:\n%s" % (self.path_id, 
            "\n".join(["%s: %s %s" % (nodeObj.node_id, event_id, str(nodeObj.mutype_counter)) for nodeObj, event_id in zip(self.nodeObjs, self.event_ids)])
            )    

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

def estimate_parameters(symbolic_equations_by_mutuple, mutuple_count_matrix, parameterObj):
    print("[+] Optimising parameters: %s ..." % (", ".join(parameterObj.boundaries.keys())))    
    start_time = timer()
    simplex_values, simplex_parameters = generate_initial_simplex(parameterObj.boundaries, parameterObj.seed)
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
            'xatol': 1e-1, 
            'fatol': 1e-1, # * block_count, # needs to be scaled by number of blocks
            'adaptive': True})
    print()
    if res.success:
        estimated_parameters = collections.OrderedDict({key: value for (key, _), value in zip(parameterObj.boundaries.items(), res.x)})
        print_params = { (param):(value if not param == 'theta' else value * 2) for param, value in estimated_parameters.items()}
        estimated_parameters_string = ", ".join(["%s=%s" % (key, round(value, 4)) for key, value in print_params.items()])
        print("[+] Parameters estimated in %ss using %s iterations (Composite Likelihood = -%s): %s" % (timer() - start_time, res.nit, res.fun, estimated_parameters_string))
    else:
        print("[-] No covergence reached after %s iterations (%ss elapsed)" % (res.nit, timer() - start_time))
    return estimated_parameters