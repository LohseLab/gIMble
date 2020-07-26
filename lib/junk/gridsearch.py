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
from lib.functions import format_percentage, create_hdf5_store
from lib.functions import check_path, check_file, check_prefix, format_count, plot_parameter_scan
import pathlib
import tempfile
import more_itertools

################################### CONSTANTS #################################

MUTYPES = ['hetA', 'fixed', 'hetB', 'hetAB']
FOUR_GAMETE_VIOLATION_IDX = [idx for idx, mutype in enumerate(MUTYPES) if mutype in set(['fixed', 'hetAB'])]
FOUR_GAMETE_VIOLATION = set(['fixed', 'hetAB'])

PARAMETERS_BY_MODEL_NAME = {
           #'model.M.txt': ['theta', 'Mig'], 
           'model.divergence.txt': ['Time', 'theta', 'C_derived'], 
           'model.IM.M_A2D.MM_D2A.txt': ['Time', 'theta', 'Migration', 'C_derived'], 
           'model.IM.M_D2A.MM_D2A.txt': ['Time', 'theta', 'Migration', 'C_derived']
           }

GRID_PARAMETER_ORDER = ['C_ancestor', 'C_derived', 'Migration', 'theta', 'Time']

###############################################################################

################################### Functions #################################

class ParameterObj(object):
    def __init__(self, args):
        # initial user parameters
        self.mode = None
        self.model_file = check_file(args.get('--model', None))
        self.model_name = pathlib.Path(args['--model']).name
        self.time_MLE = float(args['--time_MLE']) if not args['--time_MLE'] is None else None
        #self.migration_MLE = float(args['--migration_MLE']) if not args['--migration_MLE'] is None else None
        self.mutation_rate = float(args['--mu']) if not args['--mu'] is None else None
        self.derived_coalescence_MLE = float(args['--derived_MLE']) if not args['--derived_MLE'] is None else None

        self.block_size = float(args['--block_size']) if not args['--block_size'] is None else None
        self.theta_low = float(args['--theta_low']) if not args['--theta_low'] is None else None
        self.theta_high = float(args['--theta_high']) if not args['--theta_high'] is None else None
        self.migration_high = float(args['--migration_high']) if not args['--migration_high'] is None else None
        self.migration_low = float(args['--migration_low']) if not args['--migration_low'] is None else None

        self.grid_raw_file = check_file(args.get('--grid', None))

        self.ancestor_population_id = args['--ancestor_population']
        self.kmax = int(args['--kmax']) 
        self.sample_file = check_file(args.get('--sample_file', None))
        self.genome_file = check_file(args.get('--genome_file', None))
        self.windows_file = check_file(args['--windows_hd5']) if not args['--windows_hd5'] is None else None
        self.threads = int(args['--threads'])
        self.path = check_path(args.get('--prefix', None))
        self.prefix = check_prefix(args.get('--prefix', None))
        self.dataset = self.prefix if self.path is None else (self.path / self.prefix) 
        self.parameters_numeric = {}
        self.parameters_symbolic = {}
        self.check_parameters()
        self.grid_raw = None
        self.grid_gimble = None
        ##########################################################
        self.max_by_mutype = {mutype: self.kmax for mutype in MUTYPES}
        self.mutuple_space = []
        self.probability_matrix_by_parameter_tuple = {}
        self.window_pos_by_window_id = {}
        self.x_boundaries = []

    def parse_grid_file(self, entityCollection):
        print("[+] Parsing grid parameter file: %s ..." % self.grid_raw_file)
        grid_df = pd.read_csv(
            self.grid_raw_file, 
            sep="\t", 
            names=['Migration', 'theta_ancestor', 'theta_derived', 'hetA', 'hetB', 'fixed', 'hetAB','probability'], 
            skiprows=1, 
            header=None)
        shape = tuple(self.max_by_mutype[mutype] + 2 for mutype in MUTYPES)        
        # DO THEY HAVE TO BE FLIPPED?
        #grid_df.rename(columns={'hetA': 'hetB', 'hetB': 'hetA'}, inplace=True)
        for migration, theta_ancestor, theta_derived, hetA, hetB, fixed, hetAB, probability in tqdm(grid_df[['Migration', 'theta_ancestor', 'theta_derived', 'hetA', 'hetB', 'fixed', 'hetAB','probability']].values.tolist(), total=len(grid_df.index), desc="[%] ", ncols=100):
            #print(migration, theta_ancestor, theta_derived, hetA, hetB, fixed, hetAB, probability)
            parameter_tuple = (theta_ancestor, theta_derived, migration)
            if not parameter_tuple in self.probability_matrix_by_parameter_tuple:
                self.probability_matrix_by_parameter_tuple[parameter_tuple] = np.zeros(shape, np.float64)
            # this could still be source of error
            vector = tuple([int(hetA), int(fixed), int(hetB), int(hetAB)])
            self.probability_matrix_by_parameter_tuple[parameter_tuple][vector] = probability

    def parse_window_positions(self, entityCollection):
        window_hdf5_store = pd.HDFStore(self.windows_file)
        window_df = pd.read_hdf(window_hdf5_store, key='window_metrics')
        window_hdf5_store.close()
        
        position_by_window_id = {}
        offset_by_sequence_id = {}
        offset = 0
        x_boundaries = []
        for sequenceObj in entityCollection.sequenceObjs:
            offset_by_sequence_id[sequenceObj.id] = offset
            x_boundaries.append(offset)
            offset += sequenceObj.length
        x_boundaries.append(offset)
        window_df['rel_pos'] = window_df['centre'] + window_df['sequence_id'].map(offset_by_sequence_id)
        for window_id, rel_pos in window_df[['window_id', 'rel_pos']].values:
            self.window_pos_by_window_id[window_id] = rel_pos
        self.x_boundaries = x_boundaries

    def setup_grid(self):
        theta_step = (self.theta_high - self.theta_low) / 8
        migration_step = (self.migration_high - self.migration_low) / 10
        grid_raw = []
        grid_gimble = []
        # make migration_low/high
        # c_derived as float
        test_limit = 4
        i = 0
        for migration in np.arange(
                self.migration_low, 
                self.migration_high + migration_step, 
                migration_step):
            for theta_ancestral in np.arange(self.theta_low, self.theta_high, theta_step):
                for theta_derived in np.arange(
                        self.theta_low / self.derived_coalescence_MLE, 
                        self.theta_high / self.derived_coalescence_MLE, 
                        theta_step / self.derived_coalescence_MLE):
                    if i < test_limit:
                        grid_raw.append((theta_ancestral, theta_derived, migration)) 
                        i += 1
        for theta_ancestral, theta_derived, migration in grid_raw:
            theta = theta_ancestral / 2
            C_ancestor = 1
            C_derived = (theta_ancestral / theta_derived)
            Migration = ((migration * theta_ancestral) / (self.block_size * self.mutation_rate) / 2) # per branch
            Time = (self.time_MLE * 2 * self.block_size * self.mutation_rate) / theta_ancestral
            grid_gimble.append((C_ancestor, C_derived, Migration, theta, Time))
        #print(grid_raw)
        #print(grid_gimble)
        self.grid_raw = grid_raw
        self.grid_gimble = grid_gimble

    def get_mutuple_count_matrix_by_window_id(self, entityCollection):
        mutype_hdf5_store = pd.HDFStore(self.windows_file)
        mutype_df = pd.read_hdf(mutype_hdf5_store, key='mutypes')
        mutype_hdf5_store.close()
        if self.ancestor_population_id == entityCollection.populationObjs[0].id:
            # mutuples do not have to be flipped
            print("[+] Ancestor is %s ..." % self.ancestor_population_id)
        elif self.ancestor_population_id == entityCollection.populationObjs[1].id:
            mutype_df.rename(columns={'hetA': 'hetB', 'hetB': 'hetA'}, inplace=True)
            print("[+] Ancestor is %s (hetA and hetB will be flipped)... " % self.ancestor_population_id)
        mutuple_count_matrix_by_window_id = {}
        shape = tuple(self.max_by_mutype[mutype] + 2 for mutype in MUTYPES)

        lumped_df_header = ['window_id', 'count', 'hetA', 'fixed', 'hetB', 'hetAB']
        lumped_df_vals = []
        mutype_df = mutype_df.drop(mutype_df[(mutype_df.fixed > 0) & (mutype_df.hetAB > 0)].index) # exclude FGVs
        total_counts_by_window_id = mutype_df[['window_id', 'count']].groupby(['window_id']).sum().to_dict()['count']
        #print(total_counts_by_window_id)
        for window_id, count, hetA, fixed, hetB, hetAB in tqdm(mutype_df[['window_id', 'count'] + MUTYPES].values, total=len(mutype_df.index), desc="[%] ", ncols=100):
            if not window_id in mutuple_count_matrix_by_window_id:
                mutuple_count_matrix_by_window_id[window_id] = np.zeros(shape, np.float64)
            mutuple = (hetA, fixed, hetB, hetAB)
            normalised_count = count / total_counts_by_window_id[window_id]
            mutuple_vector = tuple(
                    [_count if not _count > self.max_by_mutype[mutype] else self.max_by_mutype[mutype] + 1 for _count, mutype in zip(mutuple, MUTYPES)])
                #print(mutuple, mutuple_vector)
            mutuple_count_matrix_by_window_id[window_id][mutuple_vector] += normalised_count
                #print(mutuple_count_matrix_by_window_id)
            lumped_df_vals.append([window_id, normalised_count] + list(mutuple_vector))
        #create_csv('lumped.csv', lumped_df_header, lumped_df_vals, ",")
        return mutuple_count_matrix_by_window_id

    def check_parameters(self):
        # checking params and setting symbols
        required_parameters = PARAMETERS_BY_MODEL_NAME[self.model_name]
        missing_parameters = []
        if 'Migration' in required_parameters:
            # migration rate per lineage (do the bounds also have to be divided?) 
            self.parameters_symbolic['Migration'] = sympy.Symbol('Migration', positive=True)
        if 'C_derived' in required_parameters:
            self.parameters_symbolic['C_derived'] = sympy.Symbol('C_derived', positive=True)
        if 'theta' in required_parameters:
            if (self.theta_low and self.theta_high):
                if not self.theta_low < self.theta_high:
                    sys.exit("[X] Lower bound must be smaller than upper bound: %s" % 'theta')
            self.parameters_symbolic['theta'] = sympy.Symbol('theta', positive=True)
            self.parameters_symbolic['hetA'] = sympy.Symbol('hetA', positive=True)
            self.parameters_symbolic['fixed'] = sympy.Symbol('fixed', positive=True)
            self.parameters_symbolic['hetB'] = sympy.Symbol('hetB', positive=True)
            self.parameters_symbolic['hetAB'] = sympy.Symbol('hetAB', positive=True)
        if 'Time' in required_parameters:
            self.parameters_symbolic['Time'] = sympy.Symbol('Time', positive=True)
        if missing_parameters:
            sys.exit("[X] Please specify values or lower/upper bounds for parameter(s): %s" % ",".join(missing_parameters))
        if self.ancestor_population_id is None:
            sys.exit("[X] Please specify a population id using '--ancestor_population'.")
        self.parameters_symbolic['C_ancestor'] = sympy.Symbol('C_ancestor', positive=True)
        self.parameters_symbolic['BigL'] = sympy.Symbol('BigL', positive=True)
        self.parameters_numeric['BigL'] = sympy.Symbol('BigL', positive=True)
        if self.windows_file:
            self.mode = 'windows'
        
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

#def inverse_laplace_transform(params):
#    vector, marginal_query, equation, rate_by_mutype, split_time = params
#    equation_giac = str(equation).replace("**", "^")
#    mutation_substitution = ",[%s]" % ",".join(["%s=%s" % (mutype, rate) for mutype, rate in rate_by_mutype.items()])
#    assumptions = "assume(BigL >= 0)"
#    '''
#    [ratnormal] rewrites an expression using its irreducible representation. The expression is viewed as a multivariate rational fraction with coefficients in Q (or
#        Q[i]). The variables are generalized identifiers which are assumed to be algebraically independent. Unlike with normal, an algebraic extension is considered
#        as a generalized identifier. Therefore ratnormal is faster but might miss some
#        simplifications if the expression contains radicals or algebraically dependent transcendental functions.
#    '''
#    invlaplace_string = "%s; invlaplace(ratnormal(subst(%s%s) / BigL), BigL, T).subst(T, %s)" % (assumptions, equation_giac, mutation_substitution, float(split_time))
#    try:
#        process = subprocess.run(["giac", invlaplace_string], stderr=subprocess.PIPE, stdout=subprocess.PIPE, check=True, encoding='utf-8')
#        probability = process.stdout.rstrip("\n").split(",")[1]
#    except subprocess.CalledProcessError:
#        exit("[X] giac could not run.")
#    #if probability == 'undef':
#    #    print(invlaplace_string)
#    return sympy.Rational(str(probability)), vector, marginal_query

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

def calculate_probability_matrix(params):
    grid_idx, grid_params, vectors, marginal_queries, equations_by_vector, parameters_symbolic, max_by_mutype = params 
    numeric_value_by_symbol_for_substitution = {}
    theta_by_mutype = {}
    theta = None
    Time = None
    for parameter, value in zip(GRID_PARAMETER_ORDER, grid_params):
        if parameter == 'theta':
            theta = sympy.Rational(str(value))
        else:
            symbol = parameters_symbolic[parameter]
            numeric_value_by_symbol_for_substitution[symbol] = sympy.Rational(str(value))
            if parameter == 'Time':
                Time = str(value)
    giac_strings = []
    for vector, marginal_query in zip(vectors, marginal_queries):
        print(vector)
        if marginal_query is None:
            for mutype in MUTYPES:
                theta_by_mutype[mutype] = theta
                symbol = parameters_symbolic[mutype]
                numeric_value_by_symbol_for_substitution[symbol] = theta
        else:
            for mutype, value in zip(MUTYPES, marginal_query):
                if isinstance(value, int):
                    theta_by_mutype[mutype] = theta
                    symbol = parameters_symbolic[mutype]
                    numeric_value_by_symbol_for_substitution[symbol] = theta
                else:
                    theta_by_mutype[mutype] = 0
                    symbol = parameters_symbolic[mutype]
                    numeric_value_by_symbol_for_substitution[symbol] = 0
        print(numeric_value_by_symbol_for_substitution)
        equation = equations_by_vector[vector].xreplace(numeric_value_by_symbol_for_substitution)
        giac_string = []
        giac_string.append("normal(invlaplace(normal(")
        equation_giac = str(equation).replace("**", "^")
        giac_string.append(equation_giac)
        giac_string.append(") / BigL, BigL, Time)).subst(Time, ")
        giac_string.append(Time)
        giac_string.append(");")
        giac_strings.append("".join(giac_string))
    probabilities = None
    with tempfile.NamedTemporaryFile(mode='w') as temp_fh:
        temp_fh.write("%s\n" % "\n".join(giac_strings))
        try:
            process = subprocess.run(["giac", temp_fh.name], stderr=subprocess.PIPE, stdout=subprocess.PIPE, check=True, encoding='utf-8')
        except subprocess.CalledProcessError:
            exit("[X] giac could not run.")
        probabilities = [np.float64(sympy.Rational(str(probability))) for probability in process.stdout.split(",\n")]
    shape = tuple(max_by_mutype[mutype] + 2 for mutype in MUTYPES)
    probability_matrix = np.zeros(shape, np.float64)
    lines = []
    for vector, marginal_query, probability in zip(vectors, marginal_queries, probabilities):
        if probability == 'undef':
            print(vector, marginal_query, probability)
        if marginal_query:
            probability -= sum(probability_matrix[marginal_query].flatten())
        probability_matrix[vector] = probability
        lines.append("%s\t%s\t%s" % (grid_idx, vector, probability))
    return (grid_idx, probability_matrix, lines)

def compute_composite_likelihoods(mutuple_count_matrix_by_window_id, parameterObj):
    print("[+] Calculating likelihoods of model parameters ...")
    composite_likelihood_by_parameter_tuple_by_window_id = collections.defaultdict(dict)
    for window_id in tqdm(mutuple_count_matrix_by_window_id, total=len(mutuple_count_matrix_by_window_id), desc="[%] ", ncols=100):
        mutuple_count_matrix = mutuple_count_matrix_by_window_id[window_id]
        for parameter_tuple, probability_matrix in parameterObj.probability_matrix_by_parameter_tuple.items():
            #composite_likelihood = np.sum((xlogy(np.sign(probability_matrix), probability_matrix) * mutuple_count_matrix))
            composite_likelihood = np.sum((xlogy(np.sign(probability_matrix), probability_matrix) * mutuple_count_matrix))
            composite_likelihood_by_parameter_tuple_by_window_id[window_id][parameter_tuple] = composite_likelihood
            
    return composite_likelihood_by_parameter_tuple_by_window_id

def infer_composite_likelihoods(mutuple_count_matrix_by_window_id, parameterObj):
    composite_likelihood_by_parameter_tuple_by_window_id = collections.defaultdict(dict)
    for window_id in tqdm(mutuple_count_matrix_by_window_id, total=len(mutuple_count_matrix_by_window_id), desc="[%] ", ncols=100):
        mutuple_count_matrix = mutuple_count_matrix_by_window_id[window_id]
        for grid_idx, probability_matrix in enumerate(probability_matrix_by_grid_idx):
            probability_matrix = probability_matrix_by_grid_idx[grid_idx]
            composite_likelihood = np.sum((xlogy(np.sign(probability_matrix), probability_matrix) * mutuple_count_matrix))
            composite_likelihood_by_grid_idx_by_window_id[window_id][grid_idx] = composite_likelihood

def generate_output(composite_likelihood_by_parameter_tuple_by_window_id, parameterObj):
    composite_likelihood_cols = ['window_id', 'cL', 'theta_ancestral', 'theta_derived', 'migration', 'ismax']
    composite_likelihood_vals = []
    print("[+] Generating output ...")
    for window_id in tqdm(composite_likelihood_by_parameter_tuple_by_window_id, total=len(composite_likelihood_by_parameter_tuple_by_window_id), desc="[%] ", ncols=100):
        max_parameter_tuple = max(composite_likelihood_by_parameter_tuple_by_window_id[window_id], key=composite_likelihood_by_parameter_tuple_by_window_id[window_id].get)
        for parameter_tuple, composite_likelihood in composite_likelihood_by_parameter_tuple_by_window_id[window_id].items():
            ismax = False
            if parameter_tuple == max_parameter_tuple:
                #print(window_id, parameter_tuple, composite_likelihood)
                ismax = True
            theta_ancestral = str(parameter_tuple[0])
            theta_derived = str(parameter_tuple[1])
            migration = str(parameter_tuple[2])
            composite_likelihood_vals.append([window_id, composite_likelihood, theta_ancestral, theta_derived, migration, ismax])
    composite_likelihood_df = pd.DataFrame(composite_likelihood_vals, columns=composite_likelihood_cols)
    plot_parameter_scan(composite_likelihood_df, parameterObj)
    out_f = "%s.composite_likelihoods.h5" % parameterObj.prefix
    composite_likelihood_hdf5_store = create_hdf5_store(out_f=out_f, path=parameterObj.path)
    composite_likelihood_df.to_hdf(composite_likelihood_hdf5_store, 'likelihood', append=True)
    composite_likelihood_hdf5_store.close()

# def generate_output(composite_likelihood_by_parameter_tuple_by_window_id, parameterObj):
#     composite_likelihood_cols = ['window_id', 'grid_idx', 'cL', 'raw_params', 'gimble_params', 'ismax']
#     composite_likelihood_vals = []
#     for window_id in tqdm(composite_likelihood_by_grid_idx_by_window_id, total=len(composite_likelihood_by_grid_idx_by_window_id), desc="[%] ", ncols=100):
#         max_grid_idx = max(composite_likelihood_by_grid_idx_by_window_id[window_id], key=composite_likelihood_by_grid_idx_by_window_id[window_id].get)
#         for grid_idx, composite_likelihood in composite_likelihood_by_grid_idx_by_window_id[window_id].items():
#             ismax = False
#             if grid_idx == max_grid_idx:
#                 ismax = True
#             raw_string = [str(x) for x in parameterObj.grid_raw[grid_idx]]
#             gimble_string = [str(x) for x in parameterObj.grid_gimble[grid_idx]]
#             composite_likelihood_vals.append([window_id, grid_idx, composite_likelihood, ", ".join(raw_string), ", ".join(gimble_string), ismax])
#     composite_likelihood_df = pd.DataFrame(composite_likelihood_vals, columns=composite_likelihood_cols)
#     print(composite_likelihood_df)
#     out_f = "%s.composite_likelihoods.h5" % parameterObj.prefix
#     composite_likelihood_hdf5_store = create_hdf5_store(out_f=out_f, path=parameterObj.path)
#     composite_likelihood_df.to_hdf(composite_likelihood_hdf5_store, 'likelihood', append=True)
#     composite_likelihood_hdf5_store.close()

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

def vectorise_equations(symbolic_equations_by_mutuple, max_by_mutype):
    equation_by_mutuple = {mutuple: sum(equations) for mutuple, equations in symbolic_equations_by_mutuple.items()}
    equations_by_vector = collections.defaultdict(list)
    marginal_queries = []
    vectors = []
    vector_seen = set([])
    for mutypes in reversed(list(powerset(MUTYPES))):
        mutypes_masked = set([mutype for mutype in MUTYPES if not mutype in mutypes])
        if not FOUR_GAMETE_VIOLATION.issubset(mutypes_masked):
            for data_tuple in sorted(equation_by_mutuple.keys()):               
                vector = tuple([count if not mutype in mutypes_masked else max_by_mutype[mutype] + 1 for count, mutype in zip(data_tuple, MUTYPES)])
                if not vector in vector_seen:
                    vector_seen.add(vector)
                    if not len([mutype for mutype, count in zip(MUTYPES, vector) if mutype in FOUR_GAMETE_VIOLATION and count > 0]) > 1: 
                        if not mutypes_masked:
                            marginal_query = None
                        else:
                            marginal_query = tuple(
                                [(data_tuple[idx] if not mutype in mutypes_masked else slice(0, max_by_mutype[mutype] + 2)) 
                                    for idx, mutype in enumerate(MUTYPES)])
                        vectors.append(vector)
                        marginal_queries.append(marginal_query)
                        equations_by_vector[vector] = equation_by_mutuple[data_tuple]
    return (vectors, marginal_queries, equations_by_vector)

def score_grid(symbolic_equations_by_mutuple, mutuple_count_matrix_by_window_id, parameterObj):
    # (C_ancestor, C_derived, Migration, theta, Time)
    print("[+] Vectorising equations ...")
    vectors, marginal_queries, equations_by_vector = vectorise_equations(symbolic_equations_by_mutuple, parameterObj.max_by_mutype)
    params = [(grid_idx, grid_params, vectors, marginal_queries, equations_by_vector, parameterObj.parameters_symbolic, parameterObj.max_by_mutype) for grid_idx, grid_params in enumerate(parameterObj.grid_gimble)]
    probability_matrix_by_grid_idx = {}
    print("[+] Calculating probability matrices ...")
    if parameterObj.threads < 2:
        for param in tqdm(params, total=len(params), desc="[%] ", ncols=100):
            grid_idx, probability_matrix, lines = calculate_probability_matrix(param)
            probability_matrix_by_grid_idx[grid_idx] = probability_matrix
            with open("result.%s.%s.txt" % (grid_idx, "_".join(str(param) for param in [parameterObj.grid_raw[grid_idx]])), 'w') as fh:
                fh.write("%s\n" % "\n".join(lines))
    else:
        with poolcontext(processes=parameterObj.threads) as pool:
            with tqdm(total=len(params), desc="[%] ", ncols=100, unit_scale=True) as pbar:
                for grid_idx, probability_matrix, lines in pool.imap_unordered(calculate_probability_matrix, params):
                    with open("result.%s.%s.txt" % (grid_idx, "_".join(str(param) for param in [parameterObj.grid_raw[grid_idx]])), 'w') as fh:
                        fh.write("%s\n" % "\n".join(lines))
                    probability_matrix_by_grid_idx[grid_idx] = probability_matrix
                    pbar.update()
    return probability_matrix_by_grid_idx
                
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
        print_params = { (param):(value if not param == 'theta' else value * 2) for param, value in print_params.items()}
        estimated_parameters_string = ", ".join(["%s=%s" % (key, round(value, 4)) for key, value in estimated_parameters.items()])
        print("[+] Parameters estimated in %ss using %s iterations (Composite Likelihood = -%s): %s" % (timer() - start_time, res.nit, res.fun, estimated_parameters_string))
    else:
        print("[-] No covergence reached after %s iterations (%ss elapsed)" % (res.nit, timer() - start_time))
    return estimated_parameters
