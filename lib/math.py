import itertools
import sys
import sage.all
import pandas as pd
import collections
import numpy as np
import copy
import random
import math
from tqdm import tqdm
import multiprocessing
import contextlib
import lib.gimble
from scipy.special import xlogy

#sage.all.numerical_approx(value, digits=1)

'''
mutype = ['m_1', 'm_2', 'm_3', 'm_4'] (for n2+p2)
mutation_profile = [1, 1, 0, 1]
event_tuple = (event, demography_counter, mutation_counter) # step in path of stategraph
constructor = datastructure for terms of equation (used for building equations)

To Do:
    mutype = ['m_1', 'm_2', 'm_3', 'm_4'] (for n2+p2)
    used in 
    - general way to filter FGVs in generate_mutation_profiles()

# YAML (boundaries):
- *mu*: mutation rate in mutations per base per generation [optional]
- *theta*: mutation rate in coalescence time scale (using pi as unbiased estimator from data)
- *C_A*, *C_B*, *C_AB*: population size scalers (reference population should be set to 1.0) 
- *tau* split time in numbers of generation
- *m_e* probability of lineage to migrate per generation
- *t*: generations per year (optional, would trigger output in years) [optional]


ALWAYS report coalescence-scaled parameters
- Ne_A, Ne_B, Ne_A_B, T, M  

# gimble inference grid:
- Ne = theta/4mu 
- theta_block = theta * block_length
- Ne_A, Ne_B, Ne_A_B = Ne/C_A, Ne/C_B, Ne/C_A_B
- T = 2Ne*tau
- M = 4Ne*m_e
    
# numbers of generation
- *mu* in mutations per base per generation 
- *tau* split time in numbers of generation
- *m_e* probability of lineage to migrate per generation

# on coalescence time scale
- theta = 4Ne*mu
- T = 2Ne*tau
- M = 4Ne*m_e

Goal:
    - user provided *mu* in mutations per base per generation 
    - required for grid-search
        - m_e
        - Ne_A, Ne_B, Ne_A_B
        - tau (split time in numbers of generation)

YAML:
    - add mu
    - add theta
    - remove theta_block 

get_base_rate:
    - get theta_scaled : theta per block (theta * block_length) 
    - calculate Ne
    - C_A = 1, C_B = 0.5
        - Ne_A = Ne/C_A; Ne_B = Ne/C_B
    One scalar has to be set to 1, else error ('you need a reference population', 'can;t estimate 4 values')
    - M: 4Ne*m (this is currently used in yaml)
    - internal m_e = has to be used for grid 

'''
@contextlib.contextmanager
def poolcontext(*args, **kwargs):
    pool = multiprocessing.Pool(*args, **kwargs)
    yield pool
    pool.terminate()     

def block_mutype_counter(blockcounts, k_max_by_mutype, mutypes):
    A = np.zeros(tuple(k_max_by_mutype[mutype] + 2 for mutype in mutypes), np.float64)
    #does A contain floats rather than ints?
    counts, unique = blockcounts[:,0], blockcounts[:,1:]    
    capped = np.clip(counts, [0 for i in len(k_max_by_mutype)], kmax_by_mutype)
    #need to collapse values that are no longer unique, so we can no longer use the faster
    #A[tuple(unique.transpose)] = counts approach
    for mutype, count in zip(capped, counts): 
        A[mutype]+=count
    return A 
    
def multinomial(lst):
    '''https://stackoverflow.com/questions/46374185/does-python-have-a-function-which-computes-multinomial-coefficients'''
    res, i = 1, 1
    for a in lst:
        for j in range(1,a+1):
            res *= i
            res //= j
            i += 1
    return res

def get_mutation_profiles(k_max_by_mutype, offset=0, by_mutype=True):
    '''Once generalised for FGV filtering for other than n2+p2 it should be an independent function.'''
    mutation_profiles = []
    for mutation_tuple in itertools.product(*[range(0, k_max + offset + 1) for mutype, k_max in k_max_by_mutype.items()]):
        valid = True
        if len(k_max_by_mutype) == 4:
            # This only captures FGVs for n2+p2 has to be made general
            if (mutation_tuple[2] > 0 and mutation_tuple[3] > 0):
                valid = False
        if valid:
            if by_mutype:
                mutation_profile = dict(zip(k_max_by_mutype.keys(), mutation_tuple))
            else:
                mutation_profile = mutation_tuple
            mutation_profiles.append(mutation_profile)
    return mutation_profiles

def place_mutations(parameter_batch):
    mutation_profile, constructors = parameter_batch 
    mutation_tuple = tuple(mutation_profile.values())
    equations = []
    for constructor in constructors:
        equation = 0    
        if constructor.is_mutable_by(mutation_profile):
            mutation_equation = 0
            for placement_tuples in itertools.product(*[list(itertools.combinations_with_replacement(constructor.placements_by_mutation[mutype], count)) for mutype, count in mutation_profile.items()]):
                mutypes_by_idx = collections.defaultdict(list)
                for mutype, placement_tuple in zip(mutation_profile.keys(), placement_tuples):
                    if placement_tuple:
                        for idx in list(placement_tuple):
                            mutypes_by_idx[idx].append(mutype)
                mutation_part = 1
                for idx, _mutypes in mutypes_by_idx.items():
                    mutype_counter = collections.Counter(_mutypes)
                    mutation_part *= multinomial(mutype_counter.values())
                    for mutype in _mutypes:
                        mutation_part *= (sage.all.SR.var(mutype) * constructor.mutation_counters[idx][mutype]) / constructor.denominators[idx]
                mutation_equation += mutation_part
            event_equation = 1
            for numerator, denominator in zip(constructor.numerators, constructor.denominators):
                event_equation *= numerator / denominator
            equation = mutation_equation * event_equation
        equations.append(equation)
    return (mutation_tuple, sum(equations))

def param_generator(pcentre, pmin, pmax, psamples, distr):
    starts = [pmin, pcentre]
    ends = [pcentre, pmax]
    nums = [round(psamples/2) + 1, round(psamples/2)] if psamples % 2 == 0 else [round(psamples/2) + 1, round(psamples/2) + 1]
    if distr == 'linear':
        return np.unique(np.concatenate([np.linspace(start, stop, num=num, endpoint=True, dtype=np.float64) for start, stop, num in zip(starts, ends, nums)]))
    else:
        raise NotImplmentedError

def calculate_inverse_laplace(params):
    '''
    [To Do]
    - test for errors due to unsolve equations due to wrong model? 
    - or is there a way to sort out model-coherenece as pre-flight check? 
    '''
    equationObj, rates, split_time, dummy_variable = params
    equation = (equationObj.equation).substitute(rates)
    if split_time is None:
        equationObj.result = equation
    else:
        equationObj.result = sage.all.inverse_laplace(equation / dummy_variable, dummy_variable, sage.all.SR.var('T'), algorithm='giac').substitute(T=split_time)
    
    return equationObj

def calculate_composite_likelihood(ETPs, data):
    return -np.sum((xlogy(np.sign(ETPs), ETPs) * data))

class Constructor(object):
    def __init__(self, constructor_id):
        self.constructor_id = constructor_id
        self.numerators = []
        self.denominators = []
        self.placements_by_mutation = collections.defaultdict(list)
        self.mutation_counters = []
        self.mutations = []
    
    def __str__(self):
        return "\n".join([
            '[+] constructor_id %r' % self.constructor_id,
            '[+] numerators %r' % self.numerators,
            '[+] denominators %r' % self.denominators,
            '[+] placements_by_mutation %r' % self.placements_by_mutation,
            '[+] mutation_counters %r' % self.mutation_counters,
            '[+] mutations %r' % self.mutations])

    def is_mutable_by(self, mutation_profile):
        mutypes_to_place = set([mutype for mutype, count in mutation_profile.items() if count > 0])
        placeable = all([(True if self.placements_by_mutation[mutype] else False) for mutype in mutypes_to_place])
        return placeable

class EquationObj(object):
    def __init__(self, matrix_idx, marginal_idx, equation):
        self.matrix_idx = matrix_idx
        self.marginal_idx = marginal_idx
        self.equation = equation
        self.result = None

    def __repr__(self):
        return '[+] EquationObj : %s\n%s' % (self.matrix_idx, self.equation)

class EquationSystemObj(object):
    def __init__(self, parameterObj, legacy=False):
        '''
        Step 1 : 
            - parsing of model : self._get_event_tuples_by_idx
            - building of equations: self._get_equationObjs

        Grid search:
            - build grid from config
            - for each grid point:

            - 
        '''
        #parameters
        self.threads = parameterObj.threads
        self.k_max_by_mutype = parameterObj._config['k_max'] if legacy else parameterObj.config['k_max']
        self.mutypes = sorted(self.k_max_by_mutype.keys())
        #needed to generate the equationObjs
        self.model_file = parameterObj.model_file
        self.events = []
        self.event_tuples_by_idx = self._get_event_tuples_by_idx()
        self.dummy_variable = self._get_dummy_variable() 
        self.equation_batch_by_idx = collections.defaultdict(list)
        self.ETPs=None

        if not legacy:
            self.rate_by_variable = [self._get_base_rate_by_variable(params) for params in parameterObj.grid]
            self.split_times = [self._get_split_time(params) for params in parameterObj.grid]
        
        else:  
            #user provided rates, legacy code
            self.user_rate_by_event = self._get_user_rate_by_event(parameterObj)
            self.base_rate_by_variable = self._get_base_rate_by_variable(self.user_rate_by_event)
            self.split_times = [self.user_rate_by_event.get('T', None)]
            #self.boundaries = parameterObj._config['boundaries'] #this needs to be changed
            self.rate_by_event = self._get_rate_by_variable(prefix=set(['C', 'M']))
            self.rate_by_mutation = self._get_rate_by_variable(prefix=set(['m']))
            self.probcheck_file = parameterObj.probcheck_file
            #self.grid_points = self._get_grid_points(parameterObj) #should this contain all parameter combos?

    def _get_grid_points(self, parameterObj):
        '''parameterObj should already by scaled when this point is reached'''
        config = parameterObj.config
        if not parameterObj._MODULE == 'gridsearch':

            return None



    def check_ETPs(self):
        '''
        ./gimble inference -m models/gimble.model.A_B.p2.n_1_1.J_A_B.Div.tsv -c models/gimble.model.A_B.p2.n_1_1.J_A_B.Div1.config.yaml -t 1 -b -k models/gimble.model.A_B.p2.n_1_1.J_A_B.Div1.probabilities.csv
        ./gimble inference -m models/gimble.model.A_B.p2.n_1_1.J_A_B.Div.tsv -c models/gimble.model.A_B.p2.n_1_1.J_A_B.Div2B.config.yaml -t 1 -b -k models/gimble.model.A_B.p2.n_1_1.J_A_B.Div2B.probabilities.csv

        '''
        probabilities_df = pd.read_csv(self.probcheck_file, sep=",", names=['A', 'B', 'C', 'D', 'probability'],  dtype={'A': int, 'B': int, 'C': int, 'D': int, 'probability': float}, float_precision='round_trip')
        for a, b, c, d, probability in probabilities_df.values.tolist():
            mutuple = (int(a), int(b), int(c), int(d))
            print(mutuple, " : ", probability, self.ETPs[mutuple], np.isclose(probability, self.ETPs[mutuple], rtol=1e-15))

    def setup_grid(self):
        '''LEGACY CODE'''
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

    def info(self):
        print("[=] ==================================================")
        print("[+] Parsed parameters ...")
        print("[+] K_max := %s" % self.k_max_by_mutype)
        #print("[+] User-provided rates := %s" % self.user_rate_by_event)
        #print("[+] Event rates := %s" % self.rate_by_variable)
        #print("[+] Mutation rates := %s" % self.rate_by_mutation)
        print("[+] Split time (T) := %s" % self.split_times)
        print("[+] Dummy variable := %s" % self.dummy_variable)

    def _get_rate_by_variable(self, prefix=None):
        if not prefix is None:
            return {event: rate for event, rate in self.base_rate_by_variable.items() if str(event)[0] in prefix}
        return {}

    def _get_split_time(self, params):
        time = params.get('T')
        if time:
            time = sage.all.Rational(float(time)) 
        return time

    def _get_dummy_variable(self):
        '''CAUTION: only works for 2-pop case with single J_* event'''
        dummy_variable = [sage.all.SR.var(var) for var in self.events if var.startswith("J")]
        return dummy_variable[0] if dummy_variable else sage.all.SR.var('J')

    def _get_event_tuples_by_idx(self):
        print("[+] Reading model %r" % self.model_file)
        paths_df = pd.read_csv(self.model_file, sep="\t", comment='#')
        header = paths_df.columns.values.tolist()
        self.events = [s for s in header if s[0] in set(['M', 'C', 'J'])]
        event_tuples_by_idx = collections.defaultdict(list)
        for line in paths_df.values.tolist():
            v_by_k = {k: v for k, v in zip(header, line)}
            idx = v_by_k['path_idx']
            event = v_by_k['event']
            count = v_by_k['count']
            demography_counter = {k: v for k, v in v_by_k.items() if (k[0] in set(['M', 'C', 'J']) and v > 0)}
            mutation_counter = {k: v for k, v in v_by_k.items() if k[0] == 'm' and v > 0}
            event_tuple = (event, count, demography_counter, mutation_counter)
            event_tuples_by_idx[idx].append(event_tuple)
        return event_tuples_by_idx

    def _get_user_rate_by_event(self, parameterObj):
        user_rate_by_event = {}
        for event, rate in parameterObj._config['parameters'].items():
            user_rate_by_event[event] = sage.all.Rational(float(rate))
        return user_rate_by_event

    def _get_base_rate_by_variable(self, rateDict):
        '''Substitution values used in equations. 
        Does not include ``T`` (since it gets substituted later?).
        Any adjustments to user rates happen here'''
        base_rate_by_variable = {}
        for event, rate in rateDict.items():
            if event.startswith('theta'):
                for mutype in self.mutypes:
                    base_rate_by_variable[sage.all.SR.var(mutype)] = sage.all.Rational(rate * 0.5)         
            if event.startswith("C"):
                base_rate_by_variable[sage.all.SR.var(event)] = sage.all.Rational(rate)
            if event.startswith("M"):
                base_rate_by_variable[sage.all.SR.var(event)] = sage.all.Rational(rate * 0.5)
            if event.startswith("J"):
                base_rate_by_variable[sage.all.SR.var(event)] = sage.all.SR.var(event)
        return base_rate_by_variable


    def initiate_model(self, check_monomorphic=True):
        print("[=] ==================================================")
        print("[+] Initiating model ...")
        self.equationObjs = self._get_equationObjs()

        #if check_monomorphic:
        #    rates = {
        #        **{event: random.randint(1, 4) for event, rate in self.rate_by_event.items()}, 
        #        **{mutype: 0 for mutype, rate in self.rate_by_mutation.items()}
        #        }
        #    split_time = 1
        #    dummy_variable = self.dummy_variable
        #    params = (copy.deepcopy(self.equationObjs[0]), rates, split_time, dummy_variable)
        #    equationObj = calculate_inverse_laplace(params)
        #    if not equationObj.result == 1:
        #        sys.exit("[-] Monomorphic check failed: P(monomorphic) = %s (should be 1)" % equationObj.result)
        #    else:
        #        print("[+] Monomorphic check passed: P(monomorphic) = %s" % equationObj.result)

    def calculate_all_ETPs(self):
        #iterate over zip(self.rate_by_variable, self.split_times)
        self.ETPs = [] 
        for rates, split_time in zip(self.rate_by_variable, self.split_times):
            #print("HEEEERE", rates, split_time)
            self.ETPs.append(self.calculate_ETPs(rates, split_time))
        self.ETPs = np.array(self.ETPs)

    def calculate_ETPs(self, rates=None, split_time=None, threads=1):
        print("[=] ==================================================")
        print("[+] Calculating ETPs ...")
        #print(f'rates sage vars: {rates}')
        parameter_batches = []
        for equationObj in self.equationObjs:
            if rates is None: 
                rates = {**self.rate_by_event, **self.rate_by_mutation}
                split_time = self.split_times[0]
            #print((equationObj, rates, split_time, self.dummy_variable))
            parameter_batches.append((equationObj, rates, split_time, self.dummy_variable))
        desc = "[%] Solving equations"
        equationObj_by_matrix_idx = {}
        if threads <= 1:
            for parameter_batch in tqdm(parameter_batches, desc=desc, ncols=100):
                equationObj = calculate_inverse_laplace(parameter_batch)
                equationObj_by_matrix_idx[equationObj.matrix_idx] = equationObj
        else:
            # This does not work due to problem with multiprocessing library
            # replicate with `-t 2`
            '''
            pexpect.exceptions.ExceptionPexpect: isalive() encountered condition where "terminated" is 0, but there was no child process. Did someone else call waitpid() on our process?

            Maybe the solution is to import multiprocessing library under a different name? so that it does not clash
            '''
            with poolcontext(processes=threads) as pool:
                with tqdm(parameter_batches, desc=desc, ncols=100) as pbar:
                    for resultObj in pool.imap_unordered(calculate_inverse_laplace, parameter_batches):
                        equationObj_by_matrix_idx[resultObj.matrix_idx] = resultObj
                        pbar.update()
        ETPs = np.zeros(tuple(self.k_max_by_mutype[mutype] + 2 for mutype in self.mutypes), np.float64)
        for matrix_id, equationObj in equationObj_by_matrix_idx.items():
            if equationObj.marginal_idx is None:
                ETPs[matrix_id] = equationObj.result
            else:
                ETPs[matrix_id] = equationObj.result - sum(ETPs[equationObj.marginal_idx].flatten())
            print(matrix_id, ETPs[matrix_id])
        if not math.isclose(np.sum(ETPs.flatten()), 1, rel_tol=1e-5):
            print("[-] sum(ETPs) != 1 (rel_tol=1e-5)")
        else:
            print("[+] sum(ETPs) == 1 ")
        return ETPs

    def optimise_parameters(symbolic_equations_by_mutuple, mutuple_count_matrix, parameterObj):
        '''
        search-bounds vs. real-bounds vs parameterObj.boundaries
        '''

        print("[+] Optimising parameters: %s ..." % (", ".join(parameterObj.boundaries.keys())))    
        start_time = timer()
        simplex_values, simplex_parameters = generate_initial_simplex(parameterObj.boundaries, parameterObj.seed)
        x0 = tuple([0] * len(parameterObj.boundaries.keys()))
        #block_count = mutuple_count_matrix.flatten().sum()
        res = scipy.optimize.minimize(
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
    
    def _get_equationObjs(self):
        constructors = []
        for constructor_id, event_tuples in tqdm(self.event_tuples_by_idx.items(), total=len(self.event_tuples_by_idx), desc="[%] Building equations", ncols=100):
            constructor = Constructor(constructor_id)
            for idx, event_tuple in enumerate(event_tuples):
                event, count, demography_counter, mutation_counter = event_tuple
                constructor.numerators.append(sage.all.SR.var(event) * count)
                constructor.denominators.append(sum([(sage.all.SR.var(e) * c) for e, c in {**demography_counter, **mutation_counter}.items()]))
                for mutation_event in mutation_counter:
                    constructor.placements_by_mutation[mutation_event].append(idx)
                constructor.mutation_counters.append(mutation_counter) 
            constructors.append(constructor)
        # mutation placements
        parameter_batches = []
        mutation_profiles_offset_0 = get_mutation_profiles(self.k_max_by_mutype, offset=0, by_mutype=True)
        for mutation_profile in mutation_profiles_offset_0:
            parameter_batches.append(
                (mutation_profile, copy.deepcopy(constructors))
            )
        equation_by_mutation_tuple = {}
        if self.threads <= 1:
            for parameter_batch in tqdm(parameter_batches, desc="[%] Placing mutations", ncols=100):
                mutation_tuple, equations = place_mutations(parameter_batch)
                equation_by_mutation_tuple[mutation_tuple] = equations
        else:
            with poolcontext(processes=self.threads) as pool:
                with tqdm(parameter_batches, desc="[%] Placing mutations", ncols=100) as pbar:
                    for result_batch in pool.imap_unordered(place_mutations, parameter_batches):
                        mutation_tuple, equations = result_batch
                        equation_by_mutation_tuple[mutation_tuple] = equations
                        pbar.update()
        print("[+] Bundled equations along %s mutation tuples" % (len(equation_by_mutation_tuple)))
        # marginal tuples
        equationObjs = []
        mutation_profiles_offset_1 = get_mutation_profiles(self.k_max_by_mutype, offset=1, by_mutype=True)
        for mutation_profile in tqdm(mutation_profiles_offset_1, desc="[%] Accounting for marginal probabilities", ncols=100):
            matrix_idx = tuple(mutation_profile.values())
            equation_idx = matrix_idx
            marginal_idx = None 
            mutation_rates = {}
            if not matrix_idx in equation_by_mutation_tuple:
                # booleans : False := count <= kmax / Trues := count > kmax
                boolean_by_mutype = {mutype: (False if count <= self.k_max_by_mutype[mutype] else True) for mutype, count in mutation_profile.items()}
                equation_idx = tuple([(0 if boolean_by_mutype[mutype] else count) for mutype, count in mutation_profile.items()])    
                # marginal idx is a slice obj which contains all mutuples under its marginality 
                marginal_idx = tuple([(slice(0, self.k_max_by_mutype[mutype] + 2) if boolean else mutation_profile[mutype]) for mutype, boolean in boolean_by_mutype.items()]) 
                mutation_rates = {sage.all.SR.var(mutype): 0 for mutype, boolean in boolean_by_mutype.items() if boolean} 
            # equation_batch has all except rates and split_time 
            equationObj = EquationObj(
                                matrix_idx, 
                                marginal_idx, 
                                equation_by_mutation_tuple[equation_idx].substitute(mutation_rates))
            equationObjs.append(equationObj)

        print("[+] Generated equations for %s mutation tuples" % len(equationObjs))
        return equationObjs