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
'''
mutype = ['m_1', 'm_2', 'm_3', 'm_4'] (for n2+p2)
mutation_profile = [1, 1, 0, 1]
event_tuple = (event, demography_counter, mutation_counter) # step in path of stategraph
constructor = datastructure for terms of equation (used for building equations)

To Do:
    - general way to filter FGVs in generate_mutation_profiles()
'''
@contextlib.contextmanager
def poolcontext(*args, **kwargs):
    pool = multiprocessing.Pool(*args, **kwargs)
    yield pool
    pool.terminate()

def get_data_array(parameterObj):
    store = lib.gimble.load_store(parameterObj)
    print('[++]')
    print(store.tree())
    store.check_existing_data(parameterObj, 'blocks', fail=True)
    if parameterObj.data_type == 'blocks':
        store.check_existing_data(parameterObj, 'blocks', fail=True)
        '''sum all block counts across all pairs'''
        '''visualise somehow'''
        pass
    elif parameterObj.data_type == 'windows':
        store.check_existing_data(parameterObj, 'windows', fail=True)
        '''sum all block counts across all pairs'''
    else:
        sys.exit("[X2] This should never happen.")

    #mutype_hdf5_store = pd.HDFStore(infile)
    #mutype_df = pd.read_hdf(mutype_hdf5_store, key='mutypes')
    #mutype_hdf5_store.close()
    #shape = tuple(self.max_by_mutype[mutype] + 2 for mutype in MUTYPES)
    #mutuple_count_matrix = np.zeros(shape, np.float64)
    ##print(self.ancestor_population_id, (entityCollection.populationObjs[0].id, entityCollection.populationObjs[1].id))
    ##print("before")
    ##print(mutype_df)
    #if self.ancestor_population_id == entityCollection.populationObjs[0].id:
    #    # mutuples do not have to be flipped
    #    print("[+] Ancestor is %s ..." % self.ancestor_population_id)
    #elif self.ancestor_population_id == entityCollection.populationObjs[1].id:
    #    mutype_df.rename(columns={'hetA': 'hetB', 'hetB': 'hetA'}, inplace=True)
    #    print("[+] Ancestor is %s (hetA and hetB will be flipped)... " % self.ancestor_population_id)
    ##print("before")
    ##print(mutype_df)
    ## this has to be changed if order of mutypes changes
    #FGV_count = 0
    #kmax_binned_count = 0
    #total_count = mutype_df['count'].sum()
    #for count, hetA, fixed, hetB, hetAB in tqdm(mutype_df[['count'] + MUTYPES].values, total=len(mutype_df.index), desc="[%] ", ncols=100):
    #    #print(hetA, fixed, hetB, hetAB)
    #    mutuple = (hetA, fixed, hetB, hetAB)
    #    if mutuple[1] > 0 and mutuple[3] > 0:
    #        FGV_count += count  
    #    if any([count > self.kmax for count in mutuple]):
    #        kmax_binned_count += count
    #    mutuple_vector = tuple([count if not count > self.max_by_mutype[mutype] else self.max_by_mutype[mutype] + 1 for count, mutype in zip(mutuple, MUTYPES)])
    #    
    #    mutuple_count_matrix[mutuple_vector] += count
    #    #print(count, hetA, fixed, hetB, hetAB, mutuple_vector, mutuple_count_matrix)
    #print("[=] Total mutuple count = %s" % (format_count(total_count)))
    #print("[=] Counts excluded due to four-gamete-violations = %s (%s)" % (format_count(FGV_count), format_percentage(FGV_count / total_count)))
    #print("[=] Counts binned due to kmax = %s (%s)" % (format_count(kmax_binned_count), format_percentage(kmax_binned_count / total_count)))
    #return mutuple_count_matrix


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
        if constructor.is_compatible_with_mutation_profile(mutation_profile):
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

def calculate_inverse_laplace(params):
    time = sage.all.SR.var('T')
    equationObj, rates, split_time, dummy_variable = params
    equation = (equationObj.equation / dummy_variable).substitute(rates)
    equationObj.result = sage.all.inverse_laplace(equation, dummy_variable, time, algorithm='giac').substitute(T=split_time)
    return equationObj

#sage.all.numerical_approx(value, digits=1)

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

    def is_compatible_with_mutation_profile(self, mutation_profile):
        mutypes_to_place = set([mutype for mutype, count in mutation_profile.items() if count > 0])
        placeable = all([(True if self.placements_by_mutation[mutype] else False) for mutype in mutypes_to_place])
        return placeable

class EquationObj(object):
    def __init__(self, matrix_idx, marginal_idx, equation):
        self.matrix_idx = matrix_idx
        self.marginal_idx = marginal_idx
        self.equation = equation
        self.result = None

class EquationSystemObj(object):
    def __init__(self, parameterObj):
        '''
        get list of equations by mutuple
        start with 
            1. Join-model (all Ne's == 1)
            2. Migration-model (all Ne's == 1)
            3. IM-model
        - what's the proportion of rare blocks as one add pairs
        '''
        self.events = []
        self.threads = parameterObj.threads
        self.boundaries = parameterObj._config['boundaries']
        self.k_max_by_mutype = parameterObj._config['k_max']
        self.mutypes = sorted(self.k_max_by_mutype.keys())
        self.model_file = parameterObj.model_file
        self.event_tuples_by_idx = self._get_event_tuples_by_idx()
        self.user_rate_by_event = self._get_user_rate_by_event(parameterObj)
        self.base_rate_by_variable = self._get_base_rate_by_variable()
        self.split_time = self.user_rate_by_event['T']
        self.dummy_variable = self._get_dummy_variable()
        self.rate_by_event = self._get_rate_by_variable(prefix=set(['C', 'M']))
        self.rate_by_mutation = self._get_rate_by_variable(prefix=set(['m']))
        self.equation_batch_by_idx = collections.defaultdict(list)

    def info(self):
        print("[=] ==================================================")
        print("[+] Parsed parameters ...")
        print("[+] K_max := %s" % self.k_max_by_mutype)
        print("[+] User-provided rates := %s" % self.user_rate_by_event)
        print("[+] Event rates := %s" % self.rate_by_event)
        print("[+] Mutation rates := %s" % self.rate_by_mutation)
        print("[+] Split time (T) := %s" % self.split_time)
        print("[+] Dummy variable := %s" % self.dummy_variable)

    def _get_rate_by_variable(self, prefix=None):
        if not prefix is None:
            return {event: rate for event, rate in self.base_rate_by_variable.items() if str(event)[0] in prefix}
        return {}

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

    def _get_base_rate_by_variable(self):
        '''Substitution values used in equations. 
        Does not include ``T`` (since it gets substituted later?).
        Any adjustments to user rates happen here'''
        base_rate_by_variable = {}
        for event, rate in self.user_rate_by_event.items():
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
        if check_monomorphic:
            rates = {
                **{event: random.randint(1, 4) for event, rate in self.rate_by_event.items()}, 
                **{mutype: 0 for mutype, rate in self.rate_by_mutation.items()}
                }
            split_time = 1
            dummy_variable = self.dummy_variable
            params = (copy.deepcopy(self.equationObjs[0]), rates, split_time, dummy_variable)
            equationObj = calculate_inverse_laplace(params)
            if not equationObj.result == 1:
                sys.exit("[-] Monomorphic check failed: P(monomorphic) = %s (should be 1)" % equationObj.result)
            else:
                print("[+] Monomorphic check passed: P(monomorphic) = %s" % equationObj.result)

    def calculate_PODs(self):
        print("[=] ==================================================")
        print("[+] Calculating PODs ...")
        parameter_batches = []
        for equationObj in self.equationObjs:
            rates = {**self.rate_by_event, **self.rate_by_mutation}
            parameter_batches.append((equationObj, rates, self.split_time, self.dummy_variable))
        desc = "[%] Solving equations"
        equationObj_by_matrix_idx = {}
        if self.threads <= 1:
            for parameter_batch in tqdm(parameter_batches, desc=desc, ncols=100):
                equationObj = calculate_inverse_laplace(parameter_batch)
                equationObj_by_matrix_idx[equationObj.matrix_idx] = equationObj
        else:
            with poolcontext(processes=self.threads) as pool:
                with tqdm(parameter_batches, desc=desc, ncols=100) as pbar:
                    for resultObj in pool.imap_unordered(calculate_inverse_laplace, parameter_batches):
                        equationObj_by_matrix_idx[resultObj.matrix_idx] = resultObj
                        pbar.update()
        PODs = np.zeros(tuple(self.k_max_by_mutype[mutype] + 2 for mutype in self.mutypes), np.float64)
        for matrix_id, equationObj in equationObj_by_matrix_idx.items():
            if equationObj.marginal_idx is None:
                PODs[matrix_id] = equationObj.result
            else:
                PODs[matrix_id] = equationObj.result - sum(PODs[equationObj.marginal_idx].flatten())
        if not math.isclose(np.sum(PODs.flatten()), 1, rel_tol=1e-5):
            print("[-]\t∑(PODs) != 1 (rel_tol=1e-5)")
        else:
            print("[+]\t∑(PODs) == 1 ")
        return PODs

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
            print(constructor) 
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