import copy
import collections
import contextlib
import concurrent.futures
import datetime
import itertools
import math
import multiprocessing
import nlopt
import numpy as np
import pandas as pd
import random
import sys, os
import sage.all
import logging
from functools import partial
from functools import partialmethod
from timeit import default_timer as timer
from tqdm import tqdm

import lib.gimble
from lib.GeneratingFunction.gf import togimble


# [INFERENCE.py] : RENAME AS INFERENC.PY
#logging.basicConfig(filename='log_file.txt', level=logging.DEBUG)
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

NLOPT_EXIT_CODES = {
    1: 'optimum found', 
    2: 'stopvalue reached',
    3: 'tolerance on lnCL reached',
    4: 'tolerance on parameter vector reached',
    5: 'max number of evaluations reached',
    6: 'max computation time was reached'
    }

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
"""
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
"""
def param_generator(pcentre, pmin, pmax, psamples, distr):
    # [atavism]
    starts = [pmin, pcentre]
    ends = [pcentre, pmax]
    nums = [round(psamples/2) + 1, round(psamples/2)] if psamples % 2 == 0 else [round(psamples/2) + 1, round(psamples/2) + 1]
    if distr == 'linear':
        return np.unique(np.concatenate([np.linspace(start, stop, num=num, endpoint=True, dtype=np.float64) for start, stop, num in zip(starts, ends, nums)]))
    else:
        raise ValueError('"distr" must be "linear"')

def config_to_gf(config):
    sample_list = [('a','a'),('b','b'),()] #currently hard-coded
    # config['events']['coalescence'] has same order as config['populations']['pop_ids']
    # e.g. ['C_A', 'C_s', 'C_s'] if pop_ids = ['A', 'B', 'A_B'] and sync_pop_sizes = ['A_B', 'B']
    coalescence_rates = [sage.all.var(rate) for rate in config['events']['coalescence']]
    migration_direction = config['events'].get('migration', None)
    migration_rate = sage.all.var('M') if migration_direction else None
    exodus_direction = config['events'].get('exodus', None)
    exodus_rate = sage.all.var('J') if exodus_direction else None
    mutype_labels = list(config['k_max'].keys())
    gf = togimble.get_gf(sample_list, coalescence_rates, mutype_labels, migration_direction, migration_rate, exodus_direction, exodus_rate)
    return gf

def config_to_gf_before(model, mutype_labels, sync_pops=None):
    pop_ids, pop_ids_sync, events = lib.gimble.get_model_params(model)
    sample_list = [(),('a','a'),('b','b')] #currently hard-coded
    pop_mapping = ['A_B', 'A', 'B'] #pop_mapping: should be more general sorted list of pops
    coalescence_rates = [sage.all.var(f'C_{pop}') for pop in pop_mapping]
    if sync_pops:
        syncing_to, to_be_synced = sync_pops
        syncing_to_idx = pop_mapping.index(syncing_to)
        for pop_to_sync in to_be_synced:
            coalescence_rates[pop_mapping.index(pop_to_sync)] = coalescence_rates[syncing_to_idx]
    migration_events = [event for event in events if event.startswith('M')]
    exodus_events = [event for event in events if event.startswith('J')]
    if len(migration_events)>0:
        migration_rate = sage.all.var('M')
        migration_direction = [tuple(pop_mapping.index(pop) for pop in mig.lstrip('M_').split('_')) for mig in migration_events]
        #migration_direction = [(1,2)] if migration_events[0] == 'M_A_B' else [(2,1)]
    else:
        migration_rate, migration_direction = None, None
    if len(exodus_events)>0:
        exodus_rate = sage.all.var('J')
        exodus_source = [[pop_mapping.index(pop) for pop in exodus.lstrip('J_').split('_')] for exodus in exodus_events] 
        exodus_dest = [[pop_mapping.index(exodus.lstrip('J_'))] for exodus in exodus_events]
        exodus_direction = [tuple(srce + dest) for srce, dest in zip(exodus_source, exodus_dest)]
        #exodus_direction = [(1,2,0)]
    else:
        exodus_rate, exodus_direction = None, None
    gf = togimble.get_gf(sample_list, coalescence_rates, mutype_labels, migration_direction, migration_rate, exodus_direction, exodus_rate)
    return gf
"""
def calculate_inverse_laplace(params):
    '''
    [To Do]
    - test for errors due to unsolve equations due to wrong model? 
    - or is there a way to sort out model-coherenece as pre-flight check? 
    '''
    #precision = 165 #the number of bits used to represent the mantissa of a floating-point number.
    equationObj, rates, split_time, dummy_variable, precision = params
    #print(equationObj, rates, split_time, dummy_variable)
    equation = (equationObj.equation).substitute(rates)
    # GB: for testing purposes only: this needs to be done elsewhere!
    # DRL: where? you can only tell after substituting. And that's here
    # GB: here we wanted to verify whether all parameters are correctly passed by nlopt
    # however, we could catch these errors quicker by testing once the equations are built
    # whether they match with the number of parameters provided.
    assert(len(equation.arguments())<=1), "Parameters left unspecified: %s" % str(equation.arguments())
    try:
        if split_time is None:
            equationObj.result = equation
        else:
            equationObj.result = sage.all.inverse_laplace(equation / dummy_variable, dummy_variable, sage.all.SR.var('T', domain='real'), algorithm='maxima').substitute(T=split_time)
            assert 'ilt' not in str(equationObj.result), "Inverse laplace transform is undefined for equation: %s" % equation 
    except KeyboardInterrupt:
        print("Interrupted by user.")
        exit(-1)
    except AssertionError:
        with ('log_ilt.txt', 'w') as file:
            print(equationObj.result, file=file)
            sys.exit("Inverse laplace transform undefined (using maxima)")
    equationObj.result = sage.all.RealField(precision)(equationObj.result)
    return equationObj
"""


def run_single_optimize(p0, lower, upper, specified_objective_function, maxeval, xtol_rel, ftol_rel):
    opt = nlopt.opt(nlopt.LN_NELDERMEAD, len(p0))
    opt.set_lower_bounds(lower)
    opt.set_upper_bounds(upper)
    opt.set_max_objective(specified_objective_function)
    opt.set_xtol_rel(xtol_rel)
    if xtol_rel>0:
        #assigning weights to address size difference between params
        xtol_weights = 1/(len(p0) * p0)
        opt.set_x_weights(xtol_weights)
    opt.set_ftol_rel(ftol_rel)
    opt.set_maxeval(maxeval)
    optimum = opt.optimize(p0)
    rdict = {}
    rdict['lnCL'] = opt.last_optimum_value()
    rdict['optimum'] = optimum
    rdict['exitcode'] = opt.last_optimize_result()
    return rdict

def calculate_composite_likelihood(ETPs, data):
    ETP_log = np.zeros(ETPs.shape, dtype=np.float64)
    np.log(ETPs, where=ETPs>0, out=ETP_log)
    return np.sum(ETP_log * data)

# HOW DOES PATH work?
def objective_function(paramsToOptimize, grad, paramNames, fixedParams, mu, block_length, reference_pop, gfEvaluatorObj, data, verbose=True, path=None):
    if grad.size:
        raise ValueError('no optimization with derivatives implemented')
    start_time = timer()
    
    # needs:
    #
    params_to_optimize_unscaled = {k:v for k,v in zip(paramNames, paramsToOptimize)}
    all_rates_unscaled = {**params_to_optimize_unscaled, **fixedParams}
    all_rates_scaled = scale_single_gens_to_gimble_symbolic(
        all_rates_unscaled, 
        reference_pop, 
        block_length, 
        mu
        )

    ETPs = gfEvaluatorObj.evaluate_gf(all_rates_scaled, all_rates_scaled[sage.all.SR.var('theta')])
    result = calculate_composite_likelihood(ETPs, data)
    
    # can be array
    iteration_number = -1
    if isinstance(path, list):
        iteration_number = len(path)
        toSave = np.append(paramsToOptimize, (result, iteration_number))
        path.append(toSave)
    
    if verbose:
        #scale coalescence time scale: 
        #all_rates_coal_scaled = scale_single_gimble_to_coal(parameter_dict)
        # generation time scaled
        # for time use lib.gimble.format_time(seconds)
        time = str(datetime.timedelta(seconds=round(timer()-start_time, 1)))
        print("[+] i=%s -- {%s} -- L=%s -- %s" % (
            str(iteration_number).ljust(4), 
            " ".join(["%s=%s" % (k, '{:.2e}'.format(float(v))) for k, v in params_to_optimize_unscaled.items()]), 
            '{:.2f}'.format(float(result)), 
            str(time[:-5]) if len(time) > 9 else str(time))) # rounding is hard ...
    return result

def get_scaled_values_by_symbol(values_bounded_unscaled_by_parameter, config):
    # DRL 
    reference_pop_id = config['populations']['reference_pop_id']
    block_length = config['block_length']
    mu = config['mu']['mu']
    values_fixed_unscaled_by_parameter = {
        k: config['parameters_np'][k][0] for k in config['parameters_fixed']}
    values_unscaled_by_parameter = {
        **values_bounded_unscaled_by_parameter, 
        **values_fixed_unscaled_by_parameter}
    ref_Ne = values_unscaled_by_parameter.get('Ne_%s' % reference_pop_id, None)
    if not ref_Ne:
        raise ValueError('ref_Ne not found in %s' % str(values_unscaled_by_parameter))
    values_scaled_by_symbol = {}
    scaling_factor = 2
    values_scaled_by_symbol[sage.all.SR.var('theta')] = scaling_factor * sage.all.Rational(ref_Ne * mu * block_length)
    for parameter, value in values_unscaled_by_parameter.items():
        if parameter.startswith('Ne'):
            symbol = sage.all.SR.var('C_%s' % parameter.lstrip('Ne_'))
            value_scaled = sage.all.Rational(ref_Ne / value)
        elif parameter == 'T':
            symbol = sage.all.SR.var('T')
            value_scaled = sage.all.Rational(value / (scaling_factor * ref_Ne))
        elif parameter.startswith('m'):
            symbol = sage.all.SR.var('M')
            value_scaled = sage.all.Rational(scaling_factor*ref_Ne*value)
        else:
            raise ValueError('Unknown parameter %r with value %r' % (parameter, value))
        values_scaled_by_symbol[symbol] = value_scaled
    return values_scaled_by_symbol

def optimize_function(values_bounded_unscaled, grad, config, gfEvaluatorObj, data, verbose, track_likelihoods, track_values_unscaled):
    # DRL
    start_time = timer()
    values_bounded_unscaled_by_parameter = {
        k: v for k, v in zip(config['parameters_bounded'], values_bounded_unscaled)}
    values_scaled_by_symbol = get_scaled_values_by_symbol(values_bounded_unscaled_by_parameter, config)
    ETPs = gfEvaluatorObj.evaluate_gf(values_scaled_by_symbol, values_scaled_by_symbol[sage.all.SR.var('theta')])
    likelihood = calculate_composite_likelihood(ETPs, data)
    track_likelihoods.append(likelihood)
    track_values_unscaled.append(values_bounded_unscaled)
    if verbose:
        iteration = len(track_likelihoods)
        elapsed = lib.gimble.format_time(timer() - start_time)
        print("[+] i=%s -- {%s} -- L=%s -- %s" % (
            str(iteration).ljust(4), 
            " ".join(["%s=%s" % (k, '{:.2e}'.format(float(v))) for k, v in values_bounded_unscaled_by_parameter.items()]), 
            '{:.2f}'.format(float(likelihood)), 
            elapsed))
    return likelihood

def optimize_parameters(gfEvaluatorObj, data, config, verbose=True):
    # DRL
    #log_file = os.path.join(config['CWD'], "gimble.%s.%s.%s.log" % (
    #    config['gimble']['task'], 
    #    config['gimble']['model'], 
    #    config['gimble']['label']))
    #logging.basicConfig(
    #    level=logging.DEBUG, 
    #    format='%(asctime)s %(levelname)s %(message)s',
    #    filename=log_file, 
    #    filemode='w')
    track_likelihoods = []
    track_values_unscaled = []
    optimize_function_list = [partial(
                    optimize_function,
                    config=config,
                    gfEvaluatorObj=gfEvaluatorObj,
                    data=data,
                    verbose=verbose,
                    track_likelihoods=track_likelihoods,
                    track_values_unscaled=track_values_unscaled)]
    nlopt_results = [] 
    if config['num_cores'] <= 1:
        for specified_objective_function in tqdm(optimize_function_list, desc="progress", total=len(optimize_function_list), disable=verbose):
            single_run_result = run_single_optimize(
                p0=config['starting_points'][0], # has to be changed if always only one startpoint 
                lower=config['parameter_combinations_lowest'], 
                upper=config['parameter_combinations_highest'], 
                specified_objective_function=specified_objective_function, 
                maxeval=config['max_iterations'], 
                xtol_rel=config['xtol_rel'], 
                ftol_rel=config['ftol_rel']) 
            nlopt_results.append(single_run_result)
    else:
        # data has to be tested if it's an iterator
        run_list=[partial(
                run_single_optimize,
                p0=config['starting_points'][0],
                lower=config['parameter_combinations_lowest'],
                upper=config['parameter_combinations_highest'],
                specified_objective_function=specified_objective_function,
                maxeval=config['max_iterations'],
                xtol_rel=config['xtol_rel'],
                ftol_rel=config['ftol_rel']
                ) for specified_objective_function in optimize_function_list]
        with concurrent.futures.ProcessPoolExecutor(max_workers=config['num_cores']) as outer_pool:
            for idx, single_run in enumerate(tqdm(outer_pool.map(fp_map, run_list), desc="progress", total=len(run_list), disable=verbose)):
                nlopt_results.append(single_run)
    print(nlopt_results)
    print(track_likelihoods)
    print(track_values_unscaled)
    return nlopt_results

def optimize_parameters_old(gfEvaluatorObj, data, parameterObj, trackHistory=True, verbose=False, label='', param_combo_name=''):

    fixedParams = parameterObj._get_fixed_params(as_dict=True) #adapt for new config parser
    boundaryNames = _optimize_get_boundary_names(parameterObj, fixedParams, data)
    starting_points, parameter_combinations_lowest, parameter_combinations_highest = _optimize_get_boundaries(parameterObj.parameter_combinations, boundaryNames, parameterObj.numPoints)
    trackHistoryPath = [[] for _ in range(parameterObj.numPoints)]
    #specify the objective function to be optimized
    reference_pop = parameterObj.config['populations']['reference_pop']
    mu = parameterObj.config['mu']['mu']
    block_length = parameterObj.config['mu']['blocklength']
    specified_objective_function_list = _optimize_specify_objective_function(
        gfEvaluatorObj, 
        parameterObj.data_type, 
        trackHistoryPath, 
        data, 
        boundaryNames,
        fixedParams, 
        mu, 
        block_length, 
        reference_pop, 
        verbose
        )
    
    #print to screen
    startdesc = "[+] Optimization starting for specified point."
    if parameterObj.numPoints>1:
        startdesc = f"[+] Optimization starting for specified point and {parameterObj.numPoints-1} random points."
    if verbose:
        print(startdesc)
    #make log file for simulation replicates:
    _init_optimize_log(boundaryNames, f'optimize_log_{label}_{param_combo_name}.tsv', parameterObj._CWD, meta=parameterObj._get_cmd())
    
    #parallelize over starting points or data points
    allResults=[] 
    if parameterObj.gridThreads <= 1:
        with open(os.path.join(parameterObj._CWD, f'optimize_log_{label}_{param_combo_name}.tsv'), 'a') as log_file:
            for idx, (startPos, specified_objective_function) in enumerate(tqdm(zip(itertools.cycle(starting_points), specified_objective_function_list), desc="progress", total=len(specified_objective_function_list),disable=verbose)):
                single_run_result = run_single_optimize(startPos, parameter_combinations_lowest, parameter_combinations_highest, specified_objective_function, parameterObj.max_eval, parameterObj.xtol_rel, parameterObj.ftol_rel) 
                allResults.append(single_run_result)
                _optimize_log(single_run_result, idx, log_file)
                #allResults.append(run_single_optimize(startPos, parameter_combinations_lowest, parameter_combinations_highest, specified_objective_function, parameterObj.max_eval, parameterObj.xtol_rel, parameterObj.ftol_rel))
            
    else:
        #single_runs need to be specified before passing them to the pool
        specified_run_list = _optimize_specify_run_list_(parameterObj, specified_objective_function_list, starting_points, parameter_combinations_lowest, parameter_combinations_highest)
        with open(os.path.join(parameterObj._CWD, f'optimize_log_{label}_{param_combo_name}.tsv'), 'a') as log_file:
            with concurrent.futures.ProcessPoolExecutor(max_workers=parameterObj.gridThreads) as outer_pool:
                for idx, single_run in enumerate(tqdm(outer_pool.map(fp_map, specified_run_list), desc="progress", total=len(specified_run_list), disable=verbose)):
                    _optimize_log(single_run, idx, log_file)
                    allResults.append(single_run)

    #process results found in allResults (final step including exit code), and trackHistoryPath (all steps) 
    exitcodeDict = {
                    1: 'optimum found', 
                    2: 'stopvalue reached',
                    3: 'tolerance on lnCL reached',
                    4: 'tolerance on parameter vector reached',
                    5: 'max number of evaluations reached',
                    6: 'max computation time was reached'
                }
    result = _optimize_reshape_output(allResults, trackHistoryPath, boundaryNames, exitcodeDict, trackHistory, verbose)
    return result

def _optimize_get_boundary_names(parameterObj, fixedParams, data):
    syncing_to, to_be_synced = parameterObj._get_pops_to_sync_short()
    boundaryNames = [k for k in parameterObj.config['parameters'].keys() if not (k=='mu' or k in fixedParams or k in to_be_synced)]
    if len(boundaryNames) == 0:
        sys.exit("[X] No boundaries specified.")
    return boundaryNames
    
def _optimize_get_boundaries(parameter_combinations, boundary_names, numPoints):
    parameter_combinations_lowest = np.array([parameter_combinations[k][0] for k in boundary_names])
    parameter_combinations_highest = np.array([parameter_combinations[k][1] for k in boundary_names])
    #start from midpoint
    pmid = np.mean(np.vstack((parameter_combinations_lowest, parameter_combinations_highest)), axis=0)
    if numPoints > 1:
        #np.random.seed(self.seed)
        starting_points = np.random.uniform(low=parameter_combinations_lowest, high=parameter_combinations_highest, size=(numPoints-1, len(boundary_names)))
        starting_points = np.vstack((pmid,starting_points))
    else:
        starting_points = [pmid,]

    return (starting_points, parameter_combinations_lowest, parameter_combinations_highest)
    
def _optimize_specify_objective_function(gfEvaluatorObj, data_type, trackHistoryPath, data, boundary_names, fixedParams, mu, block_length, reference_pop, verbose):
    if data_type=='simulate':
        specified_objective_function_list = [partial(
                objective_function,
                paramNames=boundary_names,
                fixedParams=fixedParams,
                mu=mu,
                block_length = block_length,
                reference_pop=reference_pop,
                gfEvaluatorObj=gfEvaluatorObj,
                data=bsfs,
                verbose=verbose,
                path=None) for bsfs in data]
    
    elif data_type=='blocks':
        specified_objective_function_list = [partial(
                    objective_function,
                    paramNames=boundary_names,
                    fixedParams=fixedParams,
                    mu=mu,
                    block_length = block_length,
                    reference_pop=reference_pop,
                    gfEvaluatorObj=gfEvaluatorObj,
                    data=data,
                    verbose=verbose,
                    path=sublist) for sublist in trackHistoryPath]
    else:
        raise ValueError("Only blocks and sims have been implemented so far.")
    return specified_objective_function_list
    
def _optimize_specify_run_list_(parameterObj, specified_objective_function_list, starting_points, parameter_combinations_lowest, parameter_combinations_highest):
    specified_run_list=[partial(
                run_single_optimize,
                p0=p0,
                lower=parameter_combinations_lowest,
                upper=parameter_combinations_highest,
                specified_objective_function=specified_objective_function,
                maxeval=parameterObj.max_eval,
                xtol_rel=parameterObj.xtol_rel,
                ftol_rel=parameterObj.ftol_rel
                ) for p0, specified_objective_function in zip(itertools.cycle(starting_points), specified_objective_function_list)]
    return specified_run_list
    
def _optimize_reshape_output(raw_output, trackHistoryPath, boundary_names, exitcodeDict, trackHistory, verbose):
    #raw_output is list of dicts with {"lnCL":value, "optimum":list, "exitcode":value}
    if verbose:
        print([exitcodeDict.get(resultd['exitcode'],'Not in exitcodeDict.') for resultd in raw_output])
        
    #process trackhistory
    if not trackHistory:
        trackHistoryPath = [list(resultd['optimum'])+[resultd['lnCL'],label, resultd['exitcode']] for label,resultd in enumerate(raw_output)]
        result = [boundary_names+["lnCL", 'iterLabel', 'exitcode'],]+trackHistoryPath 
    else:
        trackHistoryPath = [list(ar)+[idx,] for idx, sublist in enumerate(trackHistoryPath) for ar in sublist]
        result = [boundary_names+['lnCL', 'step_id' , 'iterLabel'],]+trackHistoryPath

    return result

def _init_optimize_log(boundary_names, file_name, directory, meta=None):
    if not directory:
        directory=os.getcwd()
    header = boundary_names+["lnCL", 'iterLabel', 'exitcode']
    with open(os.path.join(directory, file_name), 'w') as logfile:
        if meta:
            print(meta, file=logfile)
        print('\t'.join(header), file=logfile)

def _optimize_log(single_run, identifier, log_file):
    single_run_string=list(single_run['optimum'])+[single_run['lnCL'], identifier, single_run['exitcode']]
    single_run_string = [str(value) for value in single_run_string]
    print('\t'.join(single_run_string), file=log_file)

def fp_map(f, *args):
    #used with pool.map(fp_map, function_list, arg_list1, arg_list2,...)
    return f(*args)

def scale_single_gens_to_gimble_symbolic(parameter_dict, reference_pop, block_length, mu):
    rdict = {}
    ref_value = parameter_dict[f'Ne_{reference_pop}']
    scaling_factor = 2
    rdict[sage.all.SR.var('theta')] = scaling_factor*sage.all.Rational(ref_value*mu*block_length)
    for parameter, value in parameter_dict.items():
        if parameter.startswith('Ne'):
            label = parameter.lstrip('Ne_')
            rdict[sage.all.SR.var(f'C_{label}')] = sage.all.Rational(ref_value/value)
        elif parameter=='T':
            rdict[sage.all.SR.var('T')] = sage.all.Rational(value/(scaling_factor*ref_value))
        elif parameter.startswith('m'):
            rdict[sage.all.SR.var('M')] = sage.all.Rational(scaling_factor*ref_value*value)
        else:
            pass
    return rdict

def scale_single_parameter_dict(parameter_dict, reference_pop, block_length, mu, input_scaling='gimble', output_scaling='generations', symbolic=False):
    parameter_dict = {k:[v,] for k,v in parameter_dict.items()}
    shape = 'LOD'
    LOD = scale_parameters(parameter_dict, reference_pop, block_length, mu, input_scaling, output_scaling, symbolic, 'LOD')
    return LOD[0]

def scale_parameters(parameter_dict, reference_pop, block_length, mu, input_scaling='generations', output_scaling='coalescence', symbolic=False, shape='DOL'):
    rdict = {}
    if input_scaling=='generations':
        reference_values = parameter_dict[f'Ne_{reference_pop}']
        if output_scaling=='coalescence' or output_scaling=='gimble':
            scaling_factor = 2
            rdict['theta'] = [scaling_factor*sage.all.Rational(ref*mu*block_length) for ref in reference_values]
            for parameter, values in parameter_dict.items():
                if parameter.startswith('Ne'):
                    label = parameter.lstrip('Ne_')
                    rdict[f'C_{label}'] = [sage.all.Rational(ref/other) for other, ref in zip(values, reference_values)]
                elif parameter=='T':
                    rdict['T'] = [sage.all.Rational(t/(scaling_factor*ref)) for t, ref in zip(values, reference_values)]
                elif parameter.startswith('m'):
                    rdict['M'] = [sage.all.Rational(scaling_factor*ref*m) for m, ref in zip(values, reference_values)]
                else:
                    pass
        else:
            raise ValueError(f'Scaling generations to {output_scaling} not implemented.')
    
    elif input_scaling=='gimble' or input_scaling=='coalescence':
        reference_values = parameter_dict[f'C_{reference_pop}']
        if output_scaling=='generations':
            scaling_factor = 2
            for parameter, values in parameter_dict.items():
                if parameter.startswith('C'):
                    label = parameter.lstrip('C_')
                    rdict[f'Ne_{label}'] = [float(ref/c) for c,ref in zip(values, reference_values)]
                elif parameter=='T':
                    rdict['T'] = [float(scaling_factor*t*ref) for t,ref in zip(values, reference_values)]
                elif parameter=='M':
                    rdict['m'] = [float(M/(scaling_factor*ref)) for M,ref in zip(parameter_values, reference_values)]
                else:
                    pass
        else:
            raise ValueError(f'Scaling {input_scaling} to {output_scaling} not implemented.')
    else:
        raise ValueError(f'Scaling {input_scaling} to {output_scaling} not implemented.')

    if symbolic:
        rdict = {sage.all.SR.var(k):v for k,v in rdict.items()}
    if shape == 'DOL':
        return rdict
    elif shape == 'LOD':
        return lib.gimble.DOL_to_LOD(rdict)
    else:
        raise ValueError(f'scale_parameters() not implemented for shape: {shape}')

def new_calculate_all_ETPs(gfEvaluatorObj, parameter_combinations, reference_pop, block_length, mu, processes=1, verbose=False):
    scaled_parameter_combinations = scale_parameters(
        parameter_combinations, 
        reference_pop, 
        block_length, 
        mu, 
        input_scaling='generations', 
        output_scaling='gimble', 
        symbolic=True, 
        shape='LOD'
        )
    all_ETPs = []
    desc = f'[%] Calculating mutation configuration probabilities for {len(scaled_parameter_combinations)} gridpoints'
    if processes==1:
        all_ETPs = [gfEvaluatorObj.evaluate_gf(parameter_combination, parameter_combination[sage.all.SR.var('theta')]) for parameter_combination in tqdm(scaled_parameter_combinations, desc=desc, ncols=100, disable=True)]
    else:
        args = ((param_combo, param_combo[sage.all.SR.var('theta')]) for param_combo in scaled_parameter_combinations)
        with multiprocessing.Pool(processes=processes) as pool:
            with tqdm(total=len(scaled_parameter_combinations), desc=desc, ncols=100) as pbar:
                    for ETP in pool.starmap(gfEvaluatorObj.evaluate_gf, args):
                        all_ETPs.append(ETP)
                        pbar.update()
    return np.array(all_ETPs, dtype=np.float64)
"""
class Constructor(object):
    # [EQUATIONS.py]
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
    # [EQUATIONS.py]
    def __init__(self, matrix_idx, marginal_idx, equation):
        self.matrix_idx = matrix_idx
        self.marginal_idx = marginal_idx
        self.equation = equation
        self.result = None

    def __repr__(self):
        return '[+] EquationObj : %s\n%s' % (self.matrix_idx, self.equation)

class EquationSystemObj(object):
    def __init__(self, model_file, reference_pop, k_max_by_mutype, block_length, mu, seed=None, module=None, threads=1, precision=165):
        # [EQUATIONS.py]
        #parameters
        self.model_file = model_file
        self.reference_pop = reference_pop
        self.threads = threads
        self.k_max_by_mutype = k_max_by_mutype
        
        self.block_length = sage.all.Rational(block_length)
        self.mu = mu
        np.random.seed(seed)
        self.mutypes = sorted(self.k_max_by_mutype.keys())
        #old init steps
        #self.threads = parameterObj.threads
        #self.k_max_by_mutype = parameterObj.config['k_max']
        #self.seed = np.random.seed(parameterObj.config['gimble']['random_seed'])
        #check whether this produces the desired result for example for 
        #self.mutypes = sorted(self.k_max_by_mutype.keys())
        #needed to generate the equationObjs
        #self.model_file = parameterObj.model_file
        
        self.events = []
        self.event_tuples_by_idx = self._get_event_tuples_by_idx()
        self.dummy_variable = self._get_dummy_variable() 
        self.equation_batch_by_idx = collections.defaultdict(list)
        self.ETPs=None
        self.precision = precision

        #self.scaled_parameter_combinations = self._scale_all_parameter_combinations(parameterObj)
        #self.rate_by_variable = [self._get_base_rate_by_variable(params) for params in self.scaled_parameter_combinations] 
        #self.split_times = [self._get_split_time(params) for params in self.scaled_parameter_combinations]
        
    def check_ETPs(self):
        probabilities_df = pd.read_csv(self.probcheck_file, sep=",", names=['A', 'B', 'C', 'D', 'probability'],  dtype={'A': int, 'B': int, 'C': int, 'D': int, 'probability': float}, float_precision='round_trip')
        for a, b, c, d, probability in probabilities_df.values.tolist():
            mutuple = (int(a), int(b), int(c), int(d))
            print(mutuple, " : ", probability, self.ETPs[mutuple], np.isclose(probability, self.ETPs[mutuple], rtol=1e-15))

    def info(self):
        print("[=] ==================================================")
        print("[+] Parsed parameters ...")
        print("[+] K_max := %s" % self.k_max_by_mutype)
        #print("[+] Split time (T) := %s" % self.split_times) #self.split_times no longer attribute
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
        return dummy_variable[0] if dummy_variable else sage.all.SR.var('J', domain='real')

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
        for event, rate in parameterObj.config['parameters'].items():
            user_rate_by_event[event] = sage.all.Rational(float(rate))
        return user_rate_by_event

    def _get_base_rate_by_variable(self, rateDict):
        '''Substitution values used in equations. 
        Does not include ``T`` (since it gets substituted later?).
        Any adjustments to user rates happen here'''
        base_rate_by_variable = {}
        for event, rate in rateDict.items():
            if event.startswith('theta'):
                #for mutype in self.mutypes:
                #    base_rate_by_variable[sage.all.SR.var(mutype)] = sage.all.Rational(rate * 0.5)
                base_rate_by_variable[sage.all.SR.var(event)] = sage.all.Rational(rate * 0.5)         
            if event.startswith("C"):
                base_rate_by_variable[sage.all.SR.var(event)] = sage.all.Rational(rate)
            if event.startswith("M"):
                base_rate_by_variable[sage.all.SR.var(event)] = sage.all.Rational(rate * 0.5)
            if event.startswith("J"):
                base_rate_by_variable[sage.all.SR.var(event)] = sage.all.SR.var(event, domain='real')
        return base_rate_by_variable

    def _scale_all_parameter_combinations(self, parameterObj):
        #function redundant
        #self.reference_pop = parameterObj.config['populations']['reference_pop']
        #self.block_length = sage.all.Rational(parameterObj.config['mu']['blocklength'])
        #result = []
    
        if parameterObj._MODULE == 'optimize':
            #currently not used, can be removed once function is rewritten.
            mid = self._scale_parameter_combination(parameterObj.parameter_combinations[0], self.reference_pop, self.block_length, parameterObj._MODULE)
            min_combo, max_combo = self._scale_min_max_parameter_combination(parameterObj.parameter_combinations[1], parameterObj.parameter_combinations[2], self.reference_pop, self.block_length, parameterObj)
            result = [mid, min_combo, max_combo]
            ratesPerVariable = {k:[d[k] for d in result] for k in result[0].keys()}
            ratesPerVariableSet = {k:len(set(v)) for k,v in ratesPerVariable.items()}
            self.fixed_params = []
            if len(parameterObj.fixed_params)>0:
                translation_dict = {'Ne':'C', 'me':'M'}
                self.fixed_params = parameterObj.fixed_params[:]
                for k,v in translation_dict.items(): 
                    self.fixed_params = [p.replace(k,v) for p in self.fixed_params] 
                for param in self.fixed_params:
                    assert ratesPerVariableSet[param]==1, 'The number of values should be one if a parameter is fixed.'
                #if Ne_ref fixed, theta should be fixed too
                if f'C_{self.reference_pop}' in self.fixed_params:
                    assert ratesPerVariableSet['theta']==1, 'No boundaries provided for the reference population, theta should be fixed as well.'
                    self.fixed_params.append('theta')
                #self.fixed_params = self._get_base_rate_by_variable()
            #all other parameters should have 3 values
            #assert all(ratesPerVariableSet[k]==3 for k in ratesPerVariable.keys() if k not in self.fixed_params or k not self.reference_pop), 'All parameters with boundaries should be associated with a list of length 3.'             
            print(self.fixed_params)
        
        if parameterObj._MODULE=='makegrid':
            #now dict of lists
            parameter_combinations = lib.gimble.DOL_to_LOD(parameterObj.parameter_combinations)
            result = [self._scale_parameter_combination(combo, self.reference_pop, self.block_length, parameterObj._MODULE, parameterObj.config['mu']['mu']) for combo in parameter_combinations]

        return result

    def _scale_parameter_combination(self, combo, reference_pop, block_length, module, mu):
        rdict = {}
        if module in ['makegrid', 'inference','optimize', None]:
            Ne_ref = sage.all.Rational(combo[f"Ne_{reference_pop}"])
            rdict['theta'] = 4*sage.all.Rational(Ne_ref*mu)*block_length
            if 'Ne_A' in combo:
                rdict['C_A']=Ne_ref/sage.all.Rational(combo['Ne_A'])
            if 'Ne_B' in combo:
                rdict['C_B'] = Ne_ref/sage.all.Rational(combo['Ne_B'])
            if 'Ne_A_B' in combo:
                rdict['C_A_B'] = Ne_ref/sage.all.Rational(combo['Ne_A_B'])
            if 'T' in combo:
                rdict['T'] = sage.all.Rational(combo['T'])/(2*Ne_ref)
            mig_dir = [key for key in combo.keys() if key.startswith("me_")]
            if mig_dir:
                mig_dir = mig_dir[0]
                mig_pop = mig_dir.lstrip("me_")
                rdict[f'M_{mig_pop}'] = 4*Ne_ref*sage.all.Rational(combo[mig_dir])
            return rdict
        else:
            ValueError("equationObj._scale_parameter_combination not implemented for this module.")

    def _unscale_all_parameter_combinations(self):
        pass

    def _unscale_parameter_combination(self, combo, reference_pop, block_length, mu):
        rdict = {}
        if parameterObj._MODULE in ['optimize']:
            #check if parameter needs to be scaled, e.g. not if already provided.
            theta/=blocklength 
            Ne_ref=theta/(2*mu) #remark 2's: theta here is 2*Ne*mu 
            m_e = M/(2*Ne_ref) #remark 2's: M is here 2*Ne*m_e
            tau = T/(2*Ne_ref)
            Ne_A = Ne_ref/C_A
            Ne_B = Ne_ref/C_B
            Ne_AB = Ne_ref/C_A_B
            return rdict
        else:
            sys.exit("[X] math.EquationSystemObj._unscale_parameter_combination not implemented for this module.")

    def initiate_model(self, sync_ref=None, sync_targets=None, check_monomorphic=True, disable_tqdm=False):
        print("[=] ==================================================")
        print("[+] Initiating model ...")
        self.equationObjs = self._get_equationObjs(sync_ref=sync_ref, sync_targets=sync_targets, disable_tqdm=disable_tqdm)
        #this is an essential model check: however, processing of ini file needs to be unified first
        #currently mu is missing from parameter_combinations in case of optimize for example.
        #for equationObj in self.equationObjs:
        #    if len(equationObj.equation.arguments)>len(parameterObj.parameter_combinations[0]):
        #       sys.exit("[X] Specified model contains more parameters than provided in ini file.")
        #check self.equationObjs
        # Method 2 (has to be tested)
        #if not (parameterObj and parameterObj.reference and parameterObj.toBeSynced):
        #   self.equationObjs = self._get_equationObjs(sync_ref=parameterObj.reference, sync_targets=parameterObj.toBeSynced) 
        #else:
        #    self.equationObjs = self._get_equationObjs() 

        # Method 1 
        #syncing pop sizes (coalescence rates) in the equation, there must be a better way
        #@Dom did not want to touch self._get_equationObjs() but could probably happen there
        # if parameterObj:
        #     if parameterObj.reference and parameterObj.toBeSynced:
        #         print(parameterObj.toBeSynced)
        #         for equationObj in self.equationObjs:
        #             for tBS in parameterObj.toBeSynced:
        #                 equationObj.equation = equationObj.equation.subs(sage.all.SR.symbol(f'C_{tBS}')==sage.all.SR.symbol(f'C_{parameterObj.reference}'))
        
    def calculate_all_ETPs(self, parameter_combinations, module=None, threads=1, gridThreads=1, verbose=False):
        #verboseprint = print if verbose else lambda *a, **k: None
        parameter_combinations = lib.gimble.DOL_to_LOD(parameter_combinations)
        scaled_parameter_combinations = [self._scale_parameter_combination(combo, self.reference_pop, self.block_length, module, self.mu) for combo in parameter_combinations]
        rate_by_variable = [self._get_base_rate_by_variable(params) for params in scaled_parameter_combinations] 
        split_times = [self._get_split_time(params) for params in scaled_parameter_combinations]
        ETPs = []
        desc = f'[%] Calculating likelihood for {len(rate_by_variable)} gridpoints'
        if gridThreads <= 1:
            for rates, split_time in tqdm(zip(rate_by_variable, split_times),total=len(rate_by_variable), desc=desc, ncols=100):
                arg = (rates, split_time, threads, verbose)
                ETPs.append(self.calculate_ETPs(arg))
        else:
            args = [(rates, split_time, threads, False) for rates, split_time in zip(rate_by_variable, split_times)]
            with concurrent.futures.ProcessPoolExecutor(max_workers=gridThreads) as outer_pool:
                with tqdm(total=len(args), desc=desc, ncols=100) as pbar:
                    for ETP in outer_pool.map(self.calculate_ETPs, args):
                        ETPs.append(ETP)
                        pbar.update()
        return np.array(ETPs)
    
    def calculate_ETPs(self, args):
        rates, split_time, threads, verbose = args
        #verboseprint = print if verbose else lambda *a, **k: None
        #verboseprint("[=] ==================================================")
        #verboseprint("[+] Calculating ETPs ...")
        #print(f'rates sage vars: {rates}')
        parameter_batches = []
        for equationObj in self.equationObjs:
            parameter_batches.append((equationObj, rates, split_time, self.dummy_variable, self.precision))
        desc = "[%] Solving equations"
        equationObj_by_matrix_idx = {}
        
        if threads <= 1:
            for parameter_batch in tqdm(parameter_batches, desc=desc, ncols=100, disable=not verbose):
                equationObj = calculate_inverse_laplace(parameter_batch)
                equationObj_by_matrix_idx[equationObj.matrix_idx] = equationObj        
        else:
            parameter_batches = [((pbatch,),{}) for pbatch in parameter_batches]
            #threads spawns a number of processes
            result = sage.parallel.multiprocessing_sage.parallel_iter(threads, calculate_inverse_laplace,parameter_batches)
            equationObj_by_matrix_idx = {equationObj.matrix_idx: equationObj for equationObj in [el[-1] for el in result]}
        
        #ETPs = np.zeros(tuple(self.k_max_by_mutype[mutype] + 2 for mutype in self.mutypes), np.float64)
        ETPs = np.zeros(tuple(self.k_max_by_mutype[mutype] + 2 for mutype in self.mutypes), dtype=object) # why not float?
        for matrix_id, equationObj in sorted(equationObj_by_matrix_idx.items()):
            if equationObj.marginal_idx is None:
                ETPs[matrix_id] = equationObj.result
            else:
                ETPs[matrix_id] = equationObj.result - sum(ETPs[equationObj.marginal_idx].flatten())
            #verboseprint(matrix_id, ETPs[matrix_id])
        #casting ETPs to numpy floats
        ETPs=ETPs.astype(np.float64)
        #these test should go in the final version     
        
        '''
        @Gertjan:
            - when running makegrid the ETP_logs overwrite themselves for each error.
                - logfiles should be written outside of this function.
            - also there are (stupid) cases when increasing ETPs does not help. Error message should reflect that 
        '''
        try:
            assert math.isclose(np.sum(ETPs.flatten()), 1, rel_tol=1e-5), f"[-] sum(ETPs): {np.sum(ETPs.flatten())} != 1 (rel_tol=1e-5)"
        except AssertionError:
            with open('log_ETPs.txt', 'w') as file:
                print(f"rates: {rates}", file=file)
                print(f"split_time: {split_time}", file=file)
                for idx, value in np.ndenumerate(ETPs):
                    print(idx, value, file=file)
            sys.exit(f"[-] sum(ETPs): {np.sum(ETPs.flatten())} != 1 (rel_tol=1e-5)")
        try:
            assert np.all(np.logical_and(ETPs>=0, ETPs<=1)), 'Not all ETPs in [0,1].'
        except AssertionError:
            values = ETPs[np.logical_not(np.logical_and(ETPs>=0, ETPs<=1))]
            indices = np.where(np.logical_not(np.logical_and(ETPs>=0, ETPs<=1)))
            indices = list(zip(*indices))
            print("[-] Some ETPs are not in [0,1]. Increase machine precision in the ini file.")
            with open('log_single_ETP.txt', 'w') as file:
                for idx, value in zip(indices, values):
                    print(idx, value, file=file)
            with open('log_ETPs.txt', 'w') as file:
                print(f"rates: {rates}", file=file)
                print(f"split_time: {split_time}", file=file)
                for idx, value in np.ndenumerate(ETPs):
                    print(idx, value, file=file)

        return ETPs
        
    def optimize_parameters(self, data, parameterObj, trackHistory=True, verbose=False, label='', param_combo_name=''):

        fixedParams = parameterObj._get_fixed_params(as_dict=True)
        boundaryNames = self._optimize_get_boundary_names(parameterObj, fixedParams, data)
        starting_points, parameter_combinations_lowest, parameter_combinations_highest = self._optimize_get_boundaries(parameterObj.parameter_combinations, boundaryNames, parameterObj.numPoints)
        trackHistoryPath = [[] for _ in range(parameterObj.numPoints)]
        #specify the objective function to be optimized
        specified_objective_function_list = self._optimize_specify_objective_function(parameterObj.data_type, trackHistoryPath, data, boundaryNames,fixedParams, verbose, parameterObj.threads)
        
        #print to screen
        startdesc = "[+] Optimization starting for specified point."
        if parameterObj.numPoints>1:
            startdesc = f"[+] Optimization starting for specified point and {parameterObj.numPoints-1} random points."
        if verbose:
            print(startdesc)
            #print('iteration \t'+'\t'.join(str(name) for name in boundaryNames)+'\t lnCL')
        #make log file for simulation replicates:
        _init_optimize_log(boundaryNames, f'optimize_log_{label}_{param_combo_name}.tsv', parameterObj._CWD, meta=parameterObj._get_cmd())
        
        #parallelize over starting points or data points
        allResults=[] 
        if parameterObj.gridThreads <= 1:
            with open(os.path.join(parameterObj._CWD, f'optimize_log_{label}_{param_combo_name}.tsv'), 'a') as log_file:
                for idx, (startPos, specified_objective_function) in enumerate(tqdm(zip(itertools.cycle(starting_points), specified_objective_function_list), desc="progress", total=len(specified_objective_function_list),disable=verbose)):
                    single_run_result = run_single_optimize(startPos, parameter_combinations_lowest, parameter_combinations_highest, specified_objective_function, parameterObj.max_eval, parameterObj.xtol_rel, parameterObj.ftol_rel) 
                    allResults.append(single_run_result)
                    _optimize_log(single_run_result, idx, log_file)
                    #allResults.append(run_single_optimize(startPos, parameter_combinations_lowest, parameter_combinations_highest, specified_objective_function, parameterObj.max_eval, parameterObj.xtol_rel, parameterObj.ftol_rel))
            
        else:
            #single_runs need to be specified before passing them to the pool
            specified_run_list = self._optimize_specify_run_list_(parameterObj, specified_objective_function_list, starting_points, parameter_combinations_lowest, parameter_combinations_highest)
            with open(os.path.join(parameterObj._CWD, f'optimize_log_{label}_{param_combo_name}.tsv'), 'a') as log_file:
                with concurrent.futures.ProcessPoolExecutor(max_workers=parameterObj.gridThreads) as outer_pool:
                    for idx, single_run in enumerate(tqdm(outer_pool.map(fp_map, specified_run_list), desc="progress", total=len(specified_run_list), disable=verbose)):
                        #_optimize_log(single_run, identifier, file_name, directory)
                        _optimize_log(single_run, idx, log_file)
                        allResults.append(single_run)

        #process results found in allResults (final step including exit code), and trackHistoryPath (all steps) 
        exitcodeDict = {
                        1: 'optimum found', 
                        2: 'stopvalue reached',
                        3: 'tolerance on lnCL reached',
                        4: 'tolerance on parameter vector reached',
                        5: 'max number of evaluations reached',
                        6: 'max computation time was reached'
                    }
        result = self._optimize_reshape_output(allResults, trackHistoryPath, boundaryNames, exitcodeDict, trackHistory, verbose)
        return result

    def _optimize_get_boundary_names(self, parameterObj, fixedParams, data):
        toBeSynced_pops = [f'Ne_{pop}' for pop in parameterObj.toBeSynced] if parameterObj.toBeSynced!=None else []
        boundaryNames = [k for k in parameterObj.config['parameters'].keys() if not (k=='mu' or k in fixedParams or k in toBeSynced_pops)]
        if len(boundaryNames) == 0:
            print("[-] No boundaries specified.")
            #scale all parameters
            scaled_params = self._scale_parameter_combination(fixedParams, self.reference_pop, self.block_length, parameterObj._MODULE, parameterObj.config['mu']['mu'])
            rate_by_variable = self._get_base_rate_by_variable(scaled_params)
            split_time = self._get_split_time(scaled_params)
            args = rate_by_variable, split_time, parameterObj.threads, True 
            ETPs = self.calculate_ETPs(args)
            CL = calculate_composite_likelihood(ETPs, data)
            print(f"[+] Starting point lnCL={CL}.")
            sys.exit()
        return boundaryNames
    
    def _optimize_get_boundaries(self, parameter_combinations, boundary_names, numPoints):
        # not actually using self, should be standalone function...
        parameter_combinations_lowest = np.array([parameter_combinations[k][0] for k in boundary_names])
        parameter_combinations_highest = np.array([parameter_combinations[k][1] for k in boundary_names])
        #start from midpoint
        pmid = np.mean(np.vstack((parameter_combinations_lowest, parameter_combinations_highest)), axis=0)
        if numPoints > 1:
            #np.random.seed(self.seed)
            starting_points = np.random.uniform(low=parameter_combinations_lowest, high=parameter_combinations_highest, size=(numPoints-1, len(boundary_names)))
            starting_points = np.vstack((pmid,starting_points))
        else:
            starting_points = [pmid,]
    
        return (starting_points, parameter_combinations_lowest, parameter_combinations_highest)
    
    def _optimize_specify_objective_function(self, data_type, trackHistoryPath, data, boundary_names, fixedParams, verbose, threads=1):
        # not actually using self, should be standalone function...
        if data_type=='simulate':
            specified_objective_function_list = [partial(
                    objective_function,
                    paramNames=boundary_names,
                    fixedParams=fixedParams,
                    equationSystemObj=self,
                    data=bsfs,
                    threads=threads,
                    verbose=verbose,
                    path=None) for bsfs in data]
    
        elif data_type=='blocks':
            specified_objective_function_list = [partial(
                        objective_function,
                        paramNames=boundary_names,
                        fixedParams=fixedParams,
                        equationSystemObj=self,
                        data=data,
                        threads=threads,
                        verbose=verbose,
                        path=sublist) for sublist in trackHistoryPath]
        else:
            raise ValueError("Only blocks and sims have been implemented so far.")
        return specified_objective_function_list
    
    def _optimize_specify_run_list_(self, parameterObj, specified_objective_function_list, starting_points, parameter_combinations_lowest, parameter_combinations_highest):
        # not actually using self, should be standalone function...
        specified_run_list=[partial(
                    run_single_optimize,
                    p0=p0,
                    lower=parameter_combinations_lowest,
                    upper=parameter_combinations_highest,
                    specified_objective_function=specified_objective_function,
                    maxeval=parameterObj.max_eval,
                    xtol_rel=parameterObj.xtol_rel,
                    ftol_rel=parameterObj.ftol_rel
                    ) for p0, specified_objective_function in zip(itertools.cycle(starting_points), specified_objective_function_list)]
        return specified_run_list
    
    def _optimize_reshape_output(self, raw_output, trackHistoryPath, boundary_names, exitcodeDict, trackHistory, verbose):
        # not actually using self, should be standalone function...

        #raw_output is list of dicts with {"lnCL":value, "optimum":list, "exitcode":value}
        if verbose:
            print([exitcodeDict.get(resultd['exitcode'],'Not in exitcodeDict.') for resultd in raw_output])
        
        #process trackhistory
        if not trackHistory:
            trackHistoryPath = [list(resultd['optimum'])+[resultd['lnCL'],label, resultd['exitcode']] for label,resultd in enumerate(raw_output)]
            result = [boundary_names+["lnCL", 'iterLabel', 'exitcode'],]+trackHistoryPath 
        else:
            trackHistoryPath = [list(ar)+[idx,] for idx, sublist in enumerate(trackHistoryPath) for ar in sublist]
            result = [boundary_names+["lnCL", 'iterLabel'],]+trackHistoryPath

        return result
    
    def _get_equationObjs(self, sync_ref=None, sync_targets=None, disable_tqdm=False):
        '''
        fix multithreading
        '''
        #print('_get_equationObjs', sync_ref, sync_targets)
        constructors = []
        for constructor_id, event_tuples in tqdm(self.event_tuples_by_idx.items(), total=len(self.event_tuples_by_idx), desc="[%] Building equations", ncols=100, disable=disable_tqdm):
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
            for parameter_batch in tqdm(parameter_batches, desc="[%] Placing mutations", ncols=100, disable=disable_tqdm):
                mutation_tuple, equations = place_mutations(parameter_batch)
                equation_by_mutation_tuple[mutation_tuple] = equations
        else:
            with poolcontext(processes=self.threads) as pool:
                with tqdm(parameter_batches, desc="[%] Placing mutations", ncols=100, disable=disable_tqdm) as pbar:
                    for result_batch in pool.imap_unordered(place_mutations, parameter_batches):
                        mutation_tuple, equations = result_batch
                        equation_by_mutation_tuple[mutation_tuple] = equations
                        pbar.update()
        print("[+] Bundled equations along %s mutation tuples" % (len(equation_by_mutation_tuple)))
        # marginal tuples
        equationObjs = []
        mutation_profiles_offset_1 = get_mutation_profiles(self.k_max_by_mutype, offset=1, by_mutype=True)
        for mutation_profile in tqdm(mutation_profiles_offset_1, desc="[%] Accounting for marginal probabilities", ncols=100, disable=disable_tqdm):
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
        replaceMutypesDict = {sage.all.SR.var(mutype):sage.all.SR.var('theta') for mutype in self.mutypes}
        #substitution of mutype variables by theta: does not work when done in place_mutations
        for equationObj in equationObjs:
            equationObj.equation = equationObj.equation.subs(replaceMutypesDict)
        if sync_ref and sync_targets:
            sync_ref_var = sage.all.SR.symbol(f"C_{sync_ref}")
            for equationObj in equationObjs:
                for sync_target in sync_targets:
                    sync_target_var = sage.all.SR.symbol(f"C_{sync_target}")
                    equationObj.equation = equationObj.equation.subs(sync_target_var==sync_ref_var)
        print("[+] Generated equations for %s mutation tuples" % len(equationObjs))
        return equationObjs
"""