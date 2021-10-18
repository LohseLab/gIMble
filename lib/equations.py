import itertools
import sys, os
import sage.all
import sage.parallel.multiprocessing_sage
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
from functools import partial
from functools import partialmethod
import nlopt
import concurrent.futures
#import logging
import datetime
from timeit import default_timer as timer

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
        """function redundant """
        #self.reference_pop = parameterObj.config['populations']['reference_pop']
        #self.block_length = sage.all.Rational(parameterObj.config['mu']['blocklength'])
        result = []
        """
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
        """
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