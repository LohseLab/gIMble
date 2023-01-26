import traceback
import contextlib
import datetime
import itertools
import multiprocessing
import nlopt
import numpy as np
import pandas as pd
import sys, os
import functools
import copy
from functools import partial
from functools import partialmethod
from timeit import default_timer as timer
from tqdm import tqdm
import lib.gimble

def bsfs_to_2d(bsfs):
    '''needed for testing that tally-arrays and bsfs are identical'''
    if not np.any(bsfs):
        return None
    non_zero_idxs = np.nonzero(bsfs)
    if bsfs.ndim == 4: # blocks
        return np.concatenate([bsfs[non_zero_idxs].reshape(non_zero_idxs[0].shape[0], 1), np.array(non_zero_idxs, dtype=np.uint64).T], axis=1)
    elif bsfs.ndim == 5: # windows
        non_zero_idxs_array = np.array(non_zero_idxs, dtype=np.uint64).T
        first = non_zero_idxs_array[:,0].reshape(non_zero_idxs[0].shape[0], 1)
        second = bsfs[non_zero_idxs].reshape(non_zero_idxs[0].shape[0], 1)
        third = non_zero_idxs_array[:,1:]
        return np.concatenate([first, second, third], axis=1)
    else:
        raise ValueError('bsfs_to_2d: bsfs.ndim must be 4 (blocks) or 5 (windows)')

NLOPT_EXIT_CODE = {
     1: 'NLOPT_SUCCESS',
     2: 'NLOPT_STOPVAL_REACHED',
     3: 'NLOPT_FTOL_REACHED',
     4: 'NLOPT_XTOL_REACHED',
     5: 'NLOPT_MAXEVAL_REACHED',
     6: 'NLOPT_MAXTIME_REACHED',
     -1: 'NLOPT_FAILURE',
     -2: 'NLOPT_INVALID_ARGS',
     -3: 'NLOPT_OUT_OF_MEMORY',
     -4: 'NLOPT_ROUNDOFF_LIMITED',
     -5: 'NLOPT_FORCED_STOP'
}

@contextlib.contextmanager
def poolcontext(*args, **kwargs):
    pool = multiprocessing.Pool(*args, **kwargs)
    yield pool
    pool.terminate()     

def get_nlopt_log_iteration_tuple(dataset_idx, iteration, block_length, likelihood, scaled_values_by_parameter, unscaled_values_by_parameter):
    nlopt_log_iteration = [dataset_idx, iteration, block_length, likelihood]
    for parameter, scaled_value in scaled_values_by_parameter.items():
        if not parameter.startswith("theta"):
            nlopt_log_iteration.append(float(scaled_value)) 
    for parameter, unscaled_value in unscaled_values_by_parameter.items():
        if not parameter.startswith("theta"):
            nlopt_log_iteration.append(float(unscaled_value)) 
    return tuple(nlopt_log_iteration)

def get_nlopt_log_fn(data_idx, config):
    # optimize_meta = lib.gimble.config_to_meta(config, 'optimize') # could be used to populate header further with data
    nlopt_log_header = ["dataset_idx", 'iteration', 'block_length', 'likelihood']
    nlopt_log_header += ['%s_scaled' % parameter for parameter in config['gimbleDemographyInstance'].order_of_parameters]
    nlopt_log_header += ['%s_unscaled' % parameter for parameter in config['gimbleDemographyInstance'].order_of_parameters]
    nlopt_log_fn = "gimble.optimize.%s.%s.%s.log" % (config['optimize_label'], data_idx, datetime.datetime.now().strftime("%Y%m%d_%H%M%S"))
    print("[+] Trajectories of optimization(s) are written to %r" % nlopt_log_fn)
    with open(nlopt_log_fn, 'w') as nlopt_log_fh:
        nlopt_log_fh.write("%s\n" % ",".join(nlopt_log_header)), nlopt_log_fh.flush()
    return (nlopt_log_fn, nlopt_log_header)

NLOPT_LOG_QUEUE = multiprocessing.Queue()
NLOPT_LOG_ITERATIONS_MANAGER = multiprocessing.Manager()

def nlopt_logger(nlopt_log_fn, nlopt_log_header, nlopt_log_iterations, nlopt_call_count):
    global NLOPT_LOG_QUEUE
    progress_bar_disable_bool = (True if nlopt_call_count == 1 else False) # disable if only one dataset (i.e. blocks)
    pbar = tqdm(desc="[%] Progress", total=nlopt_call_count, position=0, disable=progress_bar_disable_bool)
    with open(nlopt_log_fn, 'a') as nlopt_log_fh:
        while 1:
            msg = NLOPT_LOG_QUEUE.get()
            if isinstance(msg, dict):
                # received from nlopt_call() upon finished optimization
                pbar.write(format_nlopt_result(msg)) 
                pbar.update()
            if isinstance(msg, tuple):
                # received from likelihood_function() after each iteration (with traceback if failed)
                if msg[0] is None:
                    pbar.update()
                    break
                else:
                    # iteration array
                    nlopt_log_iterations.append(msg[0])
                    # log file
                    nlopt_log_fh.write("%s\n" % ",".join(map(str, msg[0])))
                    nlopt_log_fh.flush()
                    # progress bar
                    nlopt_log_iteration_values_by_key = {k: v for k, v in zip(nlopt_log_header, msg[0])}
                    #print('nlopt_log_iteration_values_by_key', nlopt_log_iteration_values_by_key)
                    pbar.write(format_nlopt_log_iteration_values(nlopt_log_iteration_values_by_key))
                    if msg[1]: # this needs to be checked!
                        nlopt_log_fh.write(msg[1])
                        sys.exit('[X] Something went wrong. Check the log file %r' % nlopt_log_fn)

def format_nlopt_result(result_dict):
    return "[+] data_idx=%s --------- [%s]" % (
        int(result_dict['dataset_idx']),
        result_dict['nlopt_status'])

def format_nlopt_log_iteration_values(nlopt_log_iteration_values_by_key):
    '''determines how log and screen prints work'''
    value_suffix = '_unscaled' # '_unscaled'
    return "[+] data_idx=%s i=%s -- {%s} -- L=%s" % (
        int(nlopt_log_iteration_values_by_key['dataset_idx']),
        str(int(nlopt_log_iteration_values_by_key['iteration'])).ljust(4),
        " ".join(["%s=%s" % (k.replace(value_suffix, ''), '{:.5e}'.format(float(v))) for k, v in nlopt_log_iteration_values_by_key.items() if k.endswith(value_suffix)]),
        '{:.5f}'.format(float(nlopt_log_iteration_values_by_key['likelihood'])))

def get_agemo_nlopt_args(dataset, config):
    nlopt_params = {k: v for k, v in config.items() if k.startswith("nlopt")}
    return [(
        nlopt_params, 
        i, 
        functools.partial(
            agemo_likelihood_function,
            dataset=(dataset if config['nlopt_chains'] == 1 else dataset[i]),
            dataset_idx=i, 
            config=config, 
            gimbleDemographyInstance=copy.deepcopy(config['gimbleDemographyInstance']), 
            nlopt_iterations=np.zeros(config['nlopt_chains'], dtype=np.uint16), verbose=True)
        ) for i in range(config['nlopt_chains'])]

def agemo_calculate_composite_likelihood(ETPs, data):
    assert ETPs.shape == data.shape, "[X] Incompatible shapes in calculate_composite_likelihood(): ETPs (%s) vs data (%s)" % (ETPs.shape, data.shape)
    ETP_log = np.zeros(ETPs.shape, dtype=np.float64)
    np.log(ETPs, where=ETPs>0, out=ETP_log)
    return np.sum(ETP_log * data)

#def agemo_likelihood_function(nlopt_values, grad, gfEvaluatorObj, dataset, dataset_idx, config, nlopt_iterations, verbose):
def agemo_likelihood_function(nlopt_values, grad, gimbleDemographyInstance, dataset, dataset_idx, config, nlopt_iterations, verbose):
    global NLOPT_LOG_QUEUE
    #print('[--------------------------------------------------------] nlopt_iterations', nlopt_iterations)
    nlopt_traceback = None
    nlopt_iterations[dataset_idx] += 1
    nlopt_values_by_parameter = {parameter: value for parameter, value in zip(config['nlopt_parameters'], nlopt_values)}
    #print('[+] nlopt_values_by_parameter', nlopt_values_by_parameter)
    scaled_values_by_parameter, unscaled_values_by_parameter = gimbleDemographyInstance.scale_parameters(nlopt_values_by_parameter)
    #print('[+] scaled_values_by_parameter', scaled_values_by_parameter)
    #print('[+] unscaled_values_by_parameter', unscaled_values_by_parameter)

    theta_branch, var, time, fallback_flag = gimbleDemographyInstance.get_agemo_values(scaled_values_by_parameter, fallback=True)
    #ETPs = EVALUATOR.evaluate(theta_branch, np.positive(var), time=time)
    #print('[+] scaled_values_by_parameter', scaled_values_by_parameter)
    #print('[+] unscaled_values_by_parameter', unscaled_values_by_parameter)
    #print('[+] me=%s ; fallback=%s; EVALUATOR.evaluate(%s, np.array(%s), time=%s)' % (str(unscaled_values_by_parameter['me']), fallback_flag, str(theta_branch), str([v for v in var]), str(time)), end="")
    evaluator = EVALUATOR if not fallback_flag else FALLBACK_EVALUATOR
    #_var = str([str(x) for x in var])
    ETPs = evaluator.evaluate(theta_branch, var, time=time)
    
    #if ETPs[(ETPs[:,2] > 0) & (ETPs[:,3] > 0)] > 0:
    #    df = pd.DataFrame(bsfs_to_2d(ETPs))
    #    problematic = df.loc[(df[3] >= 1) & (df[4] >= 1)]
    #    print(problematic.to_markdown())
    ETP_sum = np.sum(ETPs)
    if not np.isclose(ETP_sum, 1, rtol=1e-05):
        print('np.sum(ETPs)=%s fallback=%s scaled_values=%s' % (np.sum(ETPs), fallback_flag, scaled_values_by_parameter))
    

    likelihood = agemo_calculate_composite_likelihood(ETPs, dataset)
    
    # scaled_values_by_symbol, scaled_values_by_parameter, unscaled_values_by_parameter, block_length = scale_nlopt_values(nlopt_values, config)
    # try:
    #     gimbleDemographyInstance.set_parameters_from_array(nlopt_values)
    #     theta_branch, var, time = gimbleDemographyInstance.get_agemo_values()
    #     ETPs = EVALUATOR.evaluate(theta_branch, var, time=time)
    #      = agemo_calculate_composite_likelihood(ETPs, dataset)
    # except Exception as exception:
    #     nlopt_log_iteration_tuple = get_nlopt_log_iteration_tuple(dataset_idx, nlopt_iterations[dataset_idx], gimbleDemographyInstance.block_length, 'N/A', scaled_values_by_parameter, unscaled_values_by_parameter)
    #     nlopt_traceback = traceback.format_exc(exception)
    #     NLOPT_LOG_QUEUE.put((nlopt_log_iteration_tuple, nlopt_traceback))
    

    
    nlopt_log_iteration_tuple = get_nlopt_log_iteration_tuple(dataset_idx, nlopt_iterations[dataset_idx], gimbleDemographyInstance.block_length, likelihood, scaled_values_by_parameter, unscaled_values_by_parameter)
    NLOPT_LOG_QUEUE.put((nlopt_log_iteration_tuple, nlopt_traceback))
    
    #if verbose:
    #    elapsed = lib.gimble.format_time(timer() - start_time)
    #    #process_idx = multiprocessing.current_process()._identity
    #    #print_nlopt_line(dataset_idx, nlopt_iterations[dataset_idx], likelihood, unscaled_values_by_parameter, elapsed, process_idx)
    #    print_nlopt_line(dataset_idx, nlopt_iterations[dataset_idx], likelihood, unscaled_values_by_parameter, elapsed)
    return likelihood

def agemo_optimize(evaluator_agemo, data_idx, data, config, fallback_evaluator=None):
    global EVALUATOR
    global FALLBACK_EVALUATOR
    np.set_printoptions(precision=19, suppress = True)
    EVALUATOR = evaluator_agemo
    FALLBACK_EVALUATOR = fallback_evaluator
    nlopt_runs = get_agemo_nlopt_args(data, config)
    nlopt_results = []
    # Setup nlopt_logger/tqdm
    nlopt_log_fn, nlopt_log_header = get_nlopt_log_fn(data_idx, config)
    nlopt_log_iteration_tuples = NLOPT_LOG_ITERATIONS_MANAGER.list()
    nlopt_log_process = multiprocessing.Process(target=nlopt_logger, args=(nlopt_log_fn, nlopt_log_header, nlopt_log_iteration_tuples, len(nlopt_runs)))
    nlopt_log_process.start()
    if config['nlopt_processes'] <= 1:
        for nlopt_run in nlopt_runs:
            nlopt_results.append(agemo_nlopt_call(nlopt_run))
    else:
        with poolcontext(processes=config['nlopt_processes']) as pool:
            for nlopt_result in pool.imap_unordered(agemo_nlopt_call, nlopt_runs):
                nlopt_results.append(nlopt_result)
    # clean up logger
    NLOPT_LOG_QUEUE.put((None, None))
    nlopt_log_process.join()
    # sanitize results
    optimize_result = get_agemo_optimize_result(config, nlopt_results, nlopt_log_header, nlopt_log_iteration_tuples)
    return optimize_result

def agemo_nlopt_call(args):
    global NLOPT_LOG_QUEUE
    NLOPT_ALGORITHMS = {
        'neldermead' : nlopt.LN_NELDERMEAD,
        'sbplx': nlopt.LN_SBPLX
    }

    nlopt_params, dataset_idx, agemo_likelihood_function = args
    num_optimization_parameters = len(nlopt_params['nlopt_start_point'])
    opt = nlopt.opt(NLOPT_ALGORITHMS['sbplx'], num_optimization_parameters)
    opt.set_lower_bounds(nlopt_params['nlopt_lower_bound'])
    opt.set_upper_bounds(nlopt_params['nlopt_upper_bound'])
    opt.set_max_objective(agemo_likelihood_function)
    opt.set_xtol_rel(nlopt_params['nlopt_xtol_rel'])
    # TODO: xtol_weights needs to be revisited ... 
    if nlopt_params['nlopt_xtol_rel'] > 0:
        # assigning weights to address size difference between params
        xtol_weights = 1/(len(nlopt_params['nlopt_start_point']) * nlopt_params['nlopt_start_point'])
        opt.set_x_weights(xtol_weights)
    opt.set_ftol_rel(nlopt_params['nlopt_ftol_rel'])
    opt.set_maxeval(nlopt_params['nlopt_maxeval'])
    optimum = opt.optimize(nlopt_params['nlopt_start_point'])
    nlopt_result = {
        'dataset_idx': dataset_idx,
        'nlopt_optimum': opt.last_optimum_value(),
        'nlopt_values': optimum,
        'nlopt_status': NLOPT_EXIT_CODE[opt.last_optimize_result()]}
    NLOPT_LOG_QUEUE.put(nlopt_result)
    return nlopt_result

def get_agemo_optimize_result(config, nlopt_results, nlopt_log_header, nlopt_log_iteration_tuples):
    # optimize results are done by dataset, i.e. blocks or windows or replicates
    optimize_result = {
        'nlopt_log_iteration_header': nlopt_log_header,
        'nlopt_log_iteration_array' : np.asarray(nlopt_log_iteration_tuples),
        'dataset_count': len(nlopt_results),
        'nlopt_values_by_dataset_idx': {},
        'nlopt_optimum_by_dataset_idx': {},
        'nlopt_status_by_dataset_idx': {}}
    for nlopt_result in nlopt_results:
        dataset_idx = nlopt_result['dataset_idx']
        nlopt_values = {parameter: value for parameter, value in zip(config['nlopt_parameters'], nlopt_result['nlopt_values'])}
        fixed_values = {parameter: getattr(config['gimbleDemographyInstance'], parameter) for parameter in config['gimbleDemographyInstance'].fixed_parameters}
        optimize_result['nlopt_values_by_dataset_idx'][dataset_idx] = {**nlopt_values, **fixed_values}
        optimize_result['nlopt_optimum_by_dataset_idx'][dataset_idx] = nlopt_result['nlopt_optimum']
        optimize_result['nlopt_status_by_dataset_idx'][dataset_idx] = nlopt_result['nlopt_status']
    return optimize_result

def fp_map(f, *args):
    return f(*args)

def calculate_probabilities(evaluator, modelObjs):
    list_of_parameters_arrays = get_parameter_arrays(modelObjs)
    print("[+] Calculating probabilities for %s points in parameter space" % len(scaled_parameter_combinations))
    if processes==1:
        for parameter_combination in tqdm(scaled_parameter_combinations, desc="[%]", ncols=100, disable=True):
            result = gfEvaluatorObj.evaluate_gf(parameter_combination, parameter_combination[sage.all.SR.var('theta')]) 
            all_ETPs.append(result)
                
    else:
        args = ((param_combo, param_combo[sage.all.SR.var('theta')]) for param_combo in scaled_parameter_combinations)
        with multiprocessing.Pool(processes=processes) as pool:
            for ETP in pool.starmap(gfEvaluatorObj.evaluate_gf, tqdm(args, total=len(scaled_parameter_combinations), ncols=100, desc="[%]")):
                all_ETPs.append(ETP)
    return np.array(all_ETPs, dtype=np.float64)
