import traceback
import contextlib
import datetime
import itertools
import multiprocessing
import nlopt
import numpy as np
import sys
import functools
import copy
from timeit import default_timer as timer
from tqdm import tqdm
import pandas as pd
import collections

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

def get_nlopt_log_iteration_tuple(
    windows_idx, iteration, block_length, likelihood, scaled_values_by_parameter, unscaled_values_by_parameter, windows_flag, is_anomaly):
    nlopt_log_iteration = [iteration, block_length, likelihood, is_anomaly]
    if windows_flag:
        nlopt_log_iteration = [windows_idx] + nlopt_log_iteration
    for parameter, scaled_value in scaled_values_by_parameter.items():
        if not parameter.startswith("theta"):
            nlopt_log_iteration.append(float(scaled_value)) 
    for parameter, unscaled_value in unscaled_values_by_parameter.items():
        if not parameter.startswith("theta"):
            nlopt_log_iteration.append(float(unscaled_value)) 
    return tuple(nlopt_log_iteration)

def setup_nlopt_log(replicate_idx, config):
    nlopt_log_header = ['iteration', 'block_length', 'lnCL', 'anomaly']
    if config['nlopt_chains'] > 1: # multiple windows
        nlopt_log_header = ["windows_idx"] + nlopt_log_header
    nlopt_log_header += ['%s_scaled' % parameter for parameter in config['gimbleDemographyInstance'].order_of_parameters]
    nlopt_log_header += ['%s_unscaled' % parameter for parameter in config['gimbleDemographyInstance'].order_of_parameters]
    nlopt_log_fn = "%s.%s.%s.log" % (config['nlopt_log_prefix'], replicate_idx, datetime.datetime.now().strftime("%Y%m%d_%H%M%S"))
    print("[+] Trajectories of optimization(s) are written to %r" % nlopt_log_fn)
    with open(nlopt_log_fn, 'w') as nlopt_log_fh:
        nlopt_log_fh.write("%s\n" % ",".join(nlopt_log_header)), nlopt_log_fh.flush()
    return (nlopt_log_fn, nlopt_log_header)

NLOPT_LOG_QUEUE = multiprocessing.Queue()
NLOPT_LOG_ITERATIONS_MANAGER = multiprocessing.Manager()
NLOPT_ANOMALY_COUNTS = collections.Counter()

def nlopt_logger(nlopt_log_fn, nlopt_log_header, nlopt_log_iterations, nlopt_call_count):
    global NLOPT_LOG_QUEUE
    progress_bar_disable_bool = (True if nlopt_call_count == 1 else False) # disable if only one dataset (i.e. blocks)
    pbar = tqdm(desc="[%] Progress", total=nlopt_call_count, position=0, disable=progress_bar_disable_bool, ncols=100)
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
                    break
                else:
                    # iteration array
                    nlopt_log_iterations.append(msg[0])
                    # log file
                    nlopt_log_fh.write("%s\n" % ",".join(map(str, msg[0])))
                    nlopt_log_fh.flush()
                    # progress bar values depend on nlopt_log_header (since they also should be written to log)
                    nlopt_log_iteration_values_by_key = {k: v for k, v in zip(nlopt_log_header, msg[0])}
                    #print('nlopt_log_iteration_values_by_key', nlopt_log_iteration_values_by_key)
                    pbar.write(format_nlopt_log_iteration_values(nlopt_log_iteration_values_by_key))
                    if msg[1]: # this needs to be checked!
                        nlopt_log_fh.write(msg[1])
                        sys.exit('[X] Something went wrong. Check the log file %r' % nlopt_log_fn)

def format_nlopt_result(result_dict):
    boundary_collisions = []
    if result_dict['nlopt_lower_boundary_collision']:
        boundary_collisions.append("%s (lower)" % ", ".join(result_dict['nlopt_lower_boundary_collision']))
    if result_dict['nlopt_upper_boundary_collision']:
        boundary_collisions.append("%s (upper)" % ", ".join(result_dict['nlopt_upper_boundary_collision']))
    status_str = result_dict['nlopt_status']
    boundary_str = "--> [BOUNDARY_COLLISION] : %s" % ("; ".join(boundary_collisions)) if boundary_collisions else ""
    windows_str = "windows_idx=%s" % str(int(result_dict['windows_idx'])) if result_dict['windows_flag'] else ""
    anomaly_str = "ANOMALIES=%s" % result_dict['anomaly_count']
    return "[+] [COMPLETED] %s --------- [%s] [%s] %s" % (windows_str, status_str, anomaly_str, boundary_str)

def format_nlopt_log_iteration_values(nlopt_log_iteration_values_by_key):
    '''determines how screen prints work... and i know... this is a mess... sorry.'''
    value_suffix = '_unscaled' # '_unscaled'
    iteration_str = str(int(nlopt_log_iteration_values_by_key['iteration'])).ljust(4)
    parameter_str = " ".join(["%s=%s" % (k.replace(value_suffix, ''), '{:.5e}'.format(float(v))) for k, v in nlopt_log_iteration_values_by_key.items() if k.endswith(value_suffix)])
    anomaly_str = "[ANOMALY]" if nlopt_log_iteration_values_by_key['anomaly'] else ""
    likelihood_str = "%s %s" % ('{:.5f}'.format(float(nlopt_log_iteration_values_by_key['lnCL'])), anomaly_str)
    windows_str = " windows_idx=%s" % str(int(nlopt_log_iteration_values_by_key['windows_idx'])) if 'windows_idx' in nlopt_log_iteration_values_by_key else ""
    return "[+]%s i=%s -- {%s} -- lnCL=%s" % (windows_str, iteration_str, parameter_str, likelihood_str)

def get_agemo_nlopt_args(dataset, config):
    nlopt_params = {k: v for k, v in config.items() if k.startswith("nlopt")}
    return [(
        nlopt_params, 
        i, 
        functools.partial(
            agemo_likelihood_function,
            dataset=(dataset if config['nlopt_chains'] == 1 else dataset[i]),
            windows_idx=i, 
            config=config, 
            gimbleDemographyInstance=copy.deepcopy(config['gimbleDemographyInstance']), 
            nlopt_iterations=np.zeros(config['nlopt_chains'], dtype=np.uint16), verbose=True),
        (False if config['nlopt_chains'] == 1 else True)
        ) for i in range(config['nlopt_chains'])]

def agemo_calculate_composite_likelihood(ETPs, data):
    ETP_log = np.zeros(ETPs.shape, dtype=np.float64)
    np.log(ETPs, where=ETPs>0, out=ETP_log)
    return np.sum(ETP_log * data)

#def agemo_likelihood_function(nlopt_values, grad, gfEvaluatorObj, dataset, windows_idx, config, nlopt_iterations, verbose):
def agemo_likelihood_function(nlopt_values, grad, gimbleDemographyInstance, dataset, windows_idx, config, nlopt_iterations, verbose, nlopt_anomaly_tol=1e-5, nlopt_anomaly_skip=True):
    global NLOPT_LOG_QUEUE
    #global NLOPT_ANOMALY_COUNTS
    #print('[--------------------------------------------------------] nlopt_iterations', nlopt_iterations)
    nlopt_traceback = None
    nlopt_iterations[windows_idx] += 1
    nlopt_values_by_parameter = {parameter: value for parameter, value in zip(config['nlopt_parameters'], nlopt_values)}
    windows_flag = False if config['nlopt_chains'] == 1 else True # for status bar/log
    scaled_values_by_parameter, unscaled_values_by_parameter = gimbleDemographyInstance.scale_parameters(nlopt_values_by_parameter)
    theta_branch, var, time, fallback_flag = gimbleDemographyInstance.get_agemo_values(scaled_values_by_parameter, fallback=True)
    evaluator = EVALUATOR if not fallback_flag else FALLBACK_EVALUATOR
    ETPs = evaluator.evaluate(theta_branch, var, time=time)
    # df = pd.DataFrame(bsfs_to_2d(ETPs))
    # problematic = df.loc[(df[3] >= 1) & (df[4] >= 1)]
    # print(problematic.to_markdown())
    #print('[+] np.sum(ETPs)=%s me=%s ; fallback=%s; EVALUATOR.evaluate(%s, np.array(%s), time=%s)' % (np.sum(ETPs), str(unscaled_values_by_parameter['me']), fallback_flag, str(theta_branch), str([v for v in var]), str(time)))
    #print('[+] nlopt_values_by_parameter', nlopt_values_by_parameter)
    #print('[+] scaled_values_by_parameter', scaled_values_by_parameter)
    # print("[+] Fixing ETPs ...")
    # ETPs[(ETPs[:,2] > 0) & (ETPs[:,3] > 0)] = 0
    # print('\nnp.sum(ETPs)=%s' % np.sum(ETPs))

    # Determine anomaly
    is_anomaly = (not np.isclose(np.sum(ETPs), 1, rtol=nlopt_anomaly_tol))
    if is_anomaly:
        NLOPT_ANOMALY_COUNTS[windows_idx] += 1
    try:
        likelihood = agemo_calculate_composite_likelihood(ETPs, dataset) if not (is_anomaly and nlopt_anomaly_skip) else -np.inf
    except Exception as exception:
        nlopt_log_iteration_tuple = get_nlopt_log_iteration_tuple(windows_idx, nlopt_iterations[windows_idx], gimbleDemographyInstance.block_length, likelihood, scaled_values_by_parameter, unscaled_values_by_parameter, windows_flag, is_anomaly)
        nlopt_traceback = traceback.format_exc(exception)
        NLOPT_LOG_QUEUE.put((nlopt_log_iteration_tuple, nlopt_traceback))
    nlopt_log_iteration_tuple = get_nlopt_log_iteration_tuple(windows_idx, nlopt_iterations[windows_idx], gimbleDemographyInstance.block_length, likelihood, scaled_values_by_parameter, unscaled_values_by_parameter, windows_flag, is_anomaly)
    NLOPT_LOG_QUEUE.put((nlopt_log_iteration_tuple, nlopt_traceback))
    return likelihood

def agemo_optimize(evaluator_agemo, replicate_idx, data, config, fallback_evaluator=None):
    global EVALUATOR
    global FALLBACK_EVALUATOR
    #np.set_printoptions(precision=19, suppress = True)
    EVALUATOR = evaluator_agemo
    FALLBACK_EVALUATOR = fallback_evaluator
    nlopt_runs = get_agemo_nlopt_args(data, config)
    #print('nlopt_runs', nlopt_runs)
    nlopt_results = []
    # Setup nlopt_logger/tqdm
    #print('### replicate_idx', replicate_idx)
    #print('### nlopt_runs', nlopt_runs)
    nlopt_log_fn, nlopt_log_header = setup_nlopt_log(replicate_idx, config)
    #
    nlopt_log_iteration_tuples = NLOPT_LOG_ITERATIONS_MANAGER.list()
    nlopt_log_process = multiprocessing.Process(target=nlopt_logger, args=(nlopt_log_fn, nlopt_log_header, nlopt_log_iteration_tuples, len(nlopt_runs)))
    nlopt_log_process.start()
    if config['processes'] <= 1:
        for nlopt_run in nlopt_runs:
            nlopt_results.append(agemo_nlopt_call(nlopt_run))
    else:
        with poolcontext(processes=config['processes']) as pool:
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
    global NLOPT_ANOMALY_COUNT
    NLOPT_ALGORITHMS = {
        'neldermead' : nlopt.LN_NELDERMEAD,
        'sbplx': nlopt.LN_SBPLX,
        'CRS2': nlopt.GN_CRS2_LM,}

    nlopt_params, windows_idx, agemo_likelihood_function, windows_flag = args
    #print(nlopt_params)
    num_optimization_parameters = len(nlopt_params['nlopt_start_point'])
    nlopt_algorithm = NLOPT_ALGORITHMS[nlopt_params.get('nlopt_algorithm', 'CRS2')]
    nlopt.srand(nlopt_params['nlopt_seed'])
    opt = nlopt.opt(nlopt_algorithm, num_optimization_parameters)
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
    # optimum[0] = nlopt_params['nlopt_lower_bound'][0] # BOUNDARY COLLISIONS
    # optimum[1] = nlopt_params['nlopt_upper_bound'][1] # BOUNDARY COLLISIONS
    nlopt_result = {
        'windows_idx': windows_idx,
        'nlopt_optimum': opt.last_optimum_value(),
        'nlopt_evals': opt.get_numevals(),
        'nlopt_values': optimum,
        'nlopt_status': NLOPT_EXIT_CODE[opt.last_optimize_result()],
        'nlopt_lower_boundary_collision': [nlopt_params['nlopt_parameters'][i] for i, value in enumerate(optimum) if value == nlopt_params['nlopt_lower_bound'][i]],
        'nlopt_upper_boundary_collision': [nlopt_params['nlopt_parameters'][i] for i, value in enumerate(optimum) if value == nlopt_params['nlopt_upper_bound'][i]],
        'windows_flag': windows_flag,
        'anomaly_count': NLOPT_ANOMALY_COUNTS[windows_idx]}
    #print("nlopt_result", nlopt_result)
    NLOPT_LOG_QUEUE.put(nlopt_result)
    return nlopt_result

def get_agemo_optimize_result(config, nlopt_results, nlopt_log_header, nlopt_log_iteration_tuples):
    # optimize results are done by dataset, i.e. blocks or windows or replicates
    windows_flag = False if config['nlopt_chains'] == 1 else True
    optimize_result = {
        'nlopt_log_iteration_header': nlopt_log_header,
        'nlopt_log_iteration_array' : np.asarray(nlopt_log_iteration_tuples),
        'dataset_count': len(nlopt_results),
        'nlopt_values_by_windows_idx': {},
        'nlopt_optimum_by_windows_idx': {},
        'nlopt_evals_by_windows_idx': {},
        'nlopt_status_by_windows_idx': {},
        'nlopt_anomalies_by_windows_idx': {},
        'nlopt_lower_boundary_collision_by_windows_idx': {},
        'nlopt_upper_boundary_collision_by_windows_idx': {},
        'windows_flag': windows_flag
        }
    for nlopt_result in nlopt_results:
        windows_idx = nlopt_result['windows_idx']
        nlopt_values = {parameter: value for parameter, value in zip(config['nlopt_parameters'], nlopt_result['nlopt_values'])}
        fixed_values = {parameter: getattr(config['gimbleDemographyInstance'], parameter) for parameter in config['gimbleDemographyInstance'].fixed_parameters}
        optimize_result['nlopt_values_by_windows_idx'][windows_idx] = {**nlopt_values, **fixed_values}
        optimize_result['nlopt_optimum_by_windows_idx'][windows_idx] = nlopt_result['nlopt_optimum']
        optimize_result['nlopt_evals_by_windows_idx'][windows_idx] = nlopt_result['nlopt_evals']
        optimize_result['nlopt_status_by_windows_idx'][windows_idx] = nlopt_result['nlopt_status']
        optimize_result['nlopt_lower_boundary_collision_by_windows_idx'][windows_idx] = nlopt_result['nlopt_lower_boundary_collision']
        optimize_result['nlopt_upper_boundary_collision_by_windows_idx'][windows_idx] = nlopt_result['nlopt_upper_boundary_collision']
        optimize_result['nlopt_anomalies_by_windows_idx'][windows_idx] = nlopt_result['anomaly_count']
    return optimize_result
