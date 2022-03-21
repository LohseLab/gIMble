import itertools
from tqdm import tqdm
from lib.functions import plot_mutuple_barchart
import allel
import demes
import ast
import math
import numpy as np
import dask
from dask.diagnostics import ProgressBar as daProgressBar
import pandas as pd
import shutil
import zarr
import os
import string
#import loguru
import collections
import sys
import warnings
import pathlib
import configparser
import matplotlib.pyplot as plt
from matplotlib.path import Path
import matplotlib.patches as patches
import cerberus
import lib.simulate
#import dask
import hashlib 
import fractions
import copy
import lib.math
import tabulate
from timeit import default_timer as timer
# np.set_printoptions(threshold=sys.maxsize)

from lib.GeneratingFunction.gf import togimble
'''
[Rules for better living]

- gimbleStore.data.attrs (meta): ZARR JSON encoder does not like numpy/pandas dtypes, have to be converted to python dtypes

- Coordinate systems:
    {GIMBLE}                    : 0-based       | {GIMBLE}  
    -----------------------------------------------------------
    {GENOME_FILE} : LENGTH      : 0-based   =>  | [0, l)
    {BED_FILE}    : START:END   : 0-based   =>  | [START, END)
    {VCF_FILE}    : POS         : 1-based   =>  | [POS-1)
'''

'''

[FASTA parsing]
- https://pypi.org/project/pyfasta/

    
[QC plots]
    - variants 
        - plot barcharts of HOMREF/HOMALT/HET/MISS/MULTI as proportion of total records
    - intervals
        - plot length
        - plot distance
    - mutuples
        - mutuple pcp
'''

PURPLE = '#4F3D63'

DFRM = '[──'
DFRT = '├──'
DPPM = '    '
DFRL = '└──'
DPPL = '│   '

SIGNS = {
        'T': '%s%s%s' % (u"\u251C", u"\u2500", u"\u2500"),  # '├──'
        'S': '%s%s%s' % (" ", " ", " "),                    # '    '
        'F': '%s%s%s' % (u"\u2514", u"\u2500", u"\u2500"),  # '└──'
        'W': '%s%s%s' % (u"\u250C", u"\u2500", u"\u2500"),  # '└──'
        'P': '%s%s%s' % (u"\u2502", " ", " "),              # '│   '
        'B': '%s%s%s' % (u"\u2500", u"\u2500",  u"\u2500")} # '───'

SPACING = 16
MUTYPES = ['m_1', 'm_2', 'm_3', 'm_4']

GRIDSEARCH_DTYPE=np.float32 # -3.4028235e+38 ... 3.4028235e+38
# GRIDSEARCH_DTYPE=np.float64 # -1.7976931348623157e+308 ... 1.7976931348623157e+308

class ReportObj(object):
    '''Report class for making reports'''

    def __init__(self, width=80):
        self.width = width
        self.out = []
    
    def add_line(self, prefix="[+]", branch="", left='', center='', right='', fill=' '):
        '''
        Writes line to Report Object
        Lines is composed of {prefix} {branch} {left} {fill} {center} {fill} {right}
        SIGNS = {'T': '├──', 'F': '└──', 'S': '    ', 'P': '│   ', 'B': '───'}
        '''
        start = "%s %s" % (prefix, "".join([SIGNS[b] for b in branch])) if branch else prefix
        w = self.width - len(left) - len(right) - len(start) 
        self.out.append("%s %s %s %s" % (start, left, (" %s " % center).center(w, fill) if center else fill * w, right))

    def __add__(self, other):
        if isinstance(other, ReportObj):
            if self.width == other.width:
                result = ReportObj(width=self.width)
                result.out = self.out + other.out
                return result
            else:
                raise ValueError("addition only supported for ReportObj of equal width")    
        else:
            raise ValueError("addition only supported for ReportObj's")

    def __radd__(self, other):
        return self.__add__(other)

    def __repr__(self):
        return "\n".join(self.out)

def get_ini_model_events(model, pop_ids):
    events = []
    if model in set(['DIV', 'IM_BA', 'IM_AB']):
        events += ['J_A_B']
    if model in set(['MIG_AB', 'IM_AB']):
        events += ['M_A_B']
    if model in set(['MIG_BA', 'IM_BA']):
        events += ['M_B_A']
    events = ['C_%s' % pop_id for pop_id in pop_ids] + events
    return events

def get_ini_model_params(model):
    pop_ids = ['A', 'B']
    if model in set(['DIV', 'IM_BA', 'IM_AB']):
        pop_ids += ['A_B']
    pop_ids_sync = sorted(set([",".join(combination) for combination 
        in itertools.combinations(pop_ids, 2)] + [','.join(pop_ids),]))
    events = get_ini_model_events(model, pop_ids)
    return pop_ids, pop_ids_sync, events

def get_ini_parameter_list(task, events):
    # [parser.py]
    l = ['# Model parameters and their values/ranges:']
    if task == 'simulate':
        l.append('# A) Fixed : value')
        l.append('# B) Range : min, max, steps, lin|log')
    elif task == 'optimize':
        l.append('# A) Fixed : value')
        l.append('# B) Range : min, max')
    elif task == 'makegrid':
        l.append('# A) Fixed : value')
        l.append('# B) Range : min, max, steps, lin|log')
    else:
        raise ValueError("%s is not a supported task." % task)  
    for event in events:
        event_type, event_name = event[:1], event[2:]
        if event_type == 'C':
            l.append('# Effective population size of %s' % event_name)
            l.append('Ne_%s' % event_name)
        if event_type == 'M':
            source, sink = event_name.split("_")
            l.append('# Migration rate (migrants/generation) from %s to %s (backwards in time)' % (source, sink))
            l.append('me')
        if event_type == 'J':
            l.append('# Split time (in generations)')
            l.append('T')
    return l

def make_ini_configparser(version, task, model, label):
    # [parser.py]
    # make configparser
    config = configparser.ConfigParser(allow_no_value=True)
    config.optionxform = str # otherwise keys are lowercase
    # get variables/strings
    random_seed, precision = '19', '165'
    pop_ids, pop_ids_sync, events = get_ini_model_params(model)
    string_pop_ids = " | ".join(pop_ids)
    string_pop_ids_sync = " | ".join(pop_ids_sync)
    parameter_list = get_ini_parameter_list(task, events)
    # gimble
    config.add_section('gimble')
    config.set('gimble', 'version', version)
    config.set('gimble', 'label', label)
    config.set('gimble', 'task', task)
    config.set('gimble', 'model', model)
    config.set('gimble', 'random_seed', random_seed)
    config.set('gimble', 'precision', precision)
    # populations
    # do only makegrid and optimize need populations?
    if task == 'makegrid' or task == 'optimize':
        config.add_section('populations')
        config.set('populations', 'pop_ids', string_pop_ids)
        if task != 'simulate': 
            config.set('populations', '# Link model to data in GimbleStore (see output of gimble info)')
            config.set('populations', 'A', "")
            config.set('populations', 'B', "")
            config.set('populations', '# Pick a reference population: %s' % string_pop_ids)
            config.set('populations', 'reference_pop_id', "")
        config.set('populations', '# Choose to simplify model by assuming equality of Ne\'s (optional): %s' % string_pop_ids_sync)
        config.set('populations', 'sync_pop_ids', "")
    # simulate
    if task == 'simulate' or task == 'makegrid':
        config.add_section('k_max')
        config.set('k_max', '# k_max sets limits for mutation-type cardinalities.')
        config.set('k_max', '# Mutation-types beyond these are calculated as marginals')
        config.set('k_max', 'm_1', '2    # hetB')
        config.set('k_max', 'm_2', '2    # hetA')
        config.set('k_max', 'm_3', '2    # hetAB')
        config.set('k_max', 'm_4', '2    # fixed')
    if task == 'simulate':
        config.add_section('simulate')
        config.set('simulate', 'pop_ids', string_pop_ids)
        config.set('simulate', '# Pick a reference population: %s' % string_pop_ids)
        config.set('simulate', 'reference_pop_id', "")
        config.set('simulate', '# Ploidy of organism')
        config.set('simulate', 'ploidy', '2')
        config.set('simulate', '# Blocks')
        config.set('simulate', 'blocks', "")
        config.set('simulate', 'block_length', "")
        config.set('simulate', 'num_linked_blocks', "")
        config.set('simulate', '# Number of replicates')
        config.set('simulate', 'replicates', "")
        config.set('simulate', '# Number of samples per population')
        config.set('simulate', 'sample_size_A', '1')
        config.set('simulate', 'sample_size_B', '1')
        config.set('simulate', '# Mutations at discrete sites (True) or infinite sites model (False)')
        config.set('simulate', 'discrete_genome', 'True')
        config.set('simulate', '# Set recombination rate (optional)')
        config.set('simulate', 'recombination_rate', "")
        config.set('simulate', '# Set path to recombination map (optional)')
        config.set('simulate', 'recombination_map', "")
        config.set('simulate', '# Number of bins, cutoff, scale (required ONLY IF recombination map set)')
        config.set('simulate', 'number_bins', "")
        config.set('simulate', 'cutoff', "")
        config.set('simulate', 'scale', "")
        config.add_section('gridbased')
        config.set('gridbased', 'grid_label', "")
        config.set('gridbased', '# Parameter to fix to global optimum when performing window-wise parametric bootstrap (optional).')
        config.set('gridbased', 'fixed_parameter', "")
    # mu
    config.add_section('mu') 
    config.set('mu', '# mutation rate (in mutations/site/generation, required)')
    # [To Do] figure out how to deal with 'no mutation-rate' ... no-scaling
    config.set('mu', 'mu', "")
    if not task == 'simulate' or not task == 'optimize':
        config.set('mu', '# block_length')
        config.set('mu', '# must be identical to block_length of the data one wants to analyse with grid')
        config.set('mu', 'block_length', "")

    # parameters
    config.add_section('parameters')
    for param in parameter_list:
        if param.startswith("#"):
            config.set('parameters', param)
        else:
            config.set('parameters', param, "")
    return config

def get_model_name(model_file):
    return model_file.rstrip('.tsv').split('/')[-1]

def write_config(version, model, task, label):
    config = make_ini_configparser(version, task, model, label)
    outfile = "gimble.%s.%s.%s.config.ini" % (task, model, label)
    with open(outfile, 'w') as fh:
        config.write(fh)
    print("[+] Wrote INI file %r. Please fill in values before starting %r." % (str(outfile), task))

def get_validator_error_string(validator_errors):
    # parameterObj file ...
    out = []
    for section, errors in validator_errors.items():
        out.append("[X] In section %r" % section)
        for error_dict in errors:
            for parameter, values in error_dict.items():
                out.append("[X] \t %s \t %s" % (parameter, " ".join(values)))
    return "\n".join(out)

def DOL_to_LOD(DOL):
    """
    converts dict of lists to list of dicts
    """
    reshape_DOL = list(zip(*DOL.values()))
    return [{k:v for k,v in zip(DOL.keys(), sublist)} for sublist in reshape_DOL]  

def LOD_to_DOL(LOD):
    """
    converts list of dicts to dict of lists
    """
    reshape_LOD = list(zip(*(d.values() for d in LOD)))
    return {k:np.array(v,dtype=np.float64) for k,v in zip(LOD[0].keys(),reshape_LOD)}

def _return_np_type(x):
    return np.min_scalar_type((1 if np.min(x) >= 0 else -1) * np.max(np.abs(x)))

def get_config_pops(config):
    # make populations_by_letter
    if config['gimble']['task'] != 'simulate':
        config['populations']['population_by_letter'] = {
            'A': config['populations']['A'],
            'B': config['populations']['B'],
            }
    return config

def config_to_meta(config, task):
    '''
    - Function extracts those fields from 'config' that are meant to be saved as ZARR meta
    - ZARR does not like np.int64 but np.float64 is ok 
    '''
    meta = {}
    if task == 'tally':
        meta['data_key'] = config['data_key']
        meta['data_source'] = config['data_source']
        meta['data_type'] = config['data_type']
        meta['tally_key'] = config['tally_key']
        meta['max_k'] = None if config['max_k'] is None else tuple([int(v) for v in config['max_k']])
        meta['sample_sets'] = config['sample_sets']
        meta['sequences'] = config['sequences']
        meta['genome_file'] = config['genome_file']
        meta['blocks'] = config['blocks']
        meta['windows'] = config['windows']
        meta['marginality'] = config['marginality']
        meta['block_length'] = config['block_length']
    if task == 'simulate':
        meta['max_idx'] = config['idx']
        meta['max_k'] = tuple([int(v) for v in config['max_k']])
        meta['replicates'] = config['replicates']
        meta['parameters_LOD'] = config['parameters_LOD']
        meta['parameters'] = {k:(tuple(v.tolist()) if isinstance(v, np.ndarray) else v) for k,v in config['parameters'].items()}
        meta['discrete_genome'] = config['simulate']['discrete_genome']
        if config['gridbased']['grid_label']!="":
            meta['grid_label'] = config['gridbased']['grid_label']
            meta['window_param_idx'] = tuple(int(i) for i in config['window_param_idx'])
    if task == 'simulate_instance':
        meta['max_k'] = tuple([int(v) for v in config['max_k']])
        meta['idx'] = config['idx']
        meta['ancestry_seeds'] = tuple([int(s) for s in config['seeds'][config['idx']][:,0]]) 
        meta['mutation_seeds'] = tuple([int(s) for s in config['seeds'][config['idx']][:,1]])
        meta['replicates'] = config['replicates']
        meta['parameters_LOD'] = config['parameters_LOD'][config['idx']]
        meta['parameters'] = {k:(tuple(v.tolist()) if isinstance(v, np.ndarray) else v) for k,v in config['parameters'].items()}
        meta['discrete_genome'] = config['simulate']['discrete_genome']
        if not(isinstance(config['parameters']['recombination_rate'], float) or isinstance(config['parameters']['recombination_rate'], int)):
            meta['parameters']['recombination_rate'] = config['parameters']['recombination_rate'][config['idx']]
    if task == 'makegrid':
        meta['makegrid_key'] = config['key']
        meta['grid_dict'] = {k:list(v) for k,v in config['parameters_expanded'].items()}
        meta['block_length'] = config['mu']['block_length']
        meta['mu'] = config['mu']['mu']
        meta['label'] = config['gimble']['label']
        meta['model'] = config['gimble']['model']
        meta['reference_pop_id'] = config['populations']['reference_pop_id']      
        meta['sync_pop_ids'] = config['populations']['sync_pop_ids']
        meta['population_by_letter'] = config['populations']['population_by_letter']
        meta['max_k'] = list([float(k) for k in config['max_k']])
    if task == 'gridsearch':
        meta['makegrid_key'] = config['makegrid_key']
        meta['batch_sites'] = config['batch_sites']
        meta['data_key'] = config['data_key']
        meta['gridsearch_key'] = config['gridsearch_key']
        meta['grid_dict'] = config['grid_dict']
        meta['block_length'] = config['block_length_data']
        meta['data_label'] = config['data_label']
        meta['data_source'] = config['data_source']
        meta['gridsearch_keys'] = config['gridsearch_keys']
    if task == 'optimize':
        meta['optimize_key'] = config['optimize_key']
        meta['optimize_keys'] = config['optimize_keys']
        meta['random_seed'] = config['gimble']['random_seed']
        meta['data_source'] = config['data_source']
        meta['data_key'] = config['data_key']
        meta['start_point_method'] = config['start_point_method']
        meta['start_point'] = list([float(k) for k in config['start_point']])
        meta['population_by_letter'] = config['populations']['population_by_letter']
        meta['reference_pop_id'] = config['populations']['reference_pop_id']
        meta['optimize_time'] = config['optimize_time']
        meta['max_iterations'] = config['max_iterations']
        meta['xtol_rel'] = config['xtol_rel']
        meta['ftol_rel'] = config['ftol_rel']
        meta['parameters'] = config['parameters']
        meta['parameters']['fixed'] = config['parameters_fixed']
        meta['parameters']['bounded'] = config['parameters_bounded']
        meta['block_length'] = config['block_length']
        meta['mu'] = config['mu']['mu']
        meta['label'] = config['gimble']['label']
        meta['model'] = config['gimble']['model']
        meta['reference_pop_id'] = config['populations']['reference_pop_id'] 
        meta['sync_pop_ids'] = config['populations']['sync_pop_ids']
    return meta

def get_config_model_parameters(config, module):
    # Convert parameters to numpy arrays
    config['parameters_np'] = {}
    config['parameters_fixed'] = []
    config['parameters_bounded'] = [] # only first sync'ed pop is added to bounded
    config['parameters_gridded'] = []
    sync_flag = False # only allows one sync'ed pop 
    for parameter, values in config['parameters'].items():
        if len(values) == 1:
            config['parameters_fixed'].append(parameter)
            config['parameters_np'][parameter] = np.array(values)
        elif len(values) == 2:
            if module == 'makegrid':
                sys.exit("[X] Module %r only supports FLOAT, or (min, max, steps, lin|log) for parameters. Not %r" % (module, values))
            if parameter.startswith('Ne'):
                pop_id = parameter.replace("Ne_", "")
                if pop_id in config['populations']['sync_pop_ids']:
                    if not sync_flag:
                        parameter_sync = 'Ne_s'
                        config['parameters_bounded'].append(parameter_sync)    
                        config['parameters_np'][parameter_sync] = np.array(values)
                        sync_flag = True
                else:
                    config['parameters_bounded'].append(parameter)
                    config['parameters_np'][parameter] = np.array(values)
            else:
                config['parameters_bounded'].append(parameter)
                config['parameters_np'][parameter] = np.array(values)
        elif len(values) == 4:
            if module == 'optimize':
                sys.exit("[X] Module %r only supports FLOAT, or (MIN, MAX) for parameters. Not %r" % (module, values))
            value_min, value_max, value_num, value_scale = values
            config['parameters_gridded'].append(parameter)
            if value_scale.startswith('lin'):
                config['parameters_np'][parameter] = np.linspace(
                    value_min, value_max, num=value_num, endpoint=True, dtype=np.float64)
                if not value_min == 0 and parameter == 'me':
                    np.insert(config['parameters_np'][parameter], 0, 0)
            elif value_scale.startswith('log'):
                if value_min == 0:
                    error_msg = ["[X] Min value for log-ranged parameter %r in config file can't be 0." % parameter]
                    if parameter == 'me':
                        error_msg.append("[X] Pick the smallest Non-zero number to include (Zero will be added automatically)")
                    sys.exit("\n".join(error_msg))
                config['parameters_np'][parameter] = np.geomspace(
                    value_min, value_max, num=value_num, endpoint=True, dtype=np.float64)
                if parameter == 'me':
                    config['parameters_np'][parameter] = np.insert(config['parameters_np'][parameter], 0, 0)
            else:
                sys.exit("[X] Config: Scale should either be lin or log. Not %r." % value_scale)
        else:
            sys.exit("[X] Config: Parameters must be FLOAT, or (MIN, MAX), or (MIN, MAX, STEPS, LIN|LOG).") 
    #print('Config',config)
    return config

def get_config_model_events(config):
    pop_path = 'simulate' if config['gimble']['task'] == 'simulate' else 'populations'
    pop_ids = config[pop_path]['pop_ids']
    model = config['gimble']['model']
    sync_pops = config[pop_path].get('sync_pop_ids', [])
    events = get_ini_model_events(model, pop_ids)
    config['events'] = {}
    config['events']['coalescence'] = []
    #config['events']['p_to_c'] = {}
    for event in events:
        if event.startswith('M'):
            config['events']['migration'] = [tuple(pop_ids.index(pop_id) for pop_id in event.replace('M_', '').split('_'))]
        if event.startswith('J'):
            config['events']['exodus'] = [(0,1,2)] # populations are moving to the last element
        if event.startswith('C'):
            #pop_size = event.replace('C_', 'Ne_')
            if event.replace('C_', '') in sync_pops:
                event = 'C_s'
            config['events']['coalescence'].append(event)
        #    config['events']['p_to_c'][pop_size] = event
    return config

def get_config_kmax(config):
    try:
        #GB: sorting seems to be the way to go here
        #mutype_labels, max_k = zip(*sorted(config[k_max].items()))
        #max_k can set to np.array afterwards
        config['max_k'] = np.array([config['k_max'][mutype] for mutype in MUTYPES])
    except KeyError as e:
        sys.exit('[X] Config: No k-max value found for mutype %r ' % e.args[0])
    return config

def expand_parameters(config):
    print(config)
    if 'populations' in config and len(config['populations']['sync_pop_ids'])>0:
        to_sync = ['Ne_%s' % pop for pop in config['populations']['sync_pop_ids'][1:]]
        sync_to = "Ne_%s" % config['populations']['sync_pop_ids'][0]
        parameters_np = {k:v for k,v in config['parameters_np'].items() if k not in to_sync}
    else:
        parameters_np = config['parameters_np']
        to_sync, sync_to = None, None
    cartesian_product = itertools.product(*parameters_np.values())
    rearranged_product = list(zip(*cartesian_product))
    config['parameters_grid_points'] = len(rearranged_product[0])
    config['parameters_expanded'] = {k: np.array(v, dtype=np.float64) for k, v in zip(parameters_np.keys(), rearranged_product)}
    if to_sync:
        for pop in to_sync:
            config['parameters_expanded'][pop] = config['parameters_expanded'][sync_to]
    return config

def get_config_optimize(config):
    # bounds 
    config['parameter_combinations_lowest'] = np.array([config['parameters_np'][k][0] for k in config['parameters_bounded']])
    config['parameter_combinations_highest'] = np.array([config['parameters_np'][k][1] for k in config['parameters_bounded']])
    return config
    
def get_config_simulate(config):
    # interpopulation sample pairs
    config['simulate']['comparisons'] = list(
        itertools.product(
            range(config['simulate']['sample_size_A']), 
            range(config['simulate']['sample_size_A'], 
                config['simulate']['sample_size_A'] + 
                config['simulate']['sample_size_B'])))

    if isinstance(config['simulate']['recombination_rate'], str):
                config['simulate']['recombination_rate'] = 0.0

    fixed_parameter = config['gridbased']['fixed_parameter'].strip()
    if fixed_parameter!='':
        if fixed_parameter not in config['parameters']:
            sys.exit('[X] Gridbased fixed_parameter should be one of the model parameters.')
    config['demographies'] = lib.simulate.make_demographies(config)
    # if config['simulate']['num_linked_blocks'] is not specified, default to config['simulate']['blocks']
    if not config['simulate']['num_linked_blocks']:
        config['simulate']['num_linked_blocks'] = config['simulate']['blocks']
    blocks_per_replicate = get_blocks_per_replicate(config['simulate']['blocks'], config['simulate']['num_linked_blocks']) 
    sequence_length_per_replicate = blocks_per_replicate * config['simulate']['block_length']
    config['simulate']['blocks_per_replicate'] = blocks_per_replicate
    config['simulate']['sequence_length'] = sequence_length_per_replicate
    config['replicates'] = blocks_per_replicate.size * config['simulate']['replicates']
    config['seeds'] = np.random.randint(1, 2 ** 32, (config.get('parameters_grid_points', 1), config['replicates'], 2))       
    config['parameters_LOD'] = DOL_to_LOD(config['parameters_expanded'])
    config['parameters'] = {**config['simulate'],**config['mu']} # is this necessary?
    return config

def get_blocks_per_replicate(num, part):
    num_entries = math.ceil(num/part)
    result = np.full(num_entries, fill_value=part, dtype=np.uint64)
    floor = num//part
    if num_entries != floor:
        result[-1] = num - floor * part 
    return result

def load_config(config_file, MODULE=None, CWD=None, VERSION=None):
    parser = configparser.ConfigParser(inline_comment_prefixes='#', allow_no_value=True)
    parser.optionxform = str # otherwise keys are lowercase
    parser.read(config_file)
    parsee = {s: dict(parser.items(s)) for s in parser.sections()}
    # BEGIN fallback if executed from tests
    MODULE = parsee['gimble']['task'] if MODULE is None else MODULE
    VERSION = parsee['gimble']['version'] if VERSION is None else VERSION
    # END fallback if executed from tests
    schema = get_config_schema(MODULE)
    validator = ConfigCustomNormalizer(schema, module=MODULE, purge_unknown=True)
    validator.validate(parsee)
    if not validator.validate(parsee):
        #print(validator.errors)
        validator_error_string = get_validator_error_string(validator.errors)
        sys.exit("[X] Problems were encountered when parsing INI config file:\n%s" % validator_error_string)
    config = validator.normalized(parsee)
    np.random.seed(config['gimble']['random_seed'])
    config = get_config_pops(config)
    config = get_config_kmax(config)
    config = get_config_model_events(config)
    config = get_config_model_parameters(config, MODULE)
    if MODULE == 'makegrid' or MODULE == 'simulate':
        config = expand_parameters(config)
    #config = expand_parameters(config) # needed for testing
    if MODULE == 'simulate':
        config = get_config_simulate(config)
    if MODULE == 'optimize':
        config = get_config_optimize(config)
    config['CWD'] = CWD
    #for k, v in config.items():
    #    print(k, '\t', v)
    if not VERSION == config['gimble']['version']:
        print("[-] Version conflict:\n\tgimble %s\n\t config INI %s" % (VERSION, config['gimble']['version']))
    return config

def config_to_demes_graph(config, idxs=None):
    #using demes: all events are defined forwards in time!
    pop_path = 'simulate' if config['gimble']['task'] == 'simulate' else 'populations'
    if idxs==None:
        idxs = range(config['parameters_grid_points'])
    for idx in idxs:
        b = demes.Builder(time_units='generations')
        if 'exodus' in config['events']:
            b.add_deme('A_B', epochs=[
                dict(
                    end_time=config['parameters_expanded']['T'][idx], 
                    start_size=config['parameters_expanded']['Ne_A_B'][idx], 
                    end_size=config['parameters_expanded']['Ne_A_B'][idx])
                ]
            )
            b.add_deme('A', ancestors=['A_B'],defaults=
                dict(epoch=dict(
                    start_size=config['parameters_expanded']['Ne_A'][idx], 
                    end_size=config['parameters_expanded']['Ne_A'][idx]))
                )
            b.add_deme('B', ancestors=['A_B'],defaults=
                dict(epoch=dict(
                    start_size=config['parameters_expanded']['Ne_B'][idx], 
                    end_size=config['parameters_expanded']['Ne_B'][idx]))
                )
        else:
            b.add_deme('A',defaults=
                dict(epoch=dict(
                    start_size=config['parameters_expanded']['Ne_A'][idx], 
                    end_size=config['parameters_expanded']['Ne_A'][idx]))
                )
            b.add_deme('B',defaults=
                dict(epoch=dict(
                    start_size=config['parameters_expanded']['Ne_B'][idx], 
                    end_size=config['parameters_expanded']['Ne_B'][idx]))
                )
        if 'migration' in config['events']:
            for mig_event in config['events']['migration']:
                destination, source = mig_event
                b.add_migration(
                    source=config[pop_path]['pop_ids'][source], 
                    dest=config[pop_path]['pop_ids'][destination], 
                    rate=config['parameters_expanded']['me'][idx])            
        yield b.resolve()

def config_to_demes_yaml(config, outputfile, idxs=None):
    CWD = config['CWD']
    outputstring = outputfile + '_{}.yaml'
    graphs = config_to_demes_graph(config, idxs)
    for idx, graph in enumerate(graphs):
        with open(os.path.join(CWD, outputstring.format(str(idx))), 'w') as file:
            demes.dump(graph, file, format='yaml')

class ConfigCustomNormalizer(cerberus.Validator):
    def __init__(self, *args, **kwargs):
        super(ConfigCustomNormalizer, self).__init__(*args, **kwargs)
        self.module = kwargs['module']

    def _normalize_coerce_pop_ids(self, value):
        pop_ids = [v for v in value.replace(" ", "").replace("|", ",").split(",") if v]
        return pop_ids if pop_ids else None

    def _normalize_coerce_reference_pop_id(self, value):
        if value in set(self.document['pop_ids']):
            return value
        self._error('reference_pop_id', 
            "invalid reference_pop_id: %r. Must be one of the following: %s" % (value, ", ".join(self.document['pop_ids'])))
        return None

    def _normalize_coerce_int(self, value):
        try:
            return int(float(value))
        except ValueError:
            return None

    def _normalize_coerce_sync_pop_ids(self, value):
        if value: 
            sync_pop_ids = frozenset(self._normalize_coerce_pop_ids(value))
            if sync_pop_ids:
                sync_pop_ids_valid_sets = set([frozenset(self._normalize_coerce_pop_ids(",".join(pops))) for pops in list(itertools.chain.from_iterable(
                    itertools.combinations(self.document['pop_ids'], r) for r in range(len(self.document['pop_ids'])+1)))[4:]])
                if sync_pop_ids in sync_pop_ids_valid_sets:
                    return sorted(sync_pop_ids)
                else:
                    sync_pop_ids_valid_strings = " or ".join([",".join(sorted(x)) for x in sync_pop_ids_valid_sets])
                    self._error('sync_pop_ids', 'invalid sync_pop_ids: %r. Must be one of the following: %s' % (str(value), sync_pop_ids_valid_strings))
        return []

    def _normalize_coerce_float_or_list(self, value):
        values = value.strip('()[]').replace(' ', '').split(",")
        try:
            if len(values) == 4 and values[-1] in set(['lin', 'log', 'LIN', 'LOG']):
                return [float(values[0]), float(values[1]), int(values[2]), values[3]]
            elif len(values) == 2:
                return [float(values[0]), float(values[1])]
            elif len(values) == 1:
                return [float(values[0])]
            else:
                return None
        except ValueError:
            return None

    def _normalize_coerce_float_or_empty(self, value):
        value = value.replace(" ", "")
        if value:
            try:
                return float(value)
            except ValueError:
                return None
        return value

    def _normalize_coerce_int_or_empty(self, value):
        value = value.replace(" ", "")
        if value:
            try:
                return int(value)
            except ValueError:
                return None
        return value

    def _normalize_coerce_path(self, value):
        if not value:
            return ''
        _path = pathlib.Path(value).resolve()
        return str(_path)
        if not _path.exists():
            if self.module == 'model':
                self._error(self.module, 'Must be a valid path to a model file. Not %r' % value)
            elif self.module == 'simulate':
                self._error(self.module, 'Must be a valid path to the recombination map. Not %r' % value)
            else:
                pass

    def _validate_notNoneInt(self, notNoneNumeric, field, value):
        """{'type':'boolean'}"""
        if value == None and notNoneNumeric:
            self._error(field, "Must be an int value or empty")
    
    def _validate_notNoneFloat(self, notNoneNumeric, field, value):
        """{'type':'boolean'}"""
        if value == None and notNoneNumeric:
            self._error(field, "Must be a float or empty")

    def _validate_notNone(self, notNone, field, value):
        """{'type':'boolean'}"""
        if not value and notNone:
            self._error(field, "Must be FLOAT, or (MIN, MAX), or (MIN, MAX, STEPS, LIN|LOG).")

    def _validate_isPath(self, isPath, field, value):
        """{'type':'boolean'}"""
        if value is None:
            return None
        _path = pathlib.Path(value).resolve()
        return str(_path)
        if not _path.exists():
            if field == 'model':
                self._error(field, 'Must be a valid path to a model file. Not %r' % value)
            else:
                self._error(field, 'Must be a valid path to the recombination map. Not %r' % value)


def get_config_schema(module):
    # [To Do] 
    # move to parameterObj file
    schema = {}
    schema['gimble'] = {
            'type': 'dict',
            'schema': {
                'version': {'required': True, 'empty':False, 'type': 'string'},
                'model': {'required': True, 'empty':False, 'type': 'string'},
                'task': {'required': True, 'empty':False, 'type': 'string'},
                'label': {'required': True, 'empty':False, 'type': 'string'},
                'precision': {'required': True, 'empty':False, 'type': 'integer', 'coerce': int},
                'random_seed': {'required': True, 'empty':False, 'type': 'integer', 'coerce': int}}}
    if not module == 'simulate':
        schema['populations'] = {
            'type': 'dict', 
            'schema': {
                'pop_ids': {'required': True, 'empty':False, 'type': 'list', 'coerce': 'pop_ids'},
                'A': {'required':True, 'empty':True, 'type':'string'},
                'B': {'required':True, 'empty':True, 'type':'string'},
                'reference_pop_id': {'required':True, 'empty':False, 'type': 'string', 'coerce': 'reference_pop_id'},
                'sync_pop_ids': {'required':False, 'empty': True, 'type': 'list', 'coerce': 'sync_pop_ids'},
                }}
    else:
        schema['populations'] = {
            'type': 'dict', 
            'schema': {
                'pop_ids': {'required': True, 'empty':False, 'type': 'list', 'coerce': 'pop_ids'},
                'sync_pop_ids': {'required':False, 'empty': True, 'type': 'list', 'coerce': 'sync_pop_ids'},
                }}
    schema['k_max'] = {
            'type': 'dict', 
            'valuesrules': {'required': True, 'empty':False, 'type': 'integer', 'min': 1, 'coerce':int}}
    if module == 'makegrid':
        schema['mu'] = {
            'type':'dict', 
            'schema':{
                'mu': {'required': True, 'empty':False, 'type': 'float', 'coerce':float},
                'block_length': {'required': True, 'empty': False, 'min': 1, 'coerce':int},
                }}
    else:
        schema['mu'] = {
            'type':'dict', 
            'schema':{
                'mu': {'required': True, 'empty':False, 'type': 'float', 'coerce':float},
                }}
    schema['parameters'] = {
            'type': 'dict', 'required':True, 'empty':False, 
            'valuesrules': {'coerce':'float_or_list', 'notNone':True}}
    if module == 'simulate':
        schema['simulate'] = {
            'type': 'dict',
            'schema': {
                'pop_ids': {'required': True, 'empty':False, 'type': 'list', 'coerce': 'pop_ids'},
                'sync_pop_ids': {'required':False, 'empty': True, 'type': 'list', 'coerce': 'sync_pop_ids'},
                'reference_pop_id': {'required':True, 'empty':False, 'type': 'string', 'coerce': 'reference_pop_id'},
                'ploidy': {'required':True,'empty': False, 'min': 1, 'coerce': int},
                'blocks': {'required':True, 'empty': False, 'type': 'integer', 'min': 1, 'coerce': 'int'},
                'block_length': {'required': True, 'empty':False, 'min': 1, 'type': 'integer', 'coerce': int},
                'num_linked_blocks': {'required':True, 'empty':True, 'min': 1, 'coerce': 'int_or_empty'},
                'replicates': {'required':True,'empty': False, 'type': 'integer', 'min': 1, 'coerce':int},
                'sample_size_A': {'required':True,'empty':False, 'type': 'integer', 'min': 1, 'coerce':int},
                'sample_size_B': {'required':True,'empty':False, 'type': 'integer', 'min': 1, 'coerce':int},
                'discrete_genome': {'empty':True, 'type': 'boolean', 'coerce':bool},
                'recombination_rate': {'empty': True, 'notNoneFloat':True, 'coerce':'float_or_empty', 'min': 0.0},
                'recombination_map': {'empty': True, 'type': 'string', 'coerce': 'path', 'dependencies':['number_bins', 'cutoff', 'scale']},
                'number_bins': {'empty': True, 'notNoneInt': True, 'coerce':'int_or_empty', 'min': 1},
                'cutoff': {'empty': True, 'notNoneFloat': True, 'coerce':'float_or_empty', 'min': 0.0},
                'scale': {'empty':True, 'type':'string', 'allowed':['lin', 'log']}
        }}
        schema['gridbased'] = {
            'type': 'dict',
            'schema': {
                'grid_label': {'required': False, 'empty':True,'type': 'string'},
                'fixed_parameter': {'required': False, 'empty':True,'type': 'string'}
            }
        }
    return schema

def _dict_product(parameter_dict):
    cartesian_product = itertools.product(*parameter_dict.values())
    rearranged_product = list(zip(*cartesian_product))
    return {k: np.array(v, dtype=np.float64) for k, v in zip(parameter_dict.keys(), rearranged_product)}

def recursive_get_size(path):
    """Gets size in bytes of the given path, recursing into directories."""
    if os.path.isfile(path):
        return os.path.getsize(path)
    if not os.path.isdir(path):
        return 0
    return sum(
        recursive_get_size(os.path.join(path, name))
        for name in os.listdir(path))

def parse_csv(csv_f='', dtype=[], usecols=[], sep=',', header=None):
    '''dtypes := "object", "int64", "float64", "bool", "datetime64", "timedelta", "category"'''
    df = pd.read_csv(csv_f, sep=sep, usecols=usecols, names=list(dtype.keys()), header=header, dtype=dtype)
    if df.isnull().values.any():
        sys.exit("[X] Bad file format %r." % csv_f)
    return df

# all formats should be in another file 

def format_query_meta(meta, ignore_long=False):
    LONG_THRESHOLD = 80
    lines = []
    for key, value in meta.items():
        if isinstance(value, int):
            formatted_value = format_count(value)
        elif isinstance(value, list):
            formatted_value = str(value)
        else:
            formatted_value = value
        if ignore_long:
            if len(formatted_value) < LONG_THRESHOLD:
                lines.append("[+]\t%s: %s" % (key, formatted_value))
        else:
            lines.append("[+]\t%s: %s" % (key, formatted_value))
    return "\n".join(lines)

def format_bases(bases):
    if bases in set(['-', 'N/A']):
        return bases
    return "%s b" % format(bases, ',d')

def format_percentage(fraction, precision=2):
    if fraction in set(['-', 'N/A']):
        return fraction
    return "{:.{}%}".format(fraction, precision)

def format_proportion(fraction, precision=2):
    if fraction in set(['-', 'N/A']):
        return fraction
    return "{:.{}f}".format(fraction, precision)

def format_count(count):
    if count in set(['-', 'N/A']):
        return count
    return "%s" % str(format(count, ',d'))

def format_time(seconds):
    hours, remainder = divmod(seconds, 3600)
    minutes, seconds = divmod(remainder, 60)
    return "{:0>2}h:{:0>2}m:{:06.3f}s".format(int(hours), int(minutes), seconds)

def format_bytes(size, precision=1):
    power = 2**10
    n = 0
    power_labels = {0 : ' B', 1: 'KB', 2: 'MB', 3: 'GB', 4: 'TB'}
    while size > power:
        size /= power
        n += 1
    return "{:.{}f} {}".format(size, precision, power_labels[n])

def get_n50_from_lengths(lengths):
    length_sorted = sorted(lengths, reverse=True)
    cum_sum = np.cumsum(length_sorted)
    half = int(sum(lengths) / 2)
    cum_sum_2 = min(cum_sum[cum_sum >= half])
    n50_idx = np.where(cum_sum == cum_sum_2)
    return length_sorted[int(n50_idx[0][0])]

def parse_intervals(bed_f, target_sequences, target_samples):
    intervals_df = parse_csv(
        csv_f=bed_f, 
        sep="\t", 
        usecols=[0, 1, 2, 4], 
        dtype={'sequence': 'category', 'start': 'int64', 'end': 'int64', 'samples': 'category'},
        header=None)
    # Subset target_sequences, sort, reindex
    intervals_df = intervals_df[intervals_df['sequence'].isin(target_sequences)].sort_values(
            ['sequence', 'start'], ascending=[True, True]).reset_index(drop=True)
    # Convert samples column to sample-matrix
    intervals_df = pd.concat([intervals_df, intervals_df.samples.str.get_dummies(sep=',').filter(
        target_samples)], axis=1).drop(columns=['samples'])
    # Subset target samples
    samples_in_df = [sample for sample in intervals_df.columns[3:]]
    target_samples_in_df = ordered_intersect(a=samples_in_df, b=target_samples, order='a')
    target_samples_not_in_df = set(target_samples).difference(set(target_samples_in_df))
    if target_samples_not_in_df:
         sys.exit("[X] Samples in SAMPLE_FILE not found in BED_FILE: %s" % ", ".join(list(target_samples_not_in_df)))
    non_target_samples_in_df = set(target_samples).difference(set(samples_in_df))
    intervals_df = intervals_df.drop(non_target_samples_in_df, axis=1)
    # Add length column
    intervals_df['length'] = (intervals_df['end'] - intervals_df['start'])
    return intervals_df

def bsfs_to_2d(bsfs):
    """Converts 4D bsfs to 2D array with counts, mutuples.
       Converts 5D bsfs to 2D array with window_idx, counts, mutuples.

    Parameters 
    ----------
    bsfs : ndarray, int, ndim (4 or 5)  
    
    Returns
    -------
    out : ndarray, int, ndim (2)
    """
    if not np.any(bsfs):
        return None
    non_zero_idxs = np.nonzero(bsfs)
    if bsfs.ndim == 4: # blocks
        return np.concatenate([bsfs[non_zero_idxs].reshape(non_zero_idxs[0].shape[0], 1), np.array(non_zero_idxs).T], axis=1)
    elif bsfs.ndim == 5: # windows
        non_zero_idxs_array = np.array(non_zero_idxs).T
        first = non_zero_idxs_array[:,0].reshape(non_zero_idxs[0].shape[0], 1)
        second = bsfs[non_zero_idxs].reshape(non_zero_idxs[0].shape[0], 1)
        third = non_zero_idxs_array[:,1:]
        return np.concatenate([first, second, third], axis=1)
    else:
        raise ValueError('bsfs_to_2d: bsfs.ndim must be 4 (blocks) or 5 (windows)')

def calculate_blocks_report_metrics(tally, sample_sets_count, block_length, intervals_span):
    BRM = collections.defaultdict(lambda: "-")
    BRM['blocks_total'] = np.sum(tally[:,0]) if np.any(tally) else 0
    effective_length = block_length * BRM['blocks_total']
    BRM['interval_coverage'] = (effective_length / sample_sets_count) / intervals_span #/ sample_sets_count
    if BRM['blocks_total']:
        BRM['blocks_invariant'] = tally[0,0] / BRM['blocks_total']
        BRM['blocks_fgv'] = np.sum(tally[(tally[:,3]>0) & (tally[:,4]>0)][:,0]) / BRM['blocks_total']
        hetB = np.sum((tally[:,0] * tally[:,1])) 
        hetA = np.sum((tally[:,0] * tally[:,2]))
        hetAB = np.sum((tally[:,0] * tally[:,3])) 
        fixed = np.sum((tally[:,0] * tally[:,4]))
        # not sure about heterozygosity ...
        BRM['heterozygosity_A'] = (hetA + hetAB) / effective_length 
        BRM['heterozygosity_B'] = (hetB + hetAB) / effective_length # hetB or hetA ?
        BRM['heterozygosity_intra'] = (float(fractions.Fraction(1, 2) * (hetA + hetB) + hetAB)) / effective_length 
        # BRM['heterozygosity_B'] = BRM['heterozygosity_A'] # this is how it was before...
        BRM['dxy'] = ((hetA + hetB + hetAB) / 2.0 + fixed) / effective_length
        mean_pi = (BRM['heterozygosity_A'] + BRM['heterozygosity_B']) / 2.0
        total_pi = (BRM['dxy'] + mean_pi) / 2.0 
        BRM['fst'] = (BRM['dxy'] - mean_pi) / (BRM['dxy'] + mean_pi) if (total_pi) else "N/A"
        total_segregating = np.sum(tally[:,0, None] * tally[:,1:])
        BRM['pi'] = float(fractions.Fraction(1, 2) * (hetA + hetB) 
            + fractions.Fraction(2, 3) * (hetAB + fixed)) / effective_length 
        BRM['watterson_theta'] = total_segregating / float(harmonic(3)) / effective_length
    return BRM

def write_info_report(version, report, prefix):
    # OUTPUTLIB
    txt = ["# %s" % version, str(report)]
    out_f = '%s.info.txt' % prefix
    print("[+] Writing info file %s ..." % out_f)
    with open(out_f, 'w') as out_fh:
        out_fh.write("\n".join(txt) + "\n")
    return out_f

def get_popgen_metrics(array, sites=0):
    '''only works for mutypes=4
    Possible inputs: 
    A) block-tally  := ndim (2); shape (n, 5)
    B) window-tally := ndim (2); shape (n, 6)
    C) window-bsfs  := ndim (4); shape (maxk_m1, maxk_m2, maxk_m3, maxk_m4)
    D) window-bsfs  := ndim (5); shape (w, maxk_m1, maxk_m2, maxk_m3, maxk_m4)
    sites = (block_length * block_count)
    '''
    assert sites > 0, "sites must be positive integer"
    array = array if array.ndim == 2 else bsfs_to_2d(array)
    if array.shape[1] == 5: # block tally
        mutype_array = np.array(np.sum(array[:,0, np.newaxis] * array[:, 1:], axis=0)).reshape(-1, 4) # must be reshaped to 2D, so that indexing works
    elif array.shape[1] == 6: # window tally
        mutype_array = np.vstack([np.bincount(array[:, 0], weights=(array[:, 1] * array[:, (2 + m_idx)])) for m_idx in range(4)]).T
    else:
        raise ValueError('get_popgen_metrics: if tally.ndim is 2, then shape must be (n, 5) (blocks) or (n, 6) (windows)')
    heterozygosity_A = (mutype_array[:,1] + mutype_array[:,2]) / sites
    heterozygosity_B = (mutype_array[:,0] + mutype_array[:,2]) / sites
    d_xy = ((mutype_array[:,1] + mutype_array[:,0] + mutype_array[:,2]) / 2.0 + mutype_array[:,3]) / sites
    mean_pi = (heterozygosity_A + heterozygosity_B) / 2.0
    f_st = np.full(mutype_array.shape[0], np.nan)
    np.true_divide((d_xy - mean_pi), (d_xy + mean_pi), out=f_st, where=((d_xy + mean_pi) > 0))
    return np.vstack([heterozygosity_A, heterozygosity_B, d_xy, f_st])

def pop_metrics_from_bsfs(bsfs, block_length=None, window_size=None):
    warnings.warn("lib.gimble.pop_metrics_from_bsfs() is deprecated. ...", DeprecationWarning)
    '''only works for mutypes=4'''
    if not bsfs.ndim == 2:
        bsfs = bsfs_to_2d(bsfs)
    mutype_array = np.vstack([np.bincount(bsfs[:, 0], weights=bsfs[:, 1] * bsfs[:, (2 + m_idx)]) for m_idx in range(4)]).T 
    heterozygosity_A = (mutype_array[:,1] + mutype_array[:,2]) / (block_length * window_size)
    heterozygosity_B = (mutype_array[:,0] + mutype_array[:,2]) / (block_length * window_size)
    d_xy = ((mutype_array[:,1] + mutype_array[:,0] + mutype_array[:,2]) / 2.0 + mutype_array[:,3]) / (block_length * window_size)
    mean_pi = (heterozygosity_A + heterozygosity_B) / 2.0
    f_st = np.full(mutype_array.shape[0], np.nan)
    np.true_divide((d_xy - mean_pi), (d_xy + mean_pi), out=f_st, where=(d_xy + mean_pi) > 0)
    pop_metrics = np.vstack([heterozygosity_A, heterozygosity_B, d_xy, f_st])
    return pop_metrics

def check_unique_pos(pos_array):
    unique_pos, counts_pos = np.unique(pos_array, return_counts=True)
    duplicates = unique_pos[counts_pos > 1]
    if duplicates.any():
        print("\n[-] %s VCF records with non-unique positions found. Rescuing records by shifting position... (abort if this is not desired)" % (len(duplicates)))
        pos_array = fix_pos_array(pos_array)
    return pos_array

def fix_pos_array(pos_array):
    '''
    De-duplicates array by shifting values forward until there aren't any collisions
    '''
    # get boolean array for first and subsequent duplicates (True) (required sorted) 
    idxs = np.insert((np.diff(pos_array)==0).astype(bool), 0, False)
    if np.any(idxs): 
        # if there are duplicates, get new values by incrementing by one
        new_values = pos_array[idxs]+1
        # get non-duplicate values
        uniq_values = pos_array[~idxs]
        # insert new_values in non-duplicated values (required sorted)
        new_idxs = np.searchsorted(uniq_values, new_values)
        # recursive call
        return fix_pos_array(np.sort(np.insert(uniq_values, new_idxs, new_values)))
    # if there are no duplicated values
    return pos_array

def szudzik_pairing(folded_minor_allele_counts):
    # adapted from: https://drhagen.com/blog/superior-pairing-function/
    return np.where(
            (folded_minor_allele_counts[:,0] >= folded_minor_allele_counts[:,1]),
            (np.square(folded_minor_allele_counts[:,0]) + 
                folded_minor_allele_counts[:,0] + folded_minor_allele_counts[:,1]),
            (np.square(folded_minor_allele_counts[:,1]) + 
                folded_minor_allele_counts[:,0])
            )

def _harmonic(a, b):
    if b-a == 1:
        return fractions.Fraction(1,a)
    m = (a+b)//2
    return _harmonic(a,m) + _harmonic(m,b)

def harmonic(n):
    '''https://fredrik-j.blogspot.com/2009/02/how-not-to-compute-harmonic-numbers.html'''
    return _harmonic(1,n+1)

def chisq(sample_set_idxs, window_samples_set_idxs):
    if window_samples_set_idxs.size == 0:
        return 0.0
    spacer = (np.max(window_samples_set_idxs)+1)
    window_count = window_samples_set_idxs.shape[0]
    window_size = window_samples_set_idxs.shape[1]
    temp_sites = window_samples_set_idxs + (spacer * np.arange(window_count, dtype=np.int64).reshape(window_count, 1))
    obs = np.bincount(temp_sites.ravel(), minlength=(window_count * spacer)).reshape(-1, spacer)[:,sample_set_idxs]
    #print(obs)
    exp = np.full(obs.shape, window_size/sample_set_idxs.shape[0])
    return np.sum((((obs-exp)**2)/exp), axis=1) # / sample_sets.shape[0]

def mse(sample_set_idxs, window_samples_set_idxs):
    '''measure of eveness'''
    if window_samples_set_idxs.size == 0:
        return 0.0
    spacer = (np.max(window_samples_set_idxs)+1)
    window_count = window_samples_set_idxs.shape[0]
    window_size = window_samples_set_idxs.shape[1]
    temp_sites = window_samples_set_idxs + (spacer * np.arange(window_count, dtype=np.int64).reshape(window_count, 1))
    obs = np.bincount(temp_sites.ravel(), minlength=(window_count * spacer)).reshape(-1, spacer)[:,sample_set_idxs]
    exp = np.full(obs.shape, window_size/sample_set_idxs.shape[0])
    # Gertjan: scale by max
    max_mse_obs = np.zeros(sample_set_idxs.shape[0])
    max_mse_obs[0] = window_size 
    max_mse = np.sum(((max_mse_obs-exp)**2), axis=1)
    return np.sum((((obs-exp)**2)), axis=1) / max_mse

def blocks_to_windows(sample_set_idxs, block_variation, start_array, end_array, block_sample_set_idxs, window_size, window_step):
    # order of blocks is defined by end_array
    # coordinate_sorted_idx is the order of blocks if one were to sort them by end_array 
    coordinate_sorted_idx = np.argsort(end_array) 
    # elements in windows are defined by window_idxs -> shape(n, window_size) 
    window_idxs = np.arange(coordinate_sorted_idx.shape[0] - window_size + 1)[::window_step, None] + np.arange(window_size)
    #window_idxs_alt = np.arange((block_variation.shape[0] - window_size) + 1)[::window_step, None] + np.arange(window_size)
    # all taking is done with coordinate_sorted_idx and window_idxs
    window_variation = block_variation.take(coordinate_sorted_idx, axis=0).take(window_idxs, axis=0)
    block_starts = start_array.take(coordinate_sorted_idx, axis=0).take(window_idxs, axis=0)
    window_starts = np.min(start_array.take(coordinate_sorted_idx, axis=0).take(window_idxs, axis=0), axis=1).T
    block_ends = end_array.take(coordinate_sorted_idx, axis=0).take(window_idxs, axis=0)
    window_ends = np.max(end_array.take(coordinate_sorted_idx, axis=0).take(window_idxs, axis=0), axis=1).T
    # needs some solution for chisq-calculation by window ...
    window_samples_set_idxs = block_sample_set_idxs.take(coordinate_sorted_idx, axis=0).take(window_idxs, axis=0)
    #np.set_printoptions(threshold=sys.maxsize)
    balance = chisq(sample_set_idxs, window_samples_set_idxs)
    mse_sample_set_cov = mse(sample_set_idxs, window_samples_set_idxs)
    block_midpoints = (block_starts / 2) + (block_ends / 2)
    window_pos_mean = np.rint(np.mean(block_midpoints, axis=1).T)
    window_pos_median = np.rint(np.median(block_midpoints, axis=1).T)
    return (window_variation, window_starts, window_ends, window_pos_mean, window_pos_median, balance, mse_sample_set_cov)

def block_sites_to_variation_arrays(block_sites, cols=np.array([1,2,3]), max_type_count=7):
    temp_sites = block_sites + (max_type_count * np.arange(block_sites.shape[0], dtype=np.int64).reshape(block_sites.shape[0], 1))
    multiallelic, missing, monomorphic, variation = np.hsplit(np.bincount(temp_sites.ravel(), minlength=(block_sites.shape[0] * max_type_count)).reshape(-1, max_type_count), cols)
    return (multiallelic, missing, monomorphic, variation)

def gt2fmac(genotype_array):
    '''turns scikit-allel/numpy-array (n,s,p) of genotypes into folded-minor-allel-count-np-array (n,s)
    [To Do] multiallelics have to be done differently if we use non-(n,2,2)-GTs ... '''
    np_genotype_array = np.array(genotype_array)
    sa_genotype_array = allel.GenotypeArray(genotype_array)
    n_samples = sa_genotype_array.n_samples
    n_variants = sa_genotype_array.n_variants
    np_allele_count_array = np.ma.masked_equal(sa_genotype_array.count_alleles(), 0, copy=True)
    missing_genotypes = sa_genotype_array.is_missing()
    if np.all(missing_genotypes):
        folded_minor_allele_counts = np.ones((n_variants, n_samples),dtype=np.int8) * -1 # -1, -1 for missing => -1
    else:
        idx_max_global_allele_count = np.nanargmax(np_allele_count_array, axis=1)
        idx_min_global_allele_count = np.nanargmin(np_allele_count_array, axis=1)
        has_major_allele = (idx_max_global_allele_count != idx_min_global_allele_count)
        idx_min_prime_allele = np.amin(np_genotype_array[:,0], axis=1)
        idx_min_global_allele = np.amin(np.amin(np_genotype_array, axis=1), axis=1)
        idx_max_global_allele = np.amax(np.amax(np_genotype_array, axis=1), axis=1)
        idx_major_allele = np.where(has_major_allele, idx_max_global_allele_count, idx_min_prime_allele)
        idx_minor_allele = np.where(has_major_allele, idx_min_global_allele_count, np.where((
            idx_min_global_allele == idx_min_prime_allele),
            np.max((idx_min_global_allele, idx_max_global_allele), axis=0), 
            np.min((idx_min_global_allele, idx_max_global_allele), axis=0)))
        allele_map = np.ones((np_allele_count_array.shape), dtype=np.int64) * np.arange(np_allele_count_array.shape[-1], dtype=np.int8) # changed to np.int8
        # for each genotype (np.arange(allele_map.shape[0])), set minor allele to 1 (1st do minor, so that overwritten if monomorphic)
        allele_map[np.arange(allele_map.shape[0]), idx_minor_allele] = 1 
        # for each genotype (np.arange(allele_map.shape[0])), set major allele to 0
        allele_map[np.arange(allele_map.shape[0]), idx_major_allele] = 0
        folded_minor_allele_counts = sa_genotype_array.map_alleles(allele_map).to_n_alt(fill=-1)
        folded_minor_allele_counts[np.any(missing_genotypes, axis=1)] = np.ones(n_samples) * -1 # -1, -1 for missing => -1
        # multiallelics take precedence over missing
        folded_minor_allele_counts[(np_allele_count_array.count(axis=1) > 2)] = np.ones(2) * (-1, -2) # -1, -2 for multiallelic => -2
    return folded_minor_allele_counts

def intervals_to_sites(intervals):
    '''starts = np.array([0, 5, 8, 11, 15])
       ends   = np.array([2, 8, 9, 13, 18])
    clens : array([2, 5, 6, 8, 11])
    _sites: array([1, 1, 1, 1, 1, 1,  1, 1, 1, 1, 1])
    _sites: array([0, 1, 1, 1, 1, 1,  1, 1, 1, 1, 1])
    _sites: array([1, 1, 4, 1, 1, 1,  3, 1, 3, 1, 1])
    sites : array([0, 1, 5, 6, 7, 8, 11,12,15,16,17])'''
    starts, ends = intervals
    if starts is None or ends is None:
        return None 
    clens = np.cumsum(ends - starts) 
    if np.any(clens):
        sites = np.ones(clens[-1], dtype=np.int64)
        sites[0] = starts[0]
        sites[clens[:-1]] = starts[1:] - ends[:-1] + 1
        sites = sites.cumsum()
        return sites
    return None

def sites_to_blocks(sites, block_length, block_span, debug=False):
    if sites is None or not block_length or not block_span:
        return None
    max_gap = (block_span - block_length - 1)
    block_sites = np.concatenate([
        x[:block_length * (x.shape[0] // block_length)].reshape(-1, block_length) 
            for x in np.split(sites, np.where(np.diff(sites) > max_gap)[0] + 1)]) 
    #block_sites_valid_mask = (((block_sites[:, -1] - block_sites[:, 0]) + 1) <= block_span)
    block_sites_valid_mask = (((block_sites[:, -1] - block_sites[:, 0])) <= block_span)
    if debug:
        success = (block_sites[block_sites_valid_mask].shape[0] * block_length)/sites.shape[0] if np.count_nonzero(block_sites[block_sites_valid_mask]) else 0
        print(block_sites[block_sites_valid_mask])
        print('[+] sites=%s : blocks=%s success=%.2f%%' % (sites.shape[0], block_sites[block_sites_valid_mask].shape[0], success))
    if np.any(block_sites_valid_mask):
        return block_sites[block_sites_valid_mask]
    return None

def subset_gt_matrix(meta_seqs, sample_set, indices, gt_matrix):
    if gt_matrix is not None:
        sample_set_vcf_idxs = np.array([meta_seqs['variants_idx_by_sample'][sample] for sample in sample_set])
        return gt_matrix.subset(indices, sample_set_vcf_idxs)
    return None

def blocks_to_arrays(blocks, gts, pos):
    starts = np.array(blocks[:, 0], dtype=np.int64)
    ends = np.array(blocks[:, -1] + 1, dtype=np.int64) # BED: end+1
    pos_in_block_sites = np.isin(pos, blocks, assume_unique=True)
    if np.any(pos_in_block_sites): # if variants in blocks
        folded_minor_allele_counts = gt2fmac(gts)
        block_sites_in_pos = np.isin(blocks, pos, assume_unique=True)
        blocks[block_sites_in_pos] = szudzik_pairing(folded_minor_allele_counts) + 2 # add 2 so that not negative for bincount
        blocks[~block_sites_in_pos] = 2 # monomorphic = 2 (0 = multiallelic, 1 = missing)
        temp_sites = blocks + (7 * np.arange(blocks.shape[0], dtype=np.int64).reshape(blocks.shape[0], 1))
        multiallelic, missing, monomorphic, variation = np.hsplit(
            np.bincount(temp_sites.ravel(), minlength=(blocks.shape[0] * 7)).reshape(-1, 7), np.array([1,2,3]))
    else:
        multiallelic = np.zeros((blocks.shape[0], 1), dtype=np.int64)
        missing = np.zeros((blocks.shape[0], 1), dtype=np.int64)
        monomorphic = np.full((blocks.shape[0], 1), blocks.shape[1], dtype=np.int64)
        variation = np.zeros((blocks.shape[0], 4), dtype=np.int64)
    return (starts, ends, multiallelic, missing, monomorphic, variation)

def sum_wbsfs(bsfs_windows):
    assert bsfs_windows.ndim == 5, "only works for bsfs_windows.ndim = 5"
    return bsfs_windows.sum(axis=0)

def tally_to_bsfs(tally, max_k, data='blocks'):
    if data == 'blocks':
        counts = tally[:,0]
        dtype = _return_np_type(counts)
        out = np.zeros((np.array(max_k,dtype='int') + 2), dtype)
        out[tuple(tally[:,1:].T)] = counts
    elif data == 'windows':
        counts = tally[:,1]
        window_idxs = tally[:,0][np.newaxis]
        num_windows = np.max(window_idxs) + 1
        out_shape = np.insert((np.array(max_k,dtype='int') + 2), 0, num_windows)
        dtype = _return_np_type(counts)
        out = np.zeros(out_shape, dtype)
        idx_and_counts = np.hstack((window_idxs.T, tally[:,2:]))
        out[tuple(idx_and_counts.T)] = counts

    else:
        raise ValueError(f'2d_to_bsfs not implemtened for {data}.')
    return out

def tally_variation(variation, form='bsfs', max_k=None):
    '''
    Parameters 
    ----------
    variation : ndarray, np.uint64, ndim (2 or 3)  

    form : 'bsfs' for bsfs tensors
           'tally' for 2D tally array

    max_k : ndarray, int, ndim (1) for capping mutypes
            defaults to np.array([8,8,8,8]) if (max_k == None & form == 'bsfs') 
    Returns
    -------
    out : ndarray, int
          ndim (4 or 5) if form == 'bsfs'
          ndim (2 or 3) if form == 'tally'
    '''
    if max_k is None:
        max_k = np.array([8,8,8,8]) if form == 'bsfs' else None
    else:
        max_k = max_k + 1 # for clipping
    if variation.ndim == 2: 
        mutuples = np.clip(variation, 0, max_k)
    elif variation.ndim == 3:
        index = np.repeat(np.arange(variation.shape[0]), variation.shape[1]).reshape(-1, 1)
        mutuples = np.concatenate((index, np.clip(variation.reshape(-1, variation.shape[-1]), 0, max_k)), axis=-1)
    else:
        raise ValueError('variation.ndim is %r, should either be 2 (blocks) or 3 (windows)' % variation.ndim)
    try:
        mutuples_unique, counts = np.unique(mutuples, return_counts=True, axis=0)
        dtype = _return_np_type(counts)
        if form == 'bsfs':
            out = np.zeros(tuple(np.max(mutuples, axis=0) + 1), dtype) 
            out[tuple(mutuples_unique.T)] = counts
        elif form == 'tally':
            out = np.concatenate((counts.reshape(counts.shape[0], 1), mutuples_unique), axis=1)
            if variation.ndim == 3:
                out[:, [0, 1]] = out[:, [1, 0]] # for window variation, order must be [idx, count, mutuple]
            out = out.astype(dtype)
        else:
            raise ValueError('form must be %r or %r, was %r' % ('bsfs', 'tally', form))    
    except MemoryError as e:
        sys.exit('[X] tally_variation() ran out of memory. Try specifying lower k-max values. %s.' % str(e))
    return out

def calculate_marginality_of_variation(data, max_k=None):
    '''to be run on variation arrays of ndim = 2 or 3'''
    assert (
        data.shape[-1] == 4 and 
            (data.ndim == 3 or data.ndim == 2)), '[X] data.ndim must be 2 or 3, data.shape[-1] must be 4'
    if max_k is None:
        return 0.0
    is_above_max_k = np.any((data- max_k) > 0, axis=-1)
    return np.sum(is_above_max_k) / is_above_max_k.flatten().shape[0]

def ordered_intersect(a=[], b=[], order='a'):
    # [GIMBLE] unless used somewhere else
    A, B = a, b
    if order == 'b' :
        B, A = a, b
    return [_a for _a in A if _a in set(B)]

def get_hash_from_dict(d):
    # probably not needed
    '''returns md5sum hash of str(dict)'''
    if isinstance(d, dict):
        return hashlib.md5(str(d).encode()).hexdigest()
    raise ValueError('must be a dict')

def grid_meta_dict_to_value_arrays_by_parameter(grid_meta_dict):
    # [INPUTLIB]
    #see function lib.gimble.LOD_to_DOL() -> should be redundant
    _values_by_parameter = collections.defaultdict(list)
    for grid_idx, grid_dict in grid_meta_dict.items():
        for key, value in grid_dict.items():
            _values_by_parameter[key].append(value)
    values_by_parameter = {}
    for key, values in _values_by_parameter.items():
        values_by_parameter[key] = np.array(values)
    return values_by_parameter

def _optimize_to_csv(results, label, parameterObj, data_type='simulate'):
    if data_type=='simulate':
        df = pd.DataFrame(results[1:])
        df.columns=results[0]
        df = df.sort_values(by='iterLabel')
        df.set_index('iterLabel', inplace=True)
        _optimize_describe_df(df, label)
    else:
        headers = results.pop(0)
        df = pd.DataFrame(results, columns=headers)
        if parameterObj.trackHistory:
            df['step_id'] = df['step_id'].astype(int)
    df_with_scaling = _optimize_scale_output(parameterObj, df)
    df_with_scaling.to_csv(f'{label}_optimize_result.csv')

def _optimize_scale_output(parameterObj, df):
    parameter_dict = df.to_dict(orient='list')
    reference_pop_id = parameterObj.config['populations']['reference_pop_id']
    block_length = parameterObj.config['mu']['blocklength']
    mu = parameterObj.config['mu']['mu']
    scaled_result = lib.math.scale_parameters(
        parameter_dict, 
        reference_pop_id, 
        block_length, 
        mu
        )
    scaled_df = pd.DataFrame(scaled_result).astype(float)
    return pd.concat([df, scaled_df], axis=1)

def _optimize_describe_df(df, label):
    summary=df.drop(labels=['lnCL', 'exitcode'], axis=1).describe(percentiles=[0.025,0.975])
    summary.to_csv(f'{label}_summary.csv')

class ParameterObj(object):
    # [INPUTLIB]
    '''Superclass ParameterObj'''
    def __init__(self, params):
        self._PATH = params['path']
        self._VERSION = params['version']
        self._MODULE = params['module']
        self._CWD = params['cwd']
        self.config = None

    def __repr__(self):
        return("[+] VER := %s\n[+] CWD := %s\n[+] CMD := %s\n" % (
            self._VERSION, 
            self._CWD, 
            self._get_cmd()
            ))

    def _dict_product(self, parameter_dict):
        cartesian_product = itertools.product(*parameter_dict.values())
        rearranged_product = list(zip(*cartesian_product))
        return {k:np.array(v,dtype=np.float64) for k,v in zip(parameter_dict.keys(), rearranged_product)}

    def _dict_zip(self, pdict):
        '''DRL: if this is only used once, no need for being separate function'''
        return [dict(zip(pdict, x)) for x in zip(*pdict.values())]

    def _get_blocks_length(self):
        warnings.warn("lib.gimble._get_blocks_length() is deprecated.", DeprecationWarning)
        blocks_length_zarr = None
        blocks_length_ini = self._get_int(self.config['mu']['blocklength'], ret_none=True)
        if self.zstore:
            gimbleStore = lib.gimble.Store(path=self.zstore, create=False, overwrite=False)
            meta_blocks = gimbleStore._get_meta('blocks')
            if isinstance(meta_blocks, dict):
                blocks_length_zarr = meta_blocks.get('length', None)
        if blocks_length_zarr and blocks_length_ini:
            if blocks_length_zarr != blocks_length_ini:
                print("[-] Block length in INI and gimbleStore differ. Using block length from INI : %s b" % blocks_length_ini)
                return blocks_length_ini
            return blocks_length_zarr
        if blocks_length_zarr:
            return blocks_length_zarr
        if blocks_length_ini:
            return blocks_length_ini
        else:
            sys.exit("[X] Blocklength needs to be specified in ini file.")

    def _get_int(self, string, ret_none=False):
        try:
            return int(string)
        except TypeError:
            if not ret_none:
                sys.exit("[X] %r can't be converted to integer." % string)
            return None
        except ValueError:
            if not ret_none:
                sys.exit("[X] %r can't be converted to integer." % string)
            return None

    def _get_float(self, string, ret_none=False):
        try:
            return float(string)
        except TypeError:
            if not ret_none:
                sys.exit("[X] %r can't be converted to float." % string)
            return None

    def _get_cmd(self):
        return "%s/gimble %s %s" % (
            self._PATH, 
            self._MODULE, 
            "".join(["--%s " % " ".join((k, str(v))) for k,v in self.__dict__.items() if not k.startswith("_")]))

    def _get_fixed_params(self, subgroup=None, as_dict=False):
        if not as_dict:
            fixed_params =  [param for param,value in self.config['parameters'].items() if isinstance(value,float) or isinstance(value,int)]
            if subgroup:
                fixed_params = [param for param in fixed_params if param.startswith(subgroup)]
            return fixed_params
        else:
            fixedParams = {k:next(iter(v)) for k, v in self.parameter_combinations.items() if len(set(v))==1}
            fixedParams['mu'] = self.config['mu']['mu']
            if subgroup:
                fixedParams = {k:v for k,v in fixedParams.items() if k.startswith(subgroup)}
            return fixedParams

    def _get_max_k(self, kmax_string):
        #mutypes = ['m_1', 'm_2', 'm_3', 'm_4']
        error = "[X] Invalid k-max string (must be a list of 4 integers)"
        if kmax_string == None:
            return None
        try:
            kmax = np.array(ast.literal_eval(kmax_string))
            return kmax if kmax.shape == (4,) else sys.exit(error)
        except ValueError:
            sys.exit(error) 

    def _get_path(self, infile, path=False):
        if infile is None:
            return None
        _path = pathlib.Path(infile).resolve()
        if not _path.exists():
            sys.exit("[X] File not found: %r" % str(infile))
        if path:
            return _path
        return str(_path)

    def _get_prefix(self, prefix):
        if prefix is None:
            return None
        path = pathlib.Path(prefix).resolve()
        new_path = pathlib.Path(prefix+'.z').resolve()
        if new_path.exists():
            sys.exit("[X] zstore already exists. Specify using -z or provide new prefix.")
        parent_path = pathlib.Path(prefix).resolve().parent
        if not parent_path.exists():
            sys.exit("[X] Path does not exist: %r" % str(parent_path))
        return str(path)

    def _get_pops_to_sync(self, config=None, valid_sync_pops=None):
        print('_get_pops_to_sync_old')
        reference, to_be_synced = None, None
        print('valid_sync_pops', valid_sync_pops)
        if config:
            syncing = config['populations']['sync_pop_sizes']
            reference_pop = config['populations']['reference_pop']
        else:
            syncing =  self.config['populations']['sync_pop_sizes']
        print('reference_pop', reference_pop)
        print('syncing', syncing)
        if syncing:
            if len(syncing)>0:
                syncing = syncing.split(',')
                reference = syncing[0]
                to_be_synced = syncing[1:]
                if self.config['populations']['reference_pop'] in to_be_synced:
                    sys.exit(f"[X] Set reference pop to {reference}.")
                reference_size = config['parameters'][f'Ne_{reference}']
                tBS_sizes = [config['parameters'][f'Ne_{pop}'] for pop in to_be_synced]
                reference_size = [s for s in tBS_sizes if s!=None and s!=reference_size and s!='']
                if len(reference_size)>0:
                    sys.exit(f"[X] Syncing pop sizes: set no value or the same value for Ne_{', Ne_'.join(to_be_synced)} as for Ne_{reference}")         
                if self._MODULE == 'optimize':
                    fixed_Nes = self._get_fixed_params(subgroup='Ne')
                    if len(fixed_Nes)>0:
                        if not f"Ne_{reference_pop}" in fixed_Nes:
                            sys.exit("[X] No. No. No. It would make much more sense to set a population with a fixed size as reference.")
        print('reference', reference)
        print('to_be_synced', to_be_synced)
        return (reference, to_be_synced)

    def _get_pops_to_sync_short(self):
        syncing_to, to_be_synced = None, []
        syncing = self.config['populations']['sync_pop_sizes']
        if syncing and syncing.strip(' ') !='':        
            syncing_to, *to_be_synced = syncing.split(',')
        return (syncing_to, to_be_synced)

    def _get_unique_hash_from_dict(self, d):
        '''active'''
        return hashlib.md5(str(d).encode()).hexdigest()

    def _get_unique_hash(self, return_dict=False, module=None):
        '''passive'''
        module = module if module else self._MODULE
        to_hash = copy.deepcopy(self.config)
        #print('to_hash', to_hash)
        if module in set(['makegrid','gridsearch', 'query']):
            del to_hash['simulations']
            if 'kmax_by_mutype' in to_hash:
                del to_hash['--kmax_by_mutype']
            if 'recombination' in to_hash['parameters']:
                del to_hash['parameters']['recombination']
        elif module == 'simulate':
            pass
        else:
            ValueError("Not implemented yet.")
        del to_hash['population_by_letter']
        for pop_name in ['A', 'B']:
            del to_hash['populations'][pop_name]
        del to_hash['gimble']
        if return_dict:
            return (hashlib.md5(str(to_hash).encode()).hexdigest(), to_hash)
        return hashlib.md5(str(to_hash).encode()).hexdigest()

    def _expand_params(self, remove=None):
        if len(self.config['parameters'])>0:
            parameter_combinations = collections.defaultdict(list)
            if remove is not None:
                for key in remove:
                    del self.config['parameters'][f'Ne_{key}']
            for key, value in self.config['parameters'].items():
                if isinstance(value, float) or isinstance(value, int):
                    parameter_combinations[key]=np.array([value,], dtype=np.float64)
                elif key=='recombination':
                    pass
                else:
                    if len(value) == 4:
                        if self._MODULE == 'optimize':
                            sys.exit(f"[X] {self._MODULE} requires a single point or boundary for all parameters.")
                        minv, maxv, n, scale = value
                        sim_range = self._expand_params_scale(scale, minv, maxv, n)
                        parameter_combinations[key] = sim_range
                    elif len(value) <= 2:
                        parameter_combinations[key] = np.unique(value)
                    else:
                        raise ValueError("Uncaught error in config file configuration.")
            return parameter_combinations
        else:
            raise ValueError("config parameters does not contain any parameters.")

    def _expand_params_scale(self, scale, minv, maxv, n):
        if scale.startswith('lin'):
            return np.linspace(minv, maxv, num=n, endpoint=True, dtype=np.float64)
        elif scale.startswith('log'):
            return np.logspace(minv, maxv, num=n, endpoint=True, dtype=np.float64)
        else:
            sys.exit("scale in config parameters should either be lin or log")
"""
    def _make_parameter_combinations(self, sync_reference=None, sync_target=None):
        parameter_combinations = self._expand_params(remove=sync_target)
        if self._MODULE == 'optimize':
            return parameter_combinations
        if self._MODULE == 'simulate':
            rec = self._set_recombination_rate()
            if rec != None:
                parameter_combinations['recombination'] = rec
        parameter_combinations =  self._dict_product(parameter_combinations)
        if sync_reference and sync_target:
            for pop in sync_target:
                parameter_combinations[f'Ne_{pop}'] = parameter_combinations[f'Ne_{sync_reference}']
        return parameter_combinations
"""
#functions for gridsearch/sim_grid

def get_slice_grid_meta_idxs(grid_meta_dict=None, lncls=None, fixed_parameter=None, parameter_value=None):
    '''
    fixed_parameter=None, parameter_value=None  => 1 idx (overall max likelihood gridkey)
    fixed_parameter=str, parameter_value=None   => list of n 1d-arrays of idxs with shape (windows,) (n=unique values of parameter)
    fixed_parameter=str, parameter_value=float  => 1d-array of idxs with shape (windows,) 
    '''
    if fixed_parameter: 
        if isinstance(grid_meta_dict, list):
            values_by_parameter = LOD_to_DOL(grid_meta_dict)
            #values_by_parameter = grid_meta_dict_to_value_arrays_by_parameter(grid_meta_dict)
        else:
            values_by_parameter = {k:np.array(v, dtype=np.float64) for k,v in grid_meta_dict.items()}
        if not fixed_parameter in values_by_parameter:
            raise ValueError("%r is not part of this model" % fixed_parameter)
        fixed_parameter_values = np.array(values_by_parameter[fixed_parameter])
        if parameter_value:
            fixed_parameter_indices = np.concatenate(np.argwhere(parameter_value==fixed_parameter_values))
            #if not np.any(fixed_parameter_indices): #will evaluate to True for np.array([0])
            if fixed_parameter_indices.size==0:
                raise ValueError("parameter_value %r not found in grid" % parameter_value)
            fixed_parameter_lncls = lncls[:, fixed_parameter_indices]
            fixed_parameter_lncls_max_idx = np.argmax(fixed_parameter_lncls, axis=1)        
            fixed_parameter_lncls_max_idx = fixed_parameter_indices[fixed_parameter_lncls_max_idx] #translate back to larger grid_meta_dict idxs
            return fixed_parameter_lncls_max_idx
        results = []
        for i in np.unique(fixed_parameter_values): #these values are sorted
            fixed_parameter_indices = np.concatenate(np.argwhere(i==fixed_parameter_values))
            fixed_parameter_lncls = lncls[:, fixed_parameter_indices]
            fixed_parameter_lncls_max_idx = np.argmax(fixed_parameter_lncls, axis=1)
            idxs = fixed_parameter_indices[fixed_parameter_lncls_max_idx]
            results.append(idxs)
        return results
    return np.argmax(lncls, axis=1)

def _get_sim_grid_config(config, lncls_global, lncls_windows, meta, window_info=None, fixed_param_grid=None):
    rec, rec_map = None, None
    if config['simulate']['recombination_map'].strip()!='':
        store_window_df = pd.DataFrame(window_info).T
        store_window_df.columns = ['chr', 'start', 'end', 'idx']
        rbins = config["simulate"]["number_bins"]
        cutoff = config["simulate"]["cutoff"]
        scale = config["simulate"]["scale"]
        rbins = 10 if rbins=='' else rbins
        scale = 'lin' if scale=='' else scale
        cutoff = 90 if cutoff=='' else cutoff
        #recmap with seq, start, end, rec_bin_value
        rec_map = _parse_recombination_map(config['simulate']['recombination_map'], cutoff, rbins, scale, store_window_df) 
    else:
        rec = config['simulate']['recombination_rate'] 
    grid_meta_dict = {k:np.array(v) for k,v in meta['grid_dict'].items()} #why are these not numpy arrays?
    grid_to_sim, r, window_param_idx = _get_sim_grid(grid_meta_dict, lncls_global, lncls_windows, fixed_param_grid, rec, rec_map)
    if isinstance(rec_map, pd.DataFrame):
        config['simulate']['recombination_rate'] = tuple(r)
    #ready to overwrite config: 
    config['parameters_expanded'] = grid_to_sim
    config['parameters_grid_points'] = len(next(iter(grid_to_sim.values())))
    config['window_param_idx'] = window_param_idx
    #do we remove entries in config that no longer make sense: parameters_np, ... ?
    return config

def _get_sim_grid(grid_meta_dict, lncls_global, lncls_windows, fixed_param_grid, rec=None, rec_map=None):
    #lncls global should be based on w_bsfs !
    global_winning_fixed_param_idx = np.argmax(lncls_global)
    #get optimal parametercombo given background for fixed parameter
    if fixed_param_grid:
        global_winning_fixed_param_value = grid_meta_dict[fixed_param_grid][global_winning_fixed_param_idx]
        local_winning_fixed_param_idx = get_slice_grid_meta_idxs(
            lncls=lncls_windows, grid_meta_dict=grid_meta_dict, 
            fixed_parameter=fixed_param_grid, 
            parameter_value=global_winning_fixed_param_value
            )
    else:
        global_winning_fixed_param_value = None
        local_winning_fixed_param_idx = get_slice_grid_meta_idxs(lncls=lncls_windows)
    
    if rec==None:
        assert(rec_map.shape[0]==len(local_winning_fixed_param_idx)), "Index recmap and windows not matching. Should have been caught."
        rec_bin_sorted = rec_map['rec_bins'].to_numpy(dtype=np.float64)
        grid_to_sim, r_to_sim, window_df = _get_sim_grid_with_rec_map(rec_bin_sorted, local_winning_fixed_param_idx, grid_meta_dict)
        return (grid_to_sim, r_to_sim, window_df)
    else:
        grid_to_sim, window_df = _get_sim_grid_fixed_rec(rec, local_winning_fixed_param_idx, grid_meta_dict)
        return (grid_to_sim, rec, window_df)
        
def _get_sim_grid_with_rec_map(rec_bin_sorted, local_winning_fixed_param_idx, grid_meta_dict):
    param_for_window = tuple(zip(local_winning_fixed_param_idx, rec_bin_sorted))
    unique, param_combo_idxs = np.unique(param_for_window, return_inverse=True, axis=0)
    r_to_sim = unique[:,1]
    unique_idxs = unique[:,0].astype(np.uint64)
    grid_to_sim = {k:grid_meta_dict[k][unique_idxs] for k in grid_meta_dict.keys()}
    return (grid_to_sim, r_to_sim, param_combo_idxs)

def _get_sim_grid_fixed_rec(rec_rate, local_winning_fixed_param_idx, grid_meta_dict):
    param_combo_idxs = np.unique(local_winning_fixed_param_idx)
    grid_to_sim = {k:grid_meta_dict[k][param_combo_idxs] for k in grid_meta_dict.keys()}
    old_to_new_idx = {old_idx:new_idx for new_idx, old_idx in enumerate(param_combo_idxs)}
    mapped_window_index = [old_to_new_idx[idx] for idx in local_winning_fixed_param_idx]
    return (grid_to_sim, mapped_window_index)

def _parse_recombination_map(path, cutoff, bins, scale, store_df):
    hapmap = pd.read_csv(path, sep='\t', 
        names=['chr', 'start', 'end', 'rec'], header=0)
    hapmap = _validate_recombination_map(store_df, hapmap)
    hapmap = hapmap.sort_values('idx')
    #from cM/Mb to rec/bp
    hapmap['rec_scaled'] = hapmap['rec']*1e-8
    return _make_bins(hapmap, scale, cutoff, bins)

def _make_bins(df, scale, cutoff=90, bins=10):
    clip_value = np.percentile(df['rec_scaled'], cutoff)
    df['rec_clipped'] = df['rec_scaled'].clip(lower=None, upper=clip_value)
    df['rec_clipped'].replace(0,np.nan,inplace=True)
    start, stop =  df['rec_clipped'].min(), df['rec_clipped'].max()
    #determine bins
    if scale.lower() == "log":
        start, stop = np.log10(start), np.log10(stop)
        bin_edges = np.logspace(start, stop, num=bins+1)
    elif scale.lower() == "lin": 
        bin_edges = np.linspace(start, stop, num=bins+1)
    else:
        sys.exit("[X] Scale of recombination values to be simulated should either be LINEAR or LOG")
    to_be_simulated  = [(bstop + bstart)/2 for bstart, bstop in zip(bin_edges[:-1],bin_edges[1:])]     
    df['rec_bins'] = pd.cut(df['rec_clipped'], bins, labels=to_be_simulated).astype(float)
    df['rec_bins'].replace(np.nan, 0.0, inplace=True)
    return df[['chr', 'start', 'end', 'rec_bins']]

def _validate_recombination_map(store_df, user_df):
    # check for consistency between window_coordinates and rec_map coordinates
    df_to_test = user_df[['chr', 'start', 'end']]
    df_merge = df_to_test.merge(store_df, how='outer', on=['chr', 'start', 'end'])
    if store_df.shape != df_merge.shape:
        sys.exit("[X] Recombination map coordinates do not match window coordinates. Use query to get window coordinates.")
    return store_df.merge(user_df, on=['chr', 'start', 'end'])

def _gridsearch_sims_single(data, grids, fixed_param_grid, gridded_params, grid_meta_dict, label, name, fixed_param_grid_value_idx, output_df=True):
    """
    -> transfered as function, should be redundant
    data=sim_bsfs, grids=lncls for each grid point, fixed_param_grid=fixed parameter, gridded_params=dict of parameters that are gridded,
    grid_meta_dict = , fixed_param_grid_value_idx=index of background value of fixed param
    determines optimal parameter combination (for each of theparameters that determine an axis along the grid)
    and returns list df containing those parameter combinations
    if a fixed_param_grid is passed: particular parameter along grid that is fixed to the global optimual value
        returns distribution of lncls for each value of that fixed parameter
    """
    assert np.product(data.shape[1:])==np.product(grids.shape[1:]), "Dimensions of sim bSFS and grid bSFS do not correspond. k_max does not correspond but not caught."
    data = np.reshape(data, (data.shape[0],-1)) #data shape: replicates * bSFS
    grids = np.reshape(grids, (grids.shape[0],-1)) #grids shape: num_grid_points * bSFS
    grids_log = np.zeros(grids.shape, dtype=np.float64)
    grids = np.log(grids, where=grids>0, out=grids_log)
    lncls = np.inner(data, grids_log) #result shape: replicates*num_grid_points
    param_idxs = np.argmax(lncls, axis=1) #returns optimal index for each replicate

    results_dict = {}
    for key in gridded_params:
        results_dict[key] = grid_meta_dict[key][param_idxs]
    df = pd.DataFrame(results_dict)
    if output_df:
        summary=df.describe(percentiles=[0.025, 0.05, 0.95, 0.975])
        summary.to_csv(f'{label}_{name}_summary.csv')
    df_fixed_param = None

    if fixed_param_grid:
        assert fixed_param_grid in gridded_params, "fixed param for bootstrap not in gridded_params list! Report this issue."
        columns = [] #shape = num_values_fixed_param * replicates
        #values are not sorted!
        for fixed_param_value_idxs in get_slice_grid_meta_idxs(grid_meta_dict=grid_meta_dict, lncls=lncls, fixed_parameter=fixed_param_grid):
            best_likelihoods = lncls[np.arange(lncls.shape[0]), fixed_param_value_idxs]
            columns.append(best_likelihoods)
        #results in column are sorted from smallest to largest param value
        columns= np.array(columns).T
        df_fixed_param = pd.DataFrame(columns, columns=[f"{fixed_param_grid}_{str(i)}" if i!=fixed_param_grid_value_idx else f'{fixed_param_grid}_background' for i in range(columns.shape[1])])
        if output_df:
            df_fixed_param.to_csv(f'{label}_{name}_lnCL_dist.csv')
    return (df, df_fixed_param)

def old_gridsearch_np(tally=None, grid=None):
    '''returns 2d array of likelihoods of shape (windows, grid)'''
    if grid is None or tally is None:
        return None
    tally = tally if isinstance(tally, np.ndarray) else np.array(tally)
    grid_log = np.zeros(grid.shape, dtype=GRIDSEARCH_DTYPE)
    print("[+] Initialising grid array ...")
    np.log(grid, where=grid>0, out=grid_log)
    print('grid.nbytes', grid.nbytes)
    if tally.ndim == 4:
        return np.squeeze(np.apply_over_axes(np.sum, (tally * grid_log), axes=[-4,-3,-2,-1]))
    return np.squeeze(np.apply_over_axes(np.sum, (tally[:, None] * grid_log), axes=[-4,-3,-2,-1]))

def gridsearch_np(tally=None, grid=None):
    '''returns 2d array of likelihoods of shape (windows, grid)'''
    if grid is None or tally is None:
        return None
    tally = tally if isinstance(tally, np.ndarray) else np.array(tally)
    if not tally.ndim == 4: # if tally.ndim == 5:
        tally = tally[:, None]
    #print('tally.shape', tally.shape)
    #print('tally.dtype', tally.dtype)
    #print('tally.nbytes', tally.nbytes)
    print("[+] Initialising grid array ...")
    grid = grid.astype(_return_np_type(grid))
    grid_log = np.zeros(grid.shape, dtype=GRIDSEARCH_DTYPE)
    np.log(grid, dtype=GRIDSEARCH_DTYPE, where=grid>0, out=grid_log)
    #print('grid.shape', grid.shape)
    #print('grid.nbytes', grid.nbytes)
    #print('grid.dtype', grid.dtype)
    product = tally * grid_log
    #print('product.shape', product.shape)
    #print('product.nbytes', product.nbytes)
    #print('product.dtype', product.dtype)
    return np.sum(product.reshape((product.shape[0], product.shape[1], np.prod(product.shape[2:]))), axis=-1)
    
def gridsearch_dask(tally=None, grid=None):
    '''returns 2d array of likelihoods of shape (windows, grid)'''
    if grid is None or tally is None:
        return None
    tally = tally if isinstance(tally, np.ndarray) else np.array(tally)
    if not tally.ndim == 4: # if tally.ndim == 5:
        tally = tally[:, None]
    grid = grid.astype(_return_np_type(grid))
    grid_log = np.zeros(grid.shape, dtype=_return_np_type(grid))
    np.log(grid, dtype=_return_np_type(grid_log), where=grid>0, out=grid_log)
    #from dask.distributed import Client
    #client = Client()
    #print(client)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        tally = dask.array.from_array(tally, chunks=(500, 1, tally.shape[-4], tally.shape[-3], tally.shape[-2], tally.shape[-1]))
        grid_log = dask.array.from_array(grid_log, chunks=(500, grid_log.shape[-4], grid_log.shape[-3], grid_log.shape[-2], grid_log.shape[-1]))
        product = dask.array.multiply(tally, grid_log)
        result = dask.array.sum(product.reshape((product.shape[0], product.shape[1], np.prod(product.shape[2:]))), axis=-1)
        with daProgressBar():
            out = result.compute()
        return out

class Store(object):
    def __init__(self, prefix=None, path=None, create=False, overwrite=False):
        self.prefix = prefix if not prefix is None else str(pathlib.Path(path).resolve().stem)
        self.path = path if not path is None else "%s.z" % prefix
        self.data = self._init_store(create, overwrite)

    def tree(self):
        print(self.data.tree())
    
    def log_action(self, module, command):
        self.data.attrs[module] = command
    
    def get_stage(self, stage):
        return self.data.attrs[stage]

    def has_stage(self, stage):
        return stage in self.data.attrs

    def setup_sim(self, parameterObj):
        print("[#] Preparing store...")
        self._init_meta(overwrite=True)

    def measure(self, genome_f=None, sample_f=None, bed_f=None, vcf_f=None):
        #measure_key = self._get_key(task='measure')
        #self._set_meta(measure_key)
        measure_key = "seqs/"
        self._set_meta(measure_key)
        print("[#] Processing GENOME_FILE %r." % genome_f)
        self._read_sequences(measure_key, genome_f)
        print("[#] Processing SAMPLE_FILE %r." % sample_f)
        self._read_samples(measure_key, sample_f)
        print("[#] Processing BED_FILE %r." % bed_f)
        self._read_intervals(measure_key, bed_f)
        print("[#] Processing VCF_FILE %r." % vcf_f)
        self._read_variants(measure_key, vcf_f)
        #print(self.data.tree())

    def _read_sequences(self, measure_key, genome_f):
        sequences_df = parse_csv(
            csv_f=genome_f, 
            sep="\t", 
            usecols=[0,1], 
            dtype={'sequence_id': 'category', 'sequence_length': 'int64'}, 
            header=None)
        #meta = self._get_meta(measure_key)
        meta = self._get_meta('seqs')
        meta['seq_names'] = sequences_df['sequence_id'].to_list()
        meta['seq_lengths'] = sequences_df['sequence_length'].to_list()
        meta['seq_n50'] = get_n50_from_lengths(meta['seq_lengths'])
        meta['genome_f'] = genome_f

    def _read_samples(self, measure_key, sample_f):
        samples_df = parse_csv(
            csv_f=sample_f, 
            sep=",", 
            usecols=[0,1], 
            dtype={'samples': 'object', 'populations': 'category'}, 
            header=None)
        #meta = self._get_meta(measure_key)
        meta = self._get_meta('seqs')
        meta['samples'] = samples_df['samples'].to_list()
        meta['populations'] = samples_df['populations'].to_list()
        meta['population_ids'] = sorted(set(samples_df['populations'].to_list()))
        meta['population_by_letter'] = {letter: population_id for population_id, letter in zip(meta['population_ids'], string.ascii_uppercase)}
        meta['population_by_sample'] = {sample: population for sample, population in zip(meta['samples'], meta['populations'])}
        meta['sample_sets'] = [
            tuple(sorted(x, key=(meta['population_by_sample'].get if meta['population_by_sample'][x[0]] != meta['population_by_sample'][x[1]] else None))) 
                for x in itertools.combinations(meta['population_by_sample'].keys(), 2)]
        meta['sample_sets_inter'] = [
            False if len(set([meta['population_by_sample'][sample] for sample in sample_set])) == 1 else True 
                for sample_set in meta['sample_sets']]
        meta['sample_sets_intra_A'] = [
            all([meta['population_by_sample'][sample] == meta['population_ids'][0] for sample in sample_set]) for sample_set in meta['sample_sets']]
        meta['sample_sets_intra_B'] = [
            all([meta['population_by_sample'][sample] == meta['population_ids'][1] for sample in sample_set]) for sample_set in meta['sample_sets']]
        meta['sample_f'] = sample_f
        # ---> reports
        #longest_sample_string = max([len(", ".join(sample_set)) for sample_set in meta['sample_sets']]) + 2
        #meta['spacing'] = longest_sample_string if longest_sample_string > meta['spacing'] else meta['spacing']

    def _read_variants(self, measure_key, vcf_f):
        meta = self._get_meta(measure_key)
        seq_names = meta['seq_names']
        samples = meta['samples']
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            gt_key, pos_key, sample_key = 'calldata/GT', 'variants/POS', 'samples'
            samples_gt_order = allel.read_vcf(vcf_f, fields=[sample_key])[sample_key]
            query_samples = ordered_intersect(a=samples_gt_order, b=samples, order='a')
            # Check if all samples were found
            if set(query_samples) != set(samples):
                sys.exit("[X] The following samples in SAMPLE_FILE were not found in VCF_FILE: %s" % (
                    ", ".join(list(set(samples).difference(set(query_samples))))
                    ))
            # Set up counts arrays
            count_shape = (len(meta['seq_names']), len(query_samples))
            count_records = np.zeros(count_shape[0], dtype=np.int64)
            count_called = np.zeros(count_shape, dtype=np.int64)
            count_hom_ref = np.zeros(count_shape, dtype=np.int64)
            count_hom_alt = np.zeros(count_shape, dtype=np.int64)
            count_het = np.zeros(count_shape, dtype=np.int64)
            count_missing = np.zeros(count_shape, dtype=np.int64)
            for idx, seq_name in tqdm(enumerate(seq_names), total=len(seq_names), desc="[%] Reading variants", ncols=100):
                vcf_data = allel.read_vcf(vcf_f, region=seq_name, samples=query_samples, fields=[gt_key, pos_key])
                if vcf_data:
                    # genotypes
                    gt_matrix_raw = vcf_data[gt_key]
                    # counts
                    intervals = self._get_interval_coordinates(seq_name=seq_name)
                    sites = intervals_to_sites(intervals)
                    if sites is not None:
                        # positions in VCF
                        pos_array_raw = check_unique_pos((vcf_data[pos_key] - 1)) # port to BED (0-based) coordinates
                        # intersection of VCF and BED intervals
                        interval_mask = np.isin(pos_array_raw, sites, assume_unique=True)
                        gt_matrix = gt_matrix_raw[interval_mask]
                        count_records[idx] = gt_matrix.shape[0]
                        pos_array = pos_array_raw[interval_mask]
                        sa_genotype_matrix = allel.GenotypeArray(gt_matrix)
                        count_called[idx,:] = sa_genotype_matrix.count_called(axis=0)
                        count_hom_ref[idx,:] = sa_genotype_matrix.count_hom_ref(axis=0)
                        count_hom_alt[idx,:] = sa_genotype_matrix.count_hom_alt(axis=0)
                        count_het[idx,:] = sa_genotype_matrix.count_het(axis=0)
                        count_missing[idx,:] = sa_genotype_matrix.count_missing(axis=0)
                        self._save_variants(seq_name, pos_array, gt_matrix)        
            meta['variants_idx_by_sample'] = {query_sample: idx for idx, query_sample in enumerate(query_samples)}
        meta['vcf_f'] = vcf_f
        meta['variants_counts'] = int(np.sum(count_records)) # ZARR JSON encoder does not like numpy dtypes
        meta['variants_counts_called'] = [int(x) for x in np.sum(count_called, axis=0)] # ZARR JSON encoder does not like numpy dtypes
        meta['variants_counts_hom_ref'] = [int(x) for x in np.sum(count_hom_ref, axis=0)] # ZARR JSON encoder does not like numpy dtypes
        meta['variants_counts_hom_alt'] = [int(x) for x in np.sum(count_hom_alt, axis=0)] # ZARR JSON encoder does not like numpy dtypes
        meta['variants_counts_het'] = [int(x) for x in np.sum(count_het, axis=0)] # ZARR JSON encoder does not like numpy dtypes
        meta['variants_counts_missing'] = [int(x) for x in np.sum(count_missing, axis=0)] # ZARR JSON encoder does not like numpy dtypes
        # QC plots 

    def _save_variants(self, sequence, pos_array, gt_matrix):
        self.data.create_dataset("seqs/%s/variants/pos" % sequence, data=pos_array, dtype=np.int64)
        self.data.create_dataset("seqs/%s/variants/matrix" % sequence, data=gt_matrix, dtype=np.int64)

    def _save_variants_meta(self):
        pass

    def _read_intervals(self, measure_key, bed_f):
        meta = self._get_meta(measure_key)
        target_sequences, target_samples = set(meta['seq_names']), set(meta['samples'])
        intervals_df = parse_intervals(bed_f, target_sequences, target_samples)
        valid_sequences = intervals_df['sequence'].unique() 
        intervals_idx_by_sample = {sample: idx for idx, sample in enumerate(intervals_df.columns[3:-1])}
        intervals_count = len(intervals_df.index)
        intervals_span = int(intervals_df['length'].sum())
        count_bases_samples = np.zeros((len(valid_sequences), len(target_samples)), dtype=np.int64)
        for idx, (sequence, _df) in tqdm(enumerate(intervals_df.groupby(['sequence'], observed=True)), total=len(valid_sequences), desc="[%] Reading intervals", ncols=100):
            interval_matrix = _df[target_samples].to_numpy()
            starts = _df['start'].to_numpy()
            ends = _df['end'].to_numpy()
            #interval_matrix_key = self._get_key(task='measure', data_label=intervals, seq_label=sequence,)
            self._set_intervals(sequence, interval_matrix, starts, ends)
            length_matrix = interval_matrix * _df['length'].to_numpy().reshape(-1, 1)
            count_bases_samples[idx,:] = np.sum(length_matrix, axis=0)
        self._set_intervals_meta(bed_f, intervals_idx_by_sample, count_bases_samples, intervals_count, intervals_span)

    def _set_intervals(self, sequence, interval_matrix, starts, ends):
        self.data.create_dataset("seqs/%s/intervals/matrix" % sequence, data=interval_matrix)
        self.data.create_dataset("seqs/%s/intervals/starts" % sequence, data=starts)
        self.data.create_dataset("seqs/%s/intervals/ends" % sequence, data=ends)

    def _set_intervals_meta(self, bed_f, intervals_idx_by_sample, count_bases_samples, intervals_count, intervals_span):
        meta_intervals = self._get_meta('seqs')
        meta_intervals['bed_f'] = bed_f
        meta_intervals['intervals_idx_by_sample'] = intervals_idx_by_sample
        meta_intervals['intervals_span_sample'] = [int(x) for x in np.sum(count_bases_samples, axis=0)] # JSON encoder does not like numpy dtypes   
        meta_intervals['intervals_count'] = intervals_count
        meta_intervals['intervals_span'] = intervals_span

    def blocks(self, block_length=64, block_span=128, block_max_multiallelic=3, block_max_missing=3, overwrite=False):
        self._preflight_blocks(overwrite=overwrite)
        print("[+] Blocking parameters = [-l %s -m %s -u %s -i %s]" % (
            block_length, 
            block_span, 
            block_max_multiallelic, 
            block_max_missing))
        self._make_blocks(block_length, block_span, block_max_missing, block_max_multiallelic)
        # after making blocks, make tally 'raw' without kmax. This is needed for heterozygosity, d_xy, F_sts metrics, etc
        self.tally('blocks', 'blocks_raw', None, 'X', None, None, overwrite=True, verbose=False)

    def windows(self, window_size=500, window_step=100, overwrite=False):
        self._preflight_windows(overwrite)
        print("[+] Window parameters = [-w %s -s %s]" % (window_size, window_step))
        self._make_windows(window_size, window_step, sample_sets='X')
        self.tally('windows', 'windows_raw', None, 'X', None, None, overwrite=True, verbose=False)
        self.tally('windows', 'windowsum_raw', None, 'X', None, None, overwrite=True, verbose=False)

    def _preflight_simulate(self, config, overwrite):
        config['simulate_key'] = self._get_key(task='simulate', analysis_label=config['gimble']['label'])
        if not overwrite and self._has_key(config['simulate_key']):
            sys.exit("[X] Simulated results with label %r already exist. Use '-f' to overwrite or change analysis label in the INI file." % 
                (config['simulate_key']))
        #deal with grid_label:
        if config['gridbased']['grid_label'].strip()!='':
            #check whether grid_label exists
            key_4D = self._get_key(task='gridsearch', data_label='windows', analysis_label=config['gridbased']['grid_label'], mod_label='4D')
            key_5D = self._get_key(task='gridsearch', data_label='windows', analysis_label=config['gridbased']['grid_label'], mod_label='5D')
            if not self._has_key(key_5D):
                sys.exit("[X] Provided grid label does not correspond to an existing grid. Run gridsearch first.")
            else:
                if config['simulate']['recombination_map']!='':
                    window_info = self._get_window_bed_columns()
                else:
                    window_info = None
                meta_5D = self._get_meta(key_5D)
                lncls_global = self._get_data(key_4D)
                lncls_windows = self._get_data(key_5D)
                fixed_param_grid = config['gridbased']['fixed_parameter'] if config['gridbased']['fixed_parameter']!='' else None
                config = _get_sim_grid_config(config, lncls_global, lncls_windows, meta_5D, window_info, fixed_param_grid)
        #print('[+] Simulating %s replicate(s) of %s block(s) for %s parameter combinations' %
        #    (config['simulate']['replicates'], config['simulate']['blocks'], config['parameters_grid_points']))
        return config

    def simulate(self, config, threads, overwrite):
        print("[#] Preflight...")
        config = self._preflight_simulate(config, overwrite)
        print('[+] Simulating %s replicate(s) of %s block(s) for %s parameter combinations' %
            (config['simulate']['replicates'], config['simulate']['blocks'], config['parameters_grid_points']))
        for idx, simulation_instance in enumerate(
            lib.simulate.run_sims(
                config, 
                threads,
                discrete=config['simulate']['discrete_genome']
                )
            ):
            self._save_simulation_instance(config, simulation_instance, idx)
        self._save_simulation_meta(config)
        #we need particular subset of details to be saved
        #self.data[f"sims/{config['gimble']['label']}"].attrs.put(global_info_ancestry)
        
        # needed for --grid flag check this
        # gimbleStore.data[f'sims/{group_name}'].attrs['fixed_param_grid'] = parameterObj.fixed_param_grid
        
        # writing report (has to be fixed/should be part of info)
        # try:
        #     print(self._get_sims_report(width=100, label=config['gimble']['label']))   
        # except UnicodeEncodeError:
        #     print(self._get_sims_report(width=100, label=config['gimble']['label']).__repr__().encode('utf-8')) #temp fix for my environment

    def _save_simulation_meta(self, config):
        meta = config_to_meta(config, 'simulate')
        self._set_meta(config['simulate_key'], meta)

    def _save_simulation_instance(self, config, simulation_instance, idx):
        simulation_instance_key = self._get_key(task='simulate', analysis_label=config['gimble']['label'], parameter_label=idx)
        config['idx'] = idx
        meta = config_to_meta(config, 'simulate_instance')
        self._set_meta_and_data(simulation_instance_key, meta, simulation_instance)

    def _preflight_query(self, data_key, version):
        if not data_key:
            '''Still needs checking whether keys are found correctly for all modules
            # blocks √
            # windows √
            # tally √ 
            # optimize √
            # makegrid ?
            # gridsearch ?
            # simulate ?
            '''
            available_keys_by_category = collections.defaultdict(set)
            max_depth_by_key = {
                'blocks': 1, 'windows': 1,
                'tally': 2, 'makegrid': 2, 'optimize': 3, 
                'gridsearch': 3, 'simulate': 3}
            category_by_module = {
                'blocks': 'measure',
                'windows': 'measure',
                'tally': 'tally',
                'makegrid': 'makegrid',
                'simulate': 'simulate',
                'optimize': 'optimize',
                'gridsearch': 'gridsearch'
                }
            def data_key_finder(path):
                key_list = str(path).split("/")
                module = key_list[0]
                if len(key_list) == max_depth_by_key.get(module, None):
                    available_keys_by_category[category_by_module[module]].add("/".join(key_list[0:max_depth_by_key.get(key_list[0], 0)]))
            self.data.visit(data_key_finder)
            print("[X] Please specify a label (-l) for which data to query. Available labels:")
            for category, available_keys in available_keys_by_category.items():
                print("# %s" % category)
                print("- %s" % "\n- ".join(sorted(available_keys)))
            sys.exit(1)
        config = {'data_key': data_key, 'version': version}
        if not self._has_key(data_key):
            sys.exit("[X] ZARR store %s has no data under the key %r." % (self.path, data_key))
        config['data_type'] = data_key.split("/")[0]
        return config

    def query(self, version, data_key):
        config = self._preflight_query(data_key, version)
        if config['data_type'] == 'tally':
            self._write_tally_tsv(config)
        elif config['data_type'] == 'optimize':
            self._write_optimize_tsv(config)
        elif config['data_type'] == 'windows':
            #sys.exit("[X] Not implemented.")
            self._write_window_bed(config)
        elif config['data_type'] == 'gridsearch':
            self._write_gridsearch_bed(config)
        elif config['data_type'] == 'makegrid':
            meta = self._get_meta(config['data_key'])
            print(dict(meta))
        else:
            sys.exit("[X] Not implemented.")

    def _write_optimize_tsv(self, config):
        optimize_meta = dict(self._get_meta(config['data_key']))
        # prints optimize_meta, could have prettier formatting ...
        print("[+] Optimize ...")
        query_meta = format_query_meta(optimize_meta)
        print(query_meta)
        single_file_flag = (len(optimize_meta['optimize_keys']) == 1)
        for idx, optimize_key in enumerate(optimize_meta['optimize_keys']):
            instance_meta = dict(self._get_meta(optimize_key))
            # optima = pd.DataFrame(instance_meta['optimize_results'])
            if single_file_flag:
                fn = "%s.%s.optimize.tsv" % (self.prefix, config['data_key'].replace("/", "_"))
            else:
                fn = "%s.%s.optimize.%s.tsv" % (self.prefix, config['data_key'].replace("/", "_"), idx)
            pd.DataFrame(
                data=instance_meta['optimize_results']).to_csv(fn, index=True, index_label='idx', sep='\t')
            print("[#] Wrote file %r." % fn)

    def _write_tally_tsv(self, config):
        config['data'] = np.array(self._get_data(config['data_key']))
        config['meta'] = dict(self._get_meta(config['data_key']))
        print("[+] Tally ...")
        query_meta = format_query_meta(config['meta'])
        print(query_meta)
        config['header'] = ['count', 'm_1', 'm_2', 'm_3', 'm_4'] if config['data'].ndim == 4 else ['window_idx', 'count', 'm_1', 'm_2', 'm_3', 'm_4']
        config['filename'] = "%s.%s.tsv" % (self.prefix, config['data_key'].replace("/", "_"))
        pd.DataFrame(
            data=bsfs_to_2d(config['data']), 
            columns=config['header'],
            dtype='int64').to_csv(config['filename'], index=False, sep='\t')
        print("[#] Wrote file %r." % config['filename'])

    def _write_window_bed(self, config):
        meta_blocks = self._get_meta('blocks')
        meta_windows = self._get_meta('windows')
        print("[+] Windows ...")
        query_meta_blocks = format_query_meta(meta_blocks, ignore_long=True)
        print(query_meta_blocks)
        query_meta_windows = format_query_meta(meta_blocks, ignore_long=True)
        print(query_meta_windows)
        #print("[+] Getting data for BED ...")
        sequences, starts, ends, index, pos_mean, pos_median, balance, mse_sample_set_cov = self._get_window_bed()
        #print("[+] Calculating popgen metrics ...")
        tally = tally_variation(self._get_variation(data_type='windows', sample_sets='X', progress=False), form='tally')
        pop_metrics = get_popgen_metrics(tally, sites=(meta_blocks['length'] * meta_windows['size']))
        int_bed = np.vstack([starts, ends, index, pos_mean, pos_median, pop_metrics, balance, mse_sample_set_cov]).T
        columns = ['sequence', 'start', 'end', 'index', 'pos_mean', 'pos_median', 'heterozygosity_A', 'heterozygosity_B', 'd_xy', 'f_st', 'balance', 'mse_sample_set_cov']
        header = ["# %s" % config['version'], "# %s" % "\t".join(columns)]
        out_f = '%s.windows.popgen.bed' % self.prefix
        print("[+] Writing BED file %s ..." % out_f)
        with open(out_f, 'w') as out_fh:
            out_fh.write("\n".join(header) + "\n")
        dtypes = {
            'start': 'int64', 'end': 'int64', 'index': 'int64', 'pos_mean': 'int64', 
            'pos_median': 'int64','heterozygosity_A': 'float64', 'heterozygosity_B': 'float64', 
            'd_xy': 'float64', 'f_st': 'float64', 'balance': 'float64'}
        bed_df = pd.DataFrame(data=int_bed, columns=columns[1:]).astype(dtype=dtypes)
        bed_df['sequence'] = sequences
        bed_df.sort_values(['sequence', 'start'], ascending=[True, True]).to_csv(out_f, na_rep='NA', mode='a', sep='\t', index=False, header=False, columns=columns, float_format='%.5f')
        return out_f

    def dump_lncls(self, parameterObj):
        unique_hash = parameterObj._get_unique_hash()
        grids, grid_meta_dict = self._get_grid(unique_hash)
        lncls_global, lncls_windows = self._get_lncls(unique_hash)
        if isinstance(grid_meta_dict, list):
            values_by_parameter = LOD_to_DOL(grid_meta_dict)
            #values_by_parameter = grid_meta_dict_to_value_arrays_by_parameter(grid_meta_dict)
        else:
            values_by_parameter = grid_meta_dict #once everything works, values_by_parameter should be simply renamed
        sequences, starts, ends, index = self._get_window_bed_columns()
        parameter_names = [name for name in values_by_parameter.keys() if name != 'mu']
        grid_meta_dict.pop('mu', None) #remove mu from grid_meta_dict
        column_headers = ['sequence', 'start', 'end', 'index', 'lnCL'] + parameter_names + ['fixed']
        dtypes = {'start': 'int64', 'end': 'int64', 'index': 'int64', 'lnCL': 'float64'}
        for param in parameter_names:
            dtypes[param] = 'float64'
        MAX_PARAM_LENGTH = max([len(param) for param in parameter_names])
        for parameter in tqdm(parameter_names, total=len(parameter_names), desc="[%] Writing output...", ncols=100, unit_scale=True): 
            bed_dfs = []
            out_f = '%s.%s.gridsearch.lnCls.%s_fixed.tsv' % (self.prefix, parameterObj.data_type, parameter)
            fixed = np.full_like(lncls_windows.shape[0], parameter, dtype='<U%s' % MAX_PARAM_LENGTH)
            for grid_meta_idxs in get_slice_grid_meta_idxs(grid_meta_dict=grid_meta_dict, lncls=lncls_windows, fixed_parameter=parameter, parameter_value=None):
                best_likelihoods = lncls_windows[np.arange(lncls_windows.shape[0]), grid_meta_idxs]
                best_params = np.array(list(zip(*grid_meta_dict.values())))[grid_meta_idxs]    
                #line above replaces code until best_params = 
                #meta_dicts = list(np.vectorize(grid_meta_dict.__getitem__)(grid_meta_idxs.astype(str)))
                #columns = []
                #for param in parameter_names:
                #    column = []
                #    for meta_dict in meta_dicts:
                #        column.append(meta_dict[param])
                #    columns.append(column)
                #best_params = np.vstack(columns).T
                int_bed = np.vstack([starts, ends, index, best_likelihoods, best_params.T]).T
                header = ["# %s" % parameterObj._VERSION]
                header += ["# %s" % "\t".join(column_headers)]
                with open(out_f, 'w') as out_fh:
                    out_fh.write("\n".join(header) + "\n")
                bed_df = pd.DataFrame(data=int_bed, columns=column_headers[1:-1]).astype(dtype=dtypes)
                bed_df['sequence'] = sequences
                bed_df['fixed'] = fixed
                # MUST be mode='a' otherwise header gets wiped ...
                bed_dfs.append(bed_df)
            df = pd.concat(bed_dfs, ignore_index=True, axis=0)
            df.sort_values(['index'], ascending=[True]).to_csv(out_f, na_rep='NA', mode='a', sep='\t', index=False, header=False, columns=column_headers)

    def _get_window_bed_columns(self):
        meta_seqs = self._get_meta('seqs')
        meta_windows = self._get_meta('windows')
        MAX_SEQNAME_LENGTH = max([len(seq_name) for seq_name in meta_seqs['seq_names']])
        sequences = np.zeros(meta_windows['count'], dtype='<U%s' % MAX_SEQNAME_LENGTH)
        starts = np.zeros(meta_windows['count'], dtype=np.int64)
        ends = np.zeros(meta_windows['count'], dtype=np.int64)
        offset = 0
        for seq_name in tqdm(meta_seqs['seq_names'], total=len(meta_seqs['seq_names']), desc="[%] Preparing output...", ncols=100, unit_scale=True): 
            start_key = 'windows/%s/starts' % (seq_name)
            end_key = 'windows/%s/ends' % (seq_name)
            if start_key in self.data:
                start_array = np.array(self.data[start_key])
                window_count = start_array.shape[0]
                starts[offset:offset+window_count] = start_array
                ends[offset:offset+window_count] = np.array(self.data[end_key])
                sequences[offset:offset+window_count] = np.full_like(window_count, seq_name, dtype='<U%s' % MAX_SEQNAME_LENGTH)
                offset += window_count
        index = np.arange(meta_windows['count'],  dtype=np.int64)
        return (sequences, starts, ends, index)

    def _write_gridsearch_bed(self, config):
        meta_gridsearch = self._get_meta(config['data_key'])
        grids = DOL_to_LOD(meta_gridsearch['grid_dict'])
        grid_params = np.array([list(subdict.values()) for subdict in grids] ,dtype=np.float64)
        sequences, starts, ends, index, pos_mean, pos_median, balance, mse_sample_set_cov = self._get_window_bed()
        params_header = list(meta_gridsearch['grid_dict'].keys())
        out_fs = []
        for gridsearch_key in meta_gridsearch['gridsearch_keys']:
            lncls = np.array(self._get_data(gridsearch_key))
            #print('lncls', lncls)
            best_lncl = np.max(lncls, axis=lncls.ndim-1)
            #print('best_lncl', best_lncl)
            best_idx = np.argmax(lncls, axis=lncls.ndim-1)
            #print('best_idx', best_idx)
            best_params = grid_params[best_idx, :]
            #print('best_params', best_params)
            if lncls.ndim == 1:
                columns = ['lnCL'] + params_header
                dtypes = {column: 'float64' for column in columns}
                int_bed = [best_lncl] + list(best_params)
            else:
                columns = ['sequence', 'start', 'end', 'index', 'pos_mean', 'pos_median', 'balance', 'mse_sample_set_cov', 'lnCL'] + params_header
                int_bed = np.vstack([starts, ends, index, pos_mean, pos_median, balance, mse_sample_set_cov, best_lncl, best_params.T]).T
                dtypes = {'start': 'int64', 'end': 'int64', 'index': 'int64', 'pos_mean': 'float64', 'pos_median': 'float64', 
                    'balance': 'float64', 'mse_sample_set_cov': 'float64', 'lnCL': 'float64'}
                for param in params_header:
                    dtypes[param] = 'float64'
            out_f = '%s.bed' % ("_".join(gridsearch_key.split("/")))
            print("[+] Sum of lnCL for winning parameters = %s" % np.sum(best_lncl))
            # write header
            header = ["# %s" % config['version']]
            header += ["# %s" % "\t".join(columns)]
            with open(out_f, 'w') as out_fh:
                out_fh.write("\n".join(header) + "\n")
            if lncls.ndim > 1:
                bed_df = pd.DataFrame(data=int_bed, columns=columns[1:]).astype(dtype=dtypes)
                bed_df['sequence'] = sequences
                # MUST be mode='a' otherwise header gets wiped ...
                bed_df.sort_values(['sequence', 'start'], ascending=[True, True])
            else:
                bed_df = pd.DataFrame(data=int_bed, columns=columns).astype(dtype=dtypes)
            bed_df.to_csv(out_f, na_rep='NA', mode='a', sep='\t', index=False, header=False, columns=columns)
            out_fs.append(out_f)
        return out_fs
        #grids = DOL_to_LOD(meta_gridsearch['grid_dict'])
        ##for grid_idx, grid_dict in grid_meta_dict.items():
        ##    grids.append(list(grid_dict.values()))
        #grid_params = np.array([list(subdict.values()) for subdict in grids] ,dtype=np.float64)
        #best_params = grid_params[np.argmax(lncls, axis=1), :]
        #best_likelihoods = np.max(lncls, axis=1)
        #delta_lncls = best_likelihoods - lncls[:, best_idx]
        #sequences, starts, ends, index = self._get_window_bed_columns() 
        #params_header = list(next(iter(grids)).keys())
        ##params_header = list(grid_dict.keys())
        #meta_blocks = self._get_meta('blocks')
        #meta_windows = self._get_meta('windows')
        #bsfs_windows_full = self.get_bsfs(
        #        data_type='windows', 
        #        population_by_letter=parameterObj.config['population_by_letter'], 
        #        sample_sets='X')
        #popgen_metrics = pop_metrics_from_bsfs(bsfs_windows_full, block_length=meta_blocks['length'], window_size=meta_windows['size'])
        #popgen_header = ['heterozygosity_A', 'heterozygosity_B', 'd_xy', 'f_st']
        #columns = ['sequence', 'start', 'end', 'index', 'lnCL', 'delta_lnCl'] + params_header + popgen_header
        #dtypes = {'start': 'int64', 'end': 'int64', 'index': 'int64', 'lnCL': 'float64', 'delta_lnCl': 'float64'}
        #for param in params_header + popgen_header:
        #    dtypes[param] = 'float64'
        #'''dtypes := "object", "int64", "float64", "bool", "datetime64", "timedelta", "category"'''
        #int_bed = np.vstack([starts, ends, index, best_likelihoods, delta_lncls, best_params.T, popgen_metrics]).T
        #out_f = '%s.%s.gridsearch.bestfit.bed' % (self.prefix, parameterObj.data_type)
        #print("[+] Sum of lnCL for winning parameters = %s" % np.sum(best_likelihoods))
        ## write header
        #header = ["# %s" % parameterObj._VERSION]
        #header += ["# %s" % "\t".join(columns)]
        #with open(out_f, 'w') as out_fh:
        #    out_fh.write("\n".join(header) + "\n")
        #bed_df = pd.DataFrame(data=int_bed, columns=columns[1:]).astype(dtype=dtypes)
        #bed_df['sequence'] = sequences
        ## MUST be mode='a' otherwise header gets wiped ...
        #bed_df.sort_values(['sequence', 'start'], ascending=[True, True]).to_csv(out_f, na_rep='NA', mode='a', sep='\t', index=False, header=False, columns=columns)
        #return out_f

    def gridsearch_preflight(self, tally_label, sim_label, grid_label, overwrite):
        config = {}
        config['gridsearch_time'] = 'None' # just so that it initialised for later
        # first deal with grid
        config['makegrid_key'] = self._get_key(task='makegrid', analysis_label=grid_label)
        grid_meta = self._get_meta(config['makegrid_key'])
        if grid_meta is None:
            sys.exit("[X] gimbleStore has no grid labelled %r." % config['makegrid_key'])
        grid = np.array(self._get_data(config['makegrid_key']), dtype=np.float64) # grid is likelihoods 
        config['grid_dict'] = grid_meta['grid_dict'] # grid_dict is params
        # Error if no data
        config['data_source'] = 'sims' if sim_label else 'meas'
        config['data_label'] = sim_label or tally_label
        if sim_label:
            config['data_key'] = self._get_key(task='simulate', analysis_label=config['data_label'])
        if tally_label:
            config['data_key'] = self._get_key(task='tally', data_label=config['data_label'])
        if not self._has_key(config['data_key']):
            sys.exit("[X] gimbleStore has no %r." % config['data_key'])
        config['gridsearch_key'] = self._get_key(task='gridsearch', data_label=config['data_label'], analysis_label=grid_label)
        if not overwrite and self._has_key(config['gridsearch_key']):
            sys.exit("[X] Gridsearch results with grid label %r on data %r already exist. Use '-f' to replace." % (grid_label, data_label))
        if config['data_source'] == 'meas':
            data = ((0, self._get_data(config['data_key'])) for _ in (0,))
            meta = self._get_meta(config['data_key'])
            config['block_length_data'] = meta['block_length']
            config['batch_sites'] = meta['block_length'] * meta['blocks']
            config['max_k'] = np.array(meta['max_k']) # INI values get overwritten by data ...
            config['gridsearch_keys'] = [self._get_key(task='gridsearch', data_label=config['data_label'], analysis_label=grid_meta['label'], parameter_label=idx) for idx in range(1)]
        else:
            data = self._get_sims_bsfs(config['data_key']) # data is an iterator across parameter combos
            meta = self._get_meta(config['data_key'])
            config['block_length_data'] = meta['parameters']['block_length']
            config['max_k'] = np.array(meta['max_k']) # INI values get overwritten by data ...
            config['gridsearch_keys'] = [self._get_key(task='gridsearch', data_label=config['data_label'], analysis_label=grid_meta['label'], parameter_label=idx) for idx in range(meta['max_idx'] + 1)]
        config['block_length_grid'] = grid_meta['block_length']
        # checking whether block_length in data and grid are compatible
        if not config['block_length_data'] == config['block_length_grid']:
            sys.exit("[X] Block lengths in data %r (%s) and grid %r (%s) are not compatible.." % (
                data_label, 
                config['block_length_data'],
                grid_label, 
                config['block_length_grid']
                ))
        return (config, data, grid)

    def gridsearch(self, tally_label, sim_label, grid_label, overwrite):
        print("[+] Gathering data for gridsearch ...")
        config, data, grid = self.gridsearch_preflight(tally_label, sim_label, grid_label, overwrite)
        print("[+] Performing global gridsearch on %r ..." % (config['data_key']))
        for key, (idx, tally) in zip(config['gridsearch_keys'], data):
            #gridsearch_instance_result = gridsearch_np(tally=tally, grid=grid)
            #print('gridsearch_instance_result.nbytes', gridsearch_instance_result.nbytes)
            #print('gridsearch_instance_result.dtype', gridsearch_instance_result.dtype)
            #print('np.max(gridsearch_instance_result)', np.max(gridsearch_instance_result))
            gridsearch_dask_result = gridsearch_dask(tally=tally, grid=grid)
            #print('gridsearch_dask_result.nbytes', gridsearch_dask_result.nbytes)
            #print('gridsearch_dask_result.dtype', gridsearch_dask_result.dtype)
            #print('np.max(gridsearch_dask_result)', np.max(gridsearch_dask_result))
            #old_gridsearch_instance_result = old_gridsearch_np(tally=tally, grid=grid)
            #print('old_gridsearch_instance_result.shape', old_gridsearch_instance_result.shape)
            #if np.array_equal(gridsearch_instance_result, gridsearch_dask_result):
            #    print("[+] correct")
            #else:
            #    print("[+] incorrect")
            self._set_data(key, gridsearch_dask_result)
        self._set_meta(config['gridsearch_key'], config_to_meta(config, 'gridsearch'))
        meta = self._get_meta(config['gridsearch_key'])
        #print(dict(meta))
        
    #def save_gridsearch(self, config, gridsearch_4D_result, gridsearch_5D_result):
    #    gridsearch_meta = config_to_meta(config, 'gridsearch')
    #    self._set_meta_and_data(config['gridsearch_key_4D'], gridsearch_meta, gridsearch_4D_result)
    #    if not gridsearch_5D_result is None:
    #        self._set_meta_and_data(config['gridsearch_key_5D'], gridsearch_meta, gridsearch_5D_result)
    #    print("[+] Saved gridsearch results.")

    def _gridsearch_sims(self, parameterObj, grids, grid_meta_dict):
        #check parameters that were fixed initially:
        #all_dicts = [grid_meta_dict[str(i)] for i in range(len(grid_meta_dict))] 
        ##this can be omitted once parameterObj.parameter_combinations is
        ##in the shape {Ne_A:[v1, v2, v3, v4, ...], Ne_B:[v1', v2' , ...], ...}
        #keys = list(all_dicts[0].keys())
        #key_all_values_dict = {}
        #for key in keys:
        #    key_all_values_dict[key] = np.array([d[key] for d in all_dicts], dtype=np.float64)
        gridded_params = sorted([key for key, items in grid_meta_dict.items() if len(set(items))>1])

        if 'fixed_param_grid' in self.data[f'sims/{parameterObj.label}'].attrs:
            fixed_param_grid = self.data[f'sims/{parameterObj.label}'].attrs['fixed_param_grid']
            fixed_param_grid_value = self.data[f'sims/{parameterObj.label}/parameter_combination_0'].attrs[fixed_param_grid]
            unique_values_fixed_param = np.unique(grid_meta_dict[fixed_param_grid])
            #unique_values_fixed_param = np.unique(key_all_values_dict[fixed_param_grid])            
            fixed_param_grid_value_idx = np.where(unique_values_fixed_param==fixed_param_grid_value)[0][0]
        else:
            fixed_param_grid = None
            fixed_param_grid_value = None
        
        param_combo_iterator = self.get_bsfs(
            data_type='simulate', 
            population_by_letter=parameterObj.config['populations'], 
            sample_sets="X", 
            kmax_by_mutype=parameterObj.config['k_max'],
            label=parameterObj.label
            ) #iterator over each parameter combination that was sim'ed
        num_param_combos = self.data[f'sims/{parameterObj.label}'].__len__()
        for name, data in tqdm(param_combo_iterator, desc='Processing all parameter combinations', total=num_param_combos,  ncols=100):
            df, df_fixed_param = _gridsearch_sims_single(data, grids, fixed_param_grid, gridded_params, grid_meta_dict, parameterObj.label, name, fixed_param_grid_value_idx)
        print(f"[+] Output written to {os.getcwd()}")    
        if fixed_param_grid:
            print("[+] Fixed param values:")
            print('\t'.join(f'{fixed_param_grid}_{i}' if i !=fixed_param_grid_value_idx else f'{fixed_param_grid}_background' for i in range(len(unique_values_fixed_param))))
            print('\t'.join("{:.3e}".format(value) for value in unique_values_fixed_param))

    def _preflight_optimize(self, config, sim_label, tally_label, num_cores, start_point, max_iterations, xtol_rel, ftol_rel, overwrite):
        config['num_cores'] = num_cores
        config['max_iterations'] = max_iterations
        config['xtol_rel'] = xtol_rel
        config['ftol_rel'] = ftol_rel
        config['optimize_time'] = 'None' # just so that it initialised for later
        # Error if no data
        config['data_source'] = 'sims' if sim_label else 'meas'
        config['data_label'] = sim_label or tally_label
        if sim_label:
            config['data_key'] = self._get_key(task='simulate', analysis_label=config['data_label'])
        if tally_label:
            config['data_key'] = self._get_key(task='tally', data_label=config['data_label'])
        if not self._has_key(config['data_key']):
            sys.exit("[X] gimbleStore has no %r." % config['data_key'])
        # Error if results clash
        config['optimize_key'] = self._get_key(task='optimize', data_label=config['data_label'], analysis_label=config['gimble']['label'])
        if not overwrite and self._has_key(config['optimize_key']):
            sys.exit("[X] Analysis with label %r on data %r already exist. Change the label in the config file or use '--force'" % (config['gimble']['label'], config['data_label']))
        # get data (this could be a separate function, also needed for gridsearch)
        if config['data_source'] == 'meas':
            data = ((0, self._get_data(config['data_key'])) for _ in (0,))
            meta = self._get_meta(config['data_key'])
            config['block_length'] = meta['block_length']
            config['max_k'] = np.array(meta['max_k']) # INI values get overwritten by data ...
            config['optimize_keys'] = [self._get_key(task='optimize', data_label=config['data_label'], analysis_label=config['gimble']['label'], parameter_label=idx) for idx in range(1)]
        else:
            data = self._get_sims_bsfs(config['data_key']) # data is an iterator across parameter combos
            meta = self._get_meta(config['data_key'])
            config['block_length'] = meta['parameters']['block_length']
            config['max_k'] = np.array(meta['max_k']) # INI values get overwritten by data ...
            config['optimize_keys'] = [self._get_key(task='optimize', data_label=config['data_label'], analysis_label=config['gimble']['label'], parameter_label=idx) for idx in range(meta['max_idx'] + 1)]
        # start point
        config['start_point_method'] = start_point
        if start_point == 'midpoint':
            config['start_point'] = np.mean(np.vstack((config['parameter_combinations_lowest'], config['parameter_combinations_highest'])), axis=0)
        if start_point == 'random':
            #np.random.seed(config['gimble']['random_seed'])
            config['start_point'] = np.random.uniform(low=config['parameter_combinations_lowest'], high=config['parameter_combinations_highest'])
        return (data, config)

    def optimize(self, config, sim_label, tally_label, num_cores, start_point, max_iterations, xtol_rel, ftol_rel, overwrite):
        data, config = self._preflight_optimize(config, sim_label, tally_label, num_cores, start_point, max_iterations, xtol_rel, ftol_rel, overwrite)
        optimize_analysis_start_time = timer() # used for saving elapsed time in meta
        print('[+] Constructing GeneratingFunction...')
        gf = lib.math.config_to_gf(config)
        gfEvaluatorObj = togimble.gfEvaluator(gf, config['max_k'], MUTYPES, config['gimble']['precision'], exclude=[(2,3),])
        print('[+] GeneratingFunctions for model %r have been generated.' % config['gimble']['model'])
        for data_idx, dataset in data:
            optimize_instance_start_time = timer()
            optimize_result = lib.math.optimize(gfEvaluatorObj, data_idx, dataset, config)
            optimize_time = format_time(timer() - optimize_instance_start_time)
            self.save_optimize_instance(data_idx, config, optimize_result, optimize_time, overwrite)
        config['optimize_time'] = format_time(timer() - optimize_analysis_start_time)
        self._set_meta(config['optimize_key'], meta=config_to_meta(config, 'optimize'))

        print("[+] Optimization results saved under label %r" % config['optimize_key'])

    def save_optimize_instance(self, data_idx, config, optimize_result, optimize_time, overwrite):
        data_idx = int(data_idx)
        optimize_key = config['optimize_keys'][data_idx]
        optimize_meta = {}
        optimize_meta['nlopt_log_iteration_header'] = optimize_result['nlopt_log_iteration_header']
        optimize_meta['optimize_results'] = []
        for dataset_idx in range(optimize_result['dataset_count']):
            result = optimize_result['nlopt_values_by_dataset_idx'][dataset_idx]
            result['nlopt_exit_code'] = optimize_result['nlopt_status_by_dataset_idx'][dataset_idx]
            result['likelihood'] = optimize_result['nlopt_optimum_by_dataset_idx'][dataset_idx]
            result['nlopt_time'] = optimize_time
            optimize_meta['optimize_results'].append(result)
        self._set_meta_and_data(optimize_key, optimize_meta, optimize_result['nlopt_log_iteration_array'])

    def _preflight_tally(self, data_source, data_label, max_k, sample_sets, sequence_ids, genome_file, overwrite):
        config = {
            'data_key': 'windows' if data_source == 'windowsum' else data_source,
            'data_source': data_source,
            'data_type': 'windows' if data_source == 'windowsum' else data_source,
            'data_label': data_label,
            'tally_key': self._get_key(task='tally', data_label=data_label),
            'max_k': max_k,
            'sample_sets': sample_sets,
            'sequences': sequence_ids,
            'genome_file': genome_file,
            'blocks': 0,
            'windows': 0,
            'marginalty': '0.0%',
            'block_length': 0
            }
        # check data is there
        if not self._has_key(config['data_key']):
            sys.exit("[X] gimbleStore has no %r." % config['data_key'])
        config['block_length'] = self._get_meta('blocks')['length'] # needs fixing if multiple block-datasets
        # check tally key
        if not overwrite and self._has_key(config['tally_key']):
            sys.exit("[X] Tally with data_label %r already exist. Change the label or use '--overwrite'" % (data_label))
        # sort out sequences (this could be prettier)
        config['sequences'] = list(self._get_meta('blocks')['count_by_sequence'].keys()) if config['sequences'] is None else [_ for _ in config['sequences'].split(",") if _]
        if config['genome_file']:
            try:
                df = parse_csv(csv_f=config['genome_file'], sep="\t", usecols=[0], dtype={'sequence_id': 'category'}, header=None)
                config['sequences'] = [item for sublist in df.values.tolist() for item in sublist] 
            except:
                sys.exit("[X] Could not parse %r. Should be list of sequence names." % config['genome_file'])
        return config
        
    def tally(self, data_source, data_label, max_k, sample_sets, sequence_ids, genome_file, overwrite, verbose=True):
        config = self._preflight_tally(data_source, data_label, max_k, sample_sets, sequence_ids, genome_file, overwrite)
        variation = self._get_variation(data_type=config['data_type'], sample_sets=config['sample_sets'], sequences=config['sequences'], progress=verbose) 
        config['windows'] = variation.shape[0] if variation.ndim == 3 else 0
        config['blocks'] = (variation.shape[0] * variation.shape[1]) if variation.ndim == 3 else variation.shape[0]
        config['marginality'] = format_percentage(calculate_marginality_of_variation(variation, max_k=config['max_k']))
        if verbose:
            if config['data_type'] == 'blocks':
                print("[+] Found %s blocks on %s sequence(s)." % (
                    format_count(config['blocks']), format_count(len(config['sequences'])))) 
            else:
                print("[+] Found %s blocks in %s windows (%s blocks per window) on %s sequence(s)." % (
                    format_count(config['blocks']), format_count(config['windows']), format_count(variation.shape[1]), format_count(len(config['sequences'])))) 
            print('[+] Percentage of blocks treated as marginals (w/ kmax = %s) = %s' % (config['max_k'], config['marginality']))
            print("[+] Tally'ing variation data ... ")
        if config['data_type'] == 'windowsum':
            variation = variation.reshape(-1, variation.shape[-1])
        variation_tally = tally_variation(variation, form='bsfs', max_k=config['max_k'])
        self.save_tally(config, variation_tally, verbose)
        return config['tally_key']

    def save_tally(self, config, tally, verbose=True):
        tally_key = config['tally_key'] 
        tally_meta = config_to_meta(config, 'tally')
        self._set_meta_and_data(tally_key, tally_meta, tally)
        if verbose:
            print("[+] Tally saved under label %r (use this for optimize/gridsearch)." % config['data_label'])

    def _get_key(self, task=None, data_label=None, grid_label=None, analysis_label=None, parameter_label=None, mod_label=None, seq_label=None):
        if task == 'tally':
            return "tally/%s" % data_label
        if task == 'measure':
            if data_label is not None and seq_label is not None and mod_label is not None:
                return "%s/%s/%s/%s" % (task, data_label, seq_label, mod_label)
            return "%s/" % (task) 
        if task == 'blocks' or task == 'windows':
            if seq_label is not None:
                return "%s/%s" % (task, seq_label)
            return "%s/" % (task)
        if task == 'simulate':
            if parameter_label is None:
                return "%s/%s" % (task, analysis_label) 
            return "%s/%s/%s" % (task, analysis_label, parameter_label)
        if task == 'makegrid':
            return "makegrid/%s" % (analysis_label)
        if task == 'optimize':
            if parameter_label is None:
                return "%s/%s/%s" % (task, data_label, analysis_label)
            return "%s/%s/%s/%s" % (task, data_label, analysis_label, parameter_label)
        if task == 'gridsearch':
            if data_label is not None and grid_label is not None:
                return "%s/%s/%s" % (task, data_label, grid_label)
            if parameter_label is not None:
                return "%s/%s/%s/%s" % (task, data_label, analysis_label, parameter_label)
            if mod_label is not None:
                return "%s/%s/%s/%s" % (task, data_label, analysis_label, mod_label)
            return "%s/%s/%s" % (task, data_label, analysis_label)
        return None

    def _has_key(self, key):
        return (key in self.data) if key else False

    def _set_data(self, key, array):
        self.data.create_dataset(key, data=array, overwrite=True)

    def _set_meta_and_data(self, key, meta, array):
        self.data.create_dataset(key, data=array, overwrite=True)
        self.data[key].attrs.put(meta)
    
    def _del_data_and_meta(self, key):
        if self._has_key(key):
            del self.data[key]

    def _get_data(self, key, dtype=None):
        if self._has_key(key):
            data = self.data[key]
            if dtype is None:
                return self.data[key]
            else:    
                if np.can_cast(data.dtype, dtype):
                    return self.data[key]
                else:
                    raise ValueError("Can't cast %r to %r" % (data.dtype, dtype))
        return None

    def _get_meta(self, key):
        if self._has_key(key):
            return self.data[key].attrs
        return None

    def _set_meta(self, key, meta={}):
        self.data.require_group(key)
        self.data[key].attrs.put(meta)

    def _set_stopping_criteria(self, data, parameterObj, label):
        set_by_user = True
        #if parameterObj.ftol_rel<0:
        #    set_by_user = False
        #    lnCL_sd, lnCL_all_data = self._get_lnCL_SD(data, parameterObj, label)
        #    parameterObj.ftol_rel = abs(1.96*lnCL_sd/lnCL_all_data)
        print("Stopping criteria for this optimization run:")
        print(f"Max number of evaluations: {parameterObj.max_eval}")
        if set_by_user and parameterObj.ftol_rel>0:
            print(f"Relative tolerance on lnCL: {parameterObj.ftol_rel}")
        else:
            pass
            #print(f"Relative tolerance on lnCL set by data resampling: {parameterObj.ftol_rel}")
        if parameterObj.xtol_rel>0:

            print(f"Relative tolerance on norm of parameter vector: {parameterObj.xtol_rel}")
        
    def _get_lnCL_SD(self, all_data, parameterObj, label):
        if parameterObj.data_type == 'windows':
            data = np.sum(all_data, axis=0)
            #if windows transform to sum_windowwise_bsfs
        elif parameterObj.data_type == 'simulate':    
            #if simulate data contains all replicates
            #if datatype=simulate all_data is an iterator! adapt this!
            all_data = self._get_sims_bsfs(label, single=True)
            data = np.sum(all_data, axis=0)
        else:
            data = all_data
        total_num_blocks = np.sum(data)
        ETPs = data/total_num_blocks
        
        if parameterObj.data_type=='simulate':
            #no need to resample, we have enough replicates
            #we treat the total dataset as approximating the truth
            #and turn those into ETPs
            lnCLs_resample = [lib.math.calculate_composite_likelihood(ETPs, replicate) for replicate in all_data]
            replicates = len(lnCLs_resample)
            resample_size = replicates
        else:
            replicates=100
            resample_size=min(1000, total_num_blocks)
            kmax_by_mutype=list(parameterObj.config['k_max'].values())
            resample_bsfs=self._resample_bsfs(data, resample_size, replicates, kmax_by_mutype)
            #get lnCL for each resampled dataset
            lnCLs_resample = [lib.math.calculate_composite_likelihood(ETPs, resample) for resample in resample_bsfs]
            #determine sd of outcome
        sd_mini = np.std(lnCLs_resample, ddof=1)
        sd = sd_mini * np.sqrt((resample_size-1)/(total_num_blocks-1))
        lnCL_all_data = lib.math.calculate_composite_likelihood(ETPs, data)
        return (sd, lnCL_all_data)

    def _resample_bsfs(self, data, resample_size, replicates, kmax_by_mutype):
        bsfs_2d = lib.gimble.bsfs_to_2d(data)
        bsfs_configs = bsfs_2d[:,1:]
        bsfs_probs = bsfs_2d[:,0]
        bsfs_probs=bsfs_probs/np.sum(bsfs_probs)
        resample=np.random.choice(np.arange(bsfs_configs.shape[0]), (replicates, resample_size), replace=True, p=bsfs_probs)
        resample_bsfs=bsfs_configs[resample]
        resample_bsfs=resample_bsfs.reshape(-1,len(kmax_by_mutype))
        index=np.repeat(np.arange(replicates), resample_size).reshape(replicates*resample_size,1)
        resample_bsfs_index=np.concatenate([index, resample_bsfs], axis=1)
        mutuples, counts = np.unique(resample_bsfs_index, return_counts=True, axis=0)
        shape_out = [replicates,]+[k+2 for k in kmax_by_mutype]
        out = np.zeros(tuple(shape_out), np.int64)
        out[tuple(mutuples.T)]=counts
        return out    

    def _preflight_makegrid(self, config, overwrite):
        key = self._get_key(task='makegrid', analysis_label=config['gimble']['label'])
        if not overwrite and self._has_key(key):
            sys.exit("[X] Grid with label %r already exist. Change the label in the config file or use '--force'" % config['gimble']['label'])
        config['key'] = key
        return config

    def makegrid(self, config, num_cores, overwrite):
        # print(self.data.tree())
        config = self._preflight_makegrid(config, overwrite)
        print("[+] Grid of %s grid points will be prepared..." % config['parameters_grid_points'])
        gf = lib.math.config_to_gf(config)
        print('[+] Generating equations...')
        gfEvaluatorObj = togimble.gfEvaluator(gf, config['max_k'], MUTYPES, config['gimble']['precision'], exclude=[(2,3),])
        print('[+] Equations for model %r have been generated.' % config['gimble']['model'])
        grid = lib.math.new_calculate_all_ETPs(
            gfEvaluatorObj, 
            config['parameters_expanded'], 
            config['populations']['reference_pop_id'], 
            config['mu']['block_length'], 
            config['mu']['mu'], 
            processes=num_cores, 
            verbose=False
            )
        self.save_grid(config, grid)

    def save_grid(self, config, grid):
        grid_meta = config_to_meta(config, 'makegrid')
        self._set_meta_and_data(config['key'], grid_meta, grid)

    def _validate_seq_names(self, sequences=None):
        """Returns valid seq_names in sequences or exits."""
        meta = self._get_meta('seqs')
        if sequences is None:
            return meta['seq_names']
        if set(sequences).issubset(set(meta['seq_names'])):
            return sequences
        sys.exit("[X] Sequence(s) %r not a subset of sequence(s) %r in ZARR store" % (", ".join(sequences), ", ".join(meta['seq_names'])))

    def _get_sample_set_idxs(self, query='X'):
        """Returns list of sample_set_idxs.

        Parameters 
        ----------
        query : string or None
                'X' - inter-population sample_sets
                'A' - intra-population sample_sets of population A
                'B' - intra-population sample_sets of population B 
                None - all sample_sets
        
        Returns
        -------
        out : list of strings
            sample_set_idxs that can be used to access data in gimble store
        """
        meta = self._get_meta('seqs')
        if query is None:
            return [str(idx) for idx in range(len(meta['sample_sets']))]
        elif query == 'X':
            return [str(idx) for (idx, is_cartesian) in enumerate(meta['sample_sets_inter']) if is_cartesian]
        elif query == 'A':
            return [str(idx) for (idx, is_intra_A) in enumerate(meta['sample_sets_intra_A']) if is_intra_A]
        elif query == 'B':
            return [str(idx) for (idx, is_intra_B) in enumerate(meta['sample_sets_intra_B']) if is_intra_B]
        else:
            raise ValueError("'query' must be 'X', 'A', 'B', or None")

    def _get_variation(self, data_type=None, sample_sets='X', sequences=None, population_by_letter=None, progress=False):
        """Returns variation array of 2 (blocks) or 3 (windows) dimensions.
        
        Parameters 
        ----------
        data_type : 'blocks' or 'windows'
        sample_sets : only needed for data_type 'blocks'. String or None
                None - all sample_sets 
                'X' - inter-population sample_sets (default)
                'A' - intra-population sample_sets of population A
                'B' - intra-population sample_sets of population B 
            If supplied, array is based only on variation in those sample_sets
        sequences : list of strings or None
            If supplied, array is based only on those sequences
        population_by_letter : dict (string -> string) or None
            Mapping of population IDs to population letter in model (from INI file)

        Returns
        -------
        out : ndarray, int, ndim (mutypes)
        """
        meta = self._get_meta('seqs')
        sequences = self._validate_seq_names(sequences)
        if population_by_letter:
            assert (set(population_by_letter.values()) == set(meta['population_by_letter'].values())), 'population_by_letter %r does not equal populations in ZARR store (%r)' % (population_by_letter, meta['population_by_letter'])
        keys = []
        if data_type == 'blocks':
            sample_set_idxs = self._get_sample_set_idxs(query=sample_sets)
            keys = ['blocks/%s/%s/variation' % (seq_name, sample_set_idx) 
                for seq_name, sample_set_idx in list(itertools.product(sequences, sample_set_idxs))]
            #print('variation keys', keys)
        elif data_type == 'windows':
            keys = ['windows/%s/variation' % (seq_name) for seq_name in sequences]
        else:
            raise ValueError("Invalid datatype: %s" % data_type)
        variations = []
        for key in tqdm(keys, total=len(keys), desc="[%] Preparing data...", ncols=100, unit_scale=True, disable=(not progress)):
            #variations.append(np.array(self.data[key], dtype=np.int64))
            variations.append(self.data[key])
        variation = np.concatenate(variations, axis=0)
        polarise_true = (
            (population_by_letter['A'] == meta['population_by_letter']['B']) and 
            (population_by_letter['B'] == meta['population_by_letter']['A'])) if population_by_letter else False
        if polarise_true:
            variation[..., [0, 1]] = variation[..., [1, 0]]
        return variation

    def _get_sims_bsfs(self, key):
        if self._has_key(key):
            return self.data[key].arrays()
        return None

    def _init_store(self, create, overwrite):
        if create:
            if os.path.isdir(self.path):
                print("[-] Gimble datastore %r already exists." % self.path)
                if not overwrite:
                    print("[X] Please specify '-f' to overwrite.")
                    sys.exit(1)
                print("[+] Deleting existing Gimble datastore %r" % self.path)
                shutil.rmtree(self.path)
            print("[+] Creating Gimble datastore in %r" % self.path)
            return zarr.open(str(self.path), mode='w')
        #print("[+] Loading GStore from %r" % self.path)
        return zarr.open(str(self.path), mode='r+')

    def _return_group_last_integer(self, name):
        try:
            all_groups = [int([namestring for namestring in groupnames.split('_')][-1]) for groupnames in list(self.data[name]) if groupnames.startswith('run')]
        except KeyError:
            return 0
        if len(all_groups):
            return max(all_groups)+1
        else:
            return 0

    def _plot_intervals(self):
        # [needs fixing]
        pass
        # QC plots
        #intervals_df['distance'] = np.where((intervals_df['sequence'] == intervals_df['sequence'].shift(-1)), (intervals_df['start'].shift(-1) - intervals_df['end']) + 1, np.nan)
        #distance_counter = intervals_df['distance'].dropna(how="any", inplace=False).value_counts()
        #length_counter = intervals_df['length'].value_counts()
        #distance_f = "%s.intervals.distance.png" % parameterObj.outprefix
        #plot_loglog(distance_counter, 'Distance to downstream BED interval', distance_f)
        #length_f = "%s.intervals.length.png" % parameterObj.outprefix
        #plot_loglog(length_counter, 'Length of BED interval', length_f)
        #count_sequences = intervals_df['sequence'].nunique()
        #count_intervals = len(intervals_df.index)
        #count_samples = len(query_samples)

    def _plot_blocks(self, parameterObj):
        # [needs fixing]
        mutypes_inter_key = 'seqs/bsfs/inter/mutypes' 
        counts_inter_key = 'seqs/bsfs/inter/counts' 
        mutypes_inter = self.data[mutypes_inter_key]
        counts_inter = self.data[counts_inter_key]
        self.plot_bsfs_pcp('%s.bsfs_pcp.png' % self.prefix, mutypes_inter, counts_inter)

    def _get_window_bed(self):
        meta_seqs = self._get_meta('seqs')
        meta_windows = self._get_meta('windows')
        MAX_SEQNAME_LENGTH = max([len(seq_name)+1 for seq_name in meta_seqs['seq_names']])
        window_count = meta_windows['count']
        index = np.arange(window_count)
        sequences = np.zeros(window_count, dtype='<U%s' % MAX_SEQNAME_LENGTH)
        starts = np.zeros(window_count, dtype=np.int64)
        ends = np.zeros(window_count, dtype=np.int64)
        pos_mean = np.zeros(window_count, dtype=np.float64)
        pos_median = np.zeros(window_count, dtype=np.float64)
        balance = np.zeros(window_count, dtype=np.float64)
        mse_sample_set_cov = np.zeros(window_count, dtype=np.float64)
        offset = 0
        for seq_name in tqdm(meta_seqs['seq_names'], total=len(meta_seqs['seq_names']), desc="[%] Preparing data...", ncols=100, unit_scale=True): 
            start_key = 'windows/%s/starts' % (seq_name)
            end_key = 'windows/%s/ends' % (seq_name)
            pos_mean_key = 'windows/%s/pos_mean' % (seq_name)
            pos_median_key = 'windows/%s/pos_median' % (seq_name)
            balance_key = 'windows/%s/balance' % (seq_name)
            mse_sample_set_cov_key = 'windows/%s/mse_sample_set_cov' % (seq_name)
            if start_key in self.data:
                start_array = np.array(self.data[start_key])
                _window_count = start_array.shape[0]
                starts[offset:offset+_window_count] = start_array
                ends[offset:offset+_window_count] = np.array(self.data[end_key])
                pos_mean[offset:offset+_window_count] = np.array(self.data[pos_mean_key])
                pos_median[offset:offset+_window_count] = np.array(self.data[pos_median_key])
                balance[offset:offset+_window_count] = np.array(self.data[balance_key])
                mse_sample_set_cov[offset:offset+_window_count] = np.array(self.data[mse_sample_set_cov_key])
                sequences[offset:offset+_window_count] = np.full_like(_window_count, seq_name, dtype='<U%s' % MAX_SEQNAME_LENGTH)
                offset += _window_count
        return (sequences, starts, ends, index, pos_mean, pos_median, balance, mse_sample_set_cov)

    def _preflight_windows(self, overwrite=False):
        '''
        [TODO]: allow for multiple window datasets
        '''
        if not self.has_stage('blocks'):
            sys.exit("[X] GStore %r has no blocks. Please run 'gimble blocks'." % self.path)
        if self.has_stage('windows'):
            if not overwrite:
                sys.exit("[X] GStore %r already contains windows.\n[X] These windows => %r\n[X] Please specify '--force' to overwrite." % (self.path, self.get_stage('windows')))
            print('[-] GStore %r already contains windows. But these will be overwritten...' % (self.path))
            self._del_data_and_meta('windows')
    
    def _preflight_blocks(self, overwrite=False):
        '''
        [TODO]: allow for multiple block datasets
        '''
        if not self.has_stage('measure'):
            sys.exit("[X] GStore %r has no data. Please run 'gimble measure'." % self.path)
        if self.has_stage('blocks'):
            if not overwrite:
                sys.exit("[X] GStore %r already contains blocks.\n[X] These blocks => %r\n[X] Please specify '--force' to overwrite." % (self.path, self.get_stage('blocks')))
            print('[-] GStore %r already contains blocks. But these will be overwritten...' % (self.path))
            # wipe bsfs, windows, AND meta, since new blocks...
            self._del_data_and_meta('blocks')
            self._del_data_and_meta('windows')
            self._del_data_and_meta('bsfs')

    def _get_interval_coordinates(self, seq_name=None, sample_set=None):
        if seq_name is None:
            raise ValueError('_get_interval_coordinates: needs seq_name')
        matrix_key = 'seqs/%s/intervals/matrix' % seq_name
        if matrix_key in self.data:
            start_key = 'seqs/%s/intervals/starts' % seq_name
            end_key = 'seqs/%s/intervals/ends' % seq_name
            if sample_set is None:
                return (np.array(self.data[start_key]), np.array(self.data[end_key])) 
            meta_seqs = self._get_meta('seqs')
            try:
                sample_set_key = np.array([meta_seqs['intervals_idx_by_sample'][sample] for sample in sample_set])
            except KeyError:
                sys.exit('_get_interval_coordinates: sample_set %s not found in store. Existing samples: %s' % (sample_set, list(meta_seqs['intervals_idx_by_sample'].keys())))
            mask = np.all(np.array(self.data[matrix_key])[:,sample_set_key], axis=1)
            return (np.array(self.data[start_key])[mask], np.array(self.data[end_key])[mask])
        return (None, None)

    def _get_variants(self, seq_name):
        pos_key = "seqs/%s/variants/pos" % (seq_name)
        gt_key = "seqs/%s/variants/matrix" % (seq_name)
        pos = np.array(self.data[pos_key], dtype=np.int64) if pos_key in self.data else None
        gt_matrix = allel.GenotypeArray(self.data[gt_key].view(read_only=True)) if gt_key in self.data else None
        if pos is not None and gt_matrix is not None:
            assert pos.shape[0] == gt_matrix.shape[0] # check whether they are the same length ...
        return (pos, gt_matrix)

    def _make_blocks(self, block_length, block_span, block_max_missing, block_max_multiallelic):
        meta_seqs = self._get_meta('seqs')
        blocks_raw_by_sample_set_idx = collections.Counter()   # all possible blocks
        blocks_by_sample_set_idx = collections.Counter()       # all valid blocks => only these get saved to store
        blocks_by_sequence = collections.Counter()             # all valid blocks 
        with tqdm(total=(len(meta_seqs['seq_names']) * len(meta_seqs['sample_sets'])), desc="[%] Building blocks ", ncols=100, unit_scale=True) as pbar:        
            for seq_name in meta_seqs['seq_names']:
                pos, gt_matrix = self._get_variants(seq_name) # arrays or None
                for sample_set_idx, sample_set in enumerate(meta_seqs['sample_sets']):
                    # get BED starts/ends from store
                    intervals = self._get_interval_coordinates(seq_name=seq_name, sample_set=sample_set)
                    # turn BED starts/ends into sites-array
                    sites = intervals_to_sites(intervals)
                    # turn sites-array into 2D np.array with block sites (or None)
                    blocks = sites_to_blocks(sites, block_length, block_span) 
                    if blocks is not None:
                        # subset gts of sample_set from gt_matrix (or None)
                        gts = subset_gt_matrix(meta_seqs, sample_set, 
                            np.isin(pos, blocks, assume_unique=True), gt_matrix)
                        # get block arrays
                        starts, ends, multiallelic, missing, monomorphic, variation = blocks_to_arrays(blocks, gts, pos)
                        # save block arrays
                        blocks_raw, blocks_valid = self._set_blocks(seq_name, sample_set_idx, starts, ends, multiallelic, missing, monomorphic, variation, block_max_missing, block_max_multiallelic)
                        # record counts
                        blocks_raw_by_sample_set_idx[sample_set_idx] += blocks_raw
                        blocks_by_sample_set_idx[sample_set_idx] += blocks_valid
                        blocks_by_sequence[seq_name] += blocks_valid
                    pbar.update(1)
        # save blocks meta
        self._set_blocks_meta(
            block_length, 
            block_span, 
            block_max_missing, 
            block_max_multiallelic, 
            blocks_raw_by_sample_set_idx, 
            blocks_by_sample_set_idx,
            blocks_by_sequence)

    def _set_blocks(self, seq_name, sample_set_idx, starts, ends, multiallelic, 
            missing, monomorphic, variation, block_max_missing, block_max_multiallelic):
        valid = (np.less_equal(missing, block_max_missing) & np.less_equal(multiallelic, block_max_multiallelic)).flatten()
        blocks_starts_key = 'blocks/%s/%s/starts' % (seq_name, sample_set_idx)
        self.data.create_dataset(blocks_starts_key, data=starts[valid], overwrite=True)
        blocks_ends_key = 'blocks/%s/%s/ends' % (seq_name, sample_set_idx)
        self.data.create_dataset(blocks_ends_key, data=ends[valid], overwrite=True)
        blocks_variation_key = 'blocks/%s/%s/variation' % (seq_name, sample_set_idx)
        self.data.create_dataset(blocks_variation_key, data=variation[valid], overwrite=True)
        blocks_missing_key = 'blocks/%s/%s/missing' % (seq_name, sample_set_idx)
        self.data.create_dataset(blocks_missing_key, data=missing[valid], overwrite=True)
        blocks_multiallelic_key = 'blocks/%s/%s/multiallelic' % (seq_name, sample_set_idx)
        self.data.create_dataset(blocks_multiallelic_key, data=multiallelic[valid], overwrite=True)
        return (valid.shape[0], valid[valid==True].shape[0])

    def _set_blocks_meta(self, block_length, block_span, block_max_missing, block_max_multiallelic, 
            blocks_raw_by_sample_set_idx, blocks_by_sample_set_idx, blocks_by_sequence):
        meta_blocks = self._get_meta('blocks')
        meta_blocks['length'] = block_length
        meta_blocks['span'] = block_span
        meta_blocks['max_missing'] = block_max_missing
        meta_blocks['max_multiallelic'] = block_max_multiallelic
        meta_blocks['count_by_sample_set_idx'] = dict(blocks_by_sample_set_idx) # keys are strings
        meta_blocks['count_by_sequence'] = dict(blocks_by_sequence) # keys are strings
        meta_blocks['count_raw_by_sample_set_idx'] = dict(blocks_raw_by_sample_set_idx) # keys are strings
        meta_blocks['count_total'] = sum([count for count in blocks_by_sample_set_idx.values()])
        meta_blocks['count_total_raw'] = sum([count for count in blocks_raw_by_sample_set_idx.values()])

    def _get_block_coordinates(self, sample_sets=None, sequences=None):
        sequences = self._validate_seq_names(sequences)
        sample_set_idxs = self._get_sample_set_idxs(query=sample_sets)
        keys_start_end = [('blocks/%s/%s/starts' % (seq_name, sample_set_idx), 'blocks/%s/%s/ends' % (seq_name, sample_set_idx))
            for seq_name, sample_set_idx in list(itertools.product(sequences, sample_set_idxs))]
        block_starts, block_ends = [], []
        for start_key, end_key in keys_start_end:
            block_starts.append(np.array(self.data[start_key], dtype=np.int64))
            block_ends.append(np.array(self.data[end_key], dtype=np.int64))
        block_start = np.concatenate(block_starts, axis=0)
        block_end = np.concatenate(block_ends, axis=0)
        return (block_start, block_end)

    def _get_block_sample_set_idxs(self, sample_sets=None, sequences=None):
        sequences = self._validate_seq_names(sequences)
        sample_set_idxs = self._get_sample_set_idxs(query=sample_sets)
        key_by_sample_set_idx = {sample_set_idx: 'blocks/%s/%s/starts' % (seq_name, sample_set_idx) 
                for seq_name, sample_set_idx in list(itertools.product(sequences, sample_set_idxs))}
        block_sample_set_idxs = []
        for sample_set_idx, key in key_by_sample_set_idx.items():
            block_sample_set_idxs.append(np.full(self.data[key].shape[0], int(sample_set_idx)))
        block_sample_set_idxs = np.concatenate(block_sample_set_idxs, axis=0)
        return block_sample_set_idxs

    def _make_windows(self, window_size, window_step, sample_sets='X'):
        meta_blocks = self._get_meta('blocks')
        sample_set_idxs = np.array(self._get_sample_set_idxs(query=sample_sets), dtype=np.int64)
        window_count = 0
        window_size_effective = window_size #* sample_set_idxs.shape[0]
        window_step_effective = window_step #* sample_set_idxs.shape[0]

        blockable_seqs, unblockable_seqs = [], []
        for seq_name, block_count in meta_blocks['count_by_sequence'].items():
            if block_count >= window_size_effective:
                blockable_seqs.append(seq_name)
            else:
                unblockable_seqs.append(seq_name)
        if not blockable_seqs:
            sys.exit("[X] Not enough blocks to make windows of this size (%s)." % (window_size_effective))
        print("[+] Making windows along %s sequences (%s sequences excluded)" % (len(blockable_seqs), len(unblockable_seqs)))
        for seq_name in tqdm(blockable_seqs, total=len(blockable_seqs), desc="[%] Making windows", ncols=100):
            block_variation = self._get_variation(data_type='blocks', sample_sets=sample_sets, sequences=[seq_name])
            #print('block_variation.shape', block_variation.shape)
            block_starts, block_ends = self._get_block_coordinates(sample_sets=sample_sets, sequences=[seq_name])
            #print('block_starts', block_starts)
            #print('block_ends', block_ends)
            block_sample_set_idxs = self._get_block_sample_set_idxs(sample_sets=sample_sets, sequences=[seq_name])
            windows = blocks_to_windows(sample_set_idxs, block_variation, block_starts, block_ends, block_sample_set_idxs, window_size_effective, window_step_effective)
            window_variation, window_starts, window_ends, window_pos_mean, window_pos_median, balance, mse_sample_set_cov = windows
            window_count += self._set_windows(seq_name, window_variation, window_starts, window_ends, window_pos_mean, window_pos_median, balance, mse_sample_set_cov)
        self._set_windows_meta(sample_set_idxs, window_size_effective, window_step_effective, window_count)

    def _set_windows(self, seq_name, window_variation, window_starts, window_ends, window_pos_mean, window_pos_median, balance, mse_sample_set_cov):
        self.data.create_dataset("windows/%s/variation" % seq_name, data=window_variation, overwrite=True)
        self.data.create_dataset("windows/%s/starts" % seq_name, data=window_starts, overwrite=True)
        self.data.create_dataset("windows/%s/ends" % seq_name, data=window_ends, overwrite=True)
        self.data.create_dataset("windows/%s/pos_mean" % seq_name, data=window_pos_mean, overwrite=True)
        self.data.create_dataset("windows/%s/pos_median" % seq_name, data=window_pos_median, overwrite=True)
        self.data.create_dataset("windows/%s/balance" % seq_name, data=balance, overwrite=True)
        self.data.create_dataset("windows/%s/mse_sample_set_cov" % seq_name, data=mse_sample_set_cov, overwrite=True)
        return window_variation.shape[0]

    def _set_windows_meta(self, sample_set_idxs, window_size, window_step, window_count):
        meta_windows = self._get_meta('windows')
        meta_windows['size_user'] = int(window_size / len(sample_set_idxs))
        meta_windows['step_user'] = int(window_step / len(sample_set_idxs))
        meta_windows['size'] = window_size 
        meta_windows['step'] = window_step
        meta_windows['count'] = window_count

####################### REPORTS ######################

    def _get_parse_report(self, width):
        meta_seqs = self._get_meta('seqs')
        reportObj = ReportObj(width=width)
        reportObj.add_line(prefix="[+]", left='[', center='Parsed data', right=']', fill='=')
        reportObj.add_line(prefix="[+]", left='seqs')
        right = "%s in %s sequence(s) (n50 = %s)" % (
            format_bases(sum(meta_seqs['seq_lengths'])), 
            format_count(len(meta_seqs['seq_lengths'])), 
            format_bases(meta_seqs['seq_n50']))
        reportObj.add_line(prefix="[+]", branch="T", fill=".", left='genome', right=right)
        right = "%s samples in %s populations" % (format_count(len(meta_seqs['samples'])), format_count(len(meta_seqs['population_ids'])))
        reportObj.add_line(prefix="[+]", branch="T", left='populations', fill=".", right=right)
        sample_counts_by_population = collections.Counter(meta_seqs['populations'])
        for idx, (letter, population_id) in enumerate(meta_seqs['population_by_letter'].items()):
            left = "%s = %s" % (letter, population_id)
            right = "%s" % format_count(sample_counts_by_population[population_id])
            branch = "P%s" % ("F" if idx == len(meta_seqs['population_by_letter']) -1 else "T")
            reportObj.add_line(prefix="[+]", branch=branch, left=left, right=right)
        reportObj.add_line(prefix="[+]", branch="T", left="sample sets", fill=".", right=format_count(len(meta_seqs['sample_sets'])))
        reportObj.add_line(prefix="[+]", branch="PT", left="INTER-population sample-sets (X)", right=format_count(len(self._get_sample_set_idxs(query='X'))))
        reportObj.add_line(prefix="[+]", branch="PT", left="INTRA-population sample-sets (A)", right=format_count(len(self._get_sample_set_idxs(query='A'))))
        reportObj.add_line(prefix="[+]", branch="PF", left="INTRA-population sample-sets (B)", right=format_count(len(self._get_sample_set_idxs(query='B'))))
        reportObj.add_line(prefix="[+]", branch="T", left="variants", fill=".", right=("%s (%s per 1 kb)" % (
               format_count(meta_seqs['variants_counts']),
               format_proportion(1000 * meta_seqs['variants_counts'] / sum(meta_seqs['seq_lengths'])))))
        reportObj.add_line(prefix="[+]", branch="PP", right="".join([c.rjust(8) for c in ["HOMREF", "HOMALT", "HET", "MISS"]]))
        for idx, sample in enumerate(meta_seqs['samples']):
            variant_idx = meta_seqs['variants_idx_by_sample'][sample]
            branch = "PF" if idx == len(meta_seqs['variants_idx_by_sample']) - 1 else "PT"
            left = sample
            right = "%s %s %s %s" % (
            (format_percentage(meta_seqs['variants_counts_hom_ref'][variant_idx] / meta_seqs['variants_counts']) if meta_seqs['variants_counts'] else format_percentage(0)).rjust(7),
            (format_percentage(meta_seqs['variants_counts_hom_alt'][variant_idx] / meta_seqs['variants_counts']) if meta_seqs['variants_counts'] else format_percentage(0)).rjust(7),
            (format_percentage(meta_seqs['variants_counts_het'][variant_idx] / meta_seqs['variants_counts']) if meta_seqs['variants_counts'] else format_percentage(0)).rjust(7),
            (format_percentage(meta_seqs['variants_counts_missing'][variant_idx] / meta_seqs['variants_counts']) if meta_seqs['variants_counts'] else format_percentage(0)).rjust(7))
            reportObj.add_line(prefix="[+]", branch=branch, left=left, right=right)
        reportObj.add_line(prefix="[+]", branch="F", left='intervals', fill=".", 
            right = "%s intervals across %s (%s of genome)" % (
                format_count(meta_seqs['intervals_count']),
                format_bases(meta_seqs['intervals_span']),
                format_percentage(meta_seqs['intervals_span'] / sum(meta_seqs['seq_lengths']))))
        return reportObj
    
    def _get_storage_report(self, width):
        reportObj = ReportObj(width=width)
        reportObj.add_line(prefix="[+]", left='[', center='Storage', right=']', fill='=')
        total_size = recursive_get_size(self.path)
        total_right = "%s | %s" % (format_bytes(total_size), format_percentage(total_size/total_size))
        reportObj.add_line(prefix="[+]", left=pathlib.Path(self.path).name, right=total_right, fill='.')
        for idx, group in enumerate(self.data):
            size = recursive_get_size(pathlib.Path(self.path) / pathlib.Path(group))
            percentage = format_percentage(size/total_size)
            branch = branch="T" if idx < len(self.data) - 1 else 'F'
            reportObj.add_line(prefix="[+]", branch=branch, left=group, right='%s | %s' % (format_bytes(size), percentage.rjust(7)), fill='.')
        return reportObj

    def _get_blocks_report_metrics(self):
        meta_blocks = self._get_meta('blocks')
        block_length = meta_blocks['length']
        meta_seqs = self._get_meta('seqs')
        intervals_span = meta_seqs['intervals_span']
        BRMs = {}
        for sample_sets in ['X', 'A', 'B']:
            sample_sets_count = len(self._get_sample_set_idxs(sample_sets))
            tally = tally_variation(self._get_variation(data_type='blocks', sample_sets=sample_sets), form='tally')
            BRM = calculate_blocks_report_metrics(tally, sample_sets_count, block_length, intervals_span)
            for k, v in BRM.items():
                if k in ['interval_coverage', 'blocks_invariant', 'blocks_fgv']:
                    BRM[k] = format_percentage(v, precision=2)
                elif k in ['blocks_total']:
                    BRM[k] = format_count(v)
                else:
                    BRM[k] = format_proportion(v, precision=5)
            if sample_sets == 'X':
                del BRM['pi']
                del BRM['watterson_theta']
            if sample_sets == 'A' or sample_sets == 'B':
                if sample_sets == 'A':
                    del BRM['heterozygosity_B']
                    BRM['heterozygosity_A'] =  BRM['heterozygosity_intra']
                if sample_sets == 'B':
                    del BRM['heterozygosity_A'] 
                    BRM['heterozygosity_B'] =  BRM['heterozygosity_intra']
                del BRM['dxy']
                del BRM['fst']
            BRMs[sample_sets] = BRM
        return BRMs

    def _get_blocks_report(self, width):
        reportObj = ReportObj(width=width)
        reportObj.add_line(prefix="[+]", left='[', center='Blocks', right=']', fill='=')
        if self.has_stage('blocks'):
            column_just = 14
            meta_blocks = self._get_meta('blocks')
            BRMs = self._get_blocks_report_metrics()
            reportObj.add_line(prefix="[+]", left='blocks')
            reportObj.add_line(prefix="[+]", branch='T', fill=".", 
                left="'-l %s -m %s -u %s -i %s'" % (
                    meta_blocks['length'], 
                    meta_blocks['span'],
                    meta_blocks['max_missing'],
                    meta_blocks['max_multiallelic']),
                right=' %s blocks (%s discarded)' % (
                    format_count(meta_blocks['count_total']), 
                    format_percentage(1 - (meta_blocks['count_total'] / meta_blocks['count_total_raw']))))
            reportObj.add_line(prefix="[+]", branch="P", right="".join(
                [c.rjust(column_just) for c in ["X", "A", "B"]]))
            for label, key, branch in [('Mean BED interval sites in blocks (%) *', 'interval_coverage', 'T'),
                                       ('Total blocks', 'blocks_total', 'T'),
                                       ('Invariant blocks', 'blocks_invariant', 'T'),
                                       ('Four-gamete-violation blocks', 'blocks_fgv', 'T'),
                                       ('Heterozygosity (population A)', 'heterozygosity_A', 'T'),
                                       ('Heterozygosity (population B)', 'heterozygosity_B', 'T'),
                                       ('D_xy', 'dxy', 'T'),
                                       ('F_st', 'fst', 'T'),
                                       ('Pi', 'pi', 'T'),
                                       ('Watterson theta', 'watterson_theta', 'F')]:
                reportObj.add_line(
                    prefix="[+]", 
                    branch=branch, 
                    left=label, 
                    right="".join([c.rjust(column_just) for c in [BRMs['X'][key], BRMs['A'][key], BRMs['B'][key]]]))
        return reportObj

    def _get_windows_report(self, width):
        meta_windows = self._get_meta('windows')
        reportObj = ReportObj(width=width)
        reportObj.add_line(prefix="[+]", left='[', center='Windows', right=']', fill='=')
        if self.has_stage('windows'):
            reportObj.add_line(prefix="[+]", left='windows')
            reportObj.add_line(prefix="[+]", branch='F', fill=".", left="'-w %s -s %s'" % (meta_windows['size_user'], meta_windows['step_user']), right=' %s windows of inter-population (X) blocks' % (
                format_count(meta_windows['count'])))
        return reportObj

    def _get_tally_report(self, width):
        reportObj = ReportObj(width=width)
        reportObj.add_line(prefix="[+]", left='[', center='Tallies', right=']', fill='=')
        reportObj.add_line(prefix="[+]", left='Tally')
        if 'tally/' in self.data:
            for tally_label, tally_array in self.data['tally/'].arrays():
                reportObj.add_line(prefix="[+]", branch='F', fill=".", left="'tally/%s'" % (tally_label), right=' shape %s ' % (str(tally_array.shape)))
        return reportObj

    def info(self, version=None, tree=False):
        width = 100
        if tree:
            return self.data.tree()
        report = self._get_storage_report(width)
        report += self._get_parse_report(width)    
        report += self._get_blocks_report(width)
        report += self._get_windows_report(width)
        report += self._get_tally_report(width)
        #report += self._get_optimize_report(width)
        #report += self._get_grids_report(width)
        #report += self._get_lncls_report(width)
        #report += self._get_bsfs_report(width)
        #report += self._get_sims_report(width)
        write_info_report(version, report, self.prefix)
        return report

    def _get_grids_report(self, width):
        reportObj = ReportObj(width=width)
        reportObj.add_line(prefix="[+]", left='[', center='Grids', right=']', fill='=')
        for grid_id in self.data['grids/']:
            grid_path = 'grids/%s' % grid_id
            grid_dict = self.data[grid_path].attrs.asdict()
            reportObj.add_line(prefix="[+]", branch='W', fill=".", left=grid_id, right=' %s grid points' % format_count(len(grid_dict)))
            # also needs kmax, populations? sync?
            table = []
            for grid_point, grid_params in grid_dict.items():
                table.append([grid_point] + list(grid_params.values()))
            rows = tabulate.tabulate(table, numalign="right", headers=['i'] + list(grid_params.keys())).split("\n")
            for i, row in enumerate(rows):
                branch = 'F' if (i + 1) == len(rows) else 'P'
                reportObj.add_line(prefix="[+]", branch=branch, fill="", left=row, right='')
        return reportObj

    def _get_lncls_report(self, width):
        reportObj = ReportObj(width=width)
        reportObj.add_line(prefix="[+]", left='[', center='lnCLs', right=']', fill='=')
        legacy = True if 'global' in self.data['lncls/'] else False
        lncls_path = 'lncls/global' if legacy else 'lncls/'
        for grid_id in self.data[lncls_path]:
            
            shape = self.data["%s/%s" % (lncls_path, grid_id)].shape
            reportObj.add_line(prefix="[+]", branch='S', fill=".", left=grid_id, right='NA')

            # grid_dict = self.data[grid_path].attrs.asdict()
            # reportObj.add_line(prefix="[+]", branch='W', fill=".", left=grid_id, right=' %s grid points' % format_count(len(grid_dict)))
            # table = []
            # for grid_point, grid_params in grid_dict.items():
            #     table.append([grid_point] + list(grid_params.values()))
            # rows = tabulate.tabulate(table, numalign="right", headers=['i'] + list(grid_params.keys())).split("\n")
            # for i, row in enumerate(rows):
            #     branch = 'F' if (i + 1) == len(rows) else 'P'
            #     reportObj.add_line(prefix="[+]", branch=branch, fill="", left=row, right='')
        return reportObj

    def _get_sims_report(self, width, label):
        simulate_key = self._get_key(task='simulate', analysis_label=label)
        meta_sims = self._get_meta('sims')
        reportObj = ReportObj(width=width)
        reportObj.add_line(prefix="[+]", left='[', center='Sims', right=']', fill='=')
        column_just = 14
        if self.has_stage('simulate'):
            for name, array in self.get_bsfs(data_type='simulate', label=label):#returns an iterator over parametercombination names and arrays
                bsfs_X = bsfs_to_2d(np.array(array))
                fgv_blocks_X = np.sum(bsfs_X[(bsfs_X[:,4]>0) & (bsfs_X[:,5]>0)][:,0]) / np.sum(bsfs_X[:,1]) if np.any(bsfs_X) else "N/A"
                reportObj.add_line(prefix="[+]", branch='T', left=f'four-gamete-violation {name}', right="".join([format_percentage(c, precision=2).rjust(column_just) for c in [fgv_blocks_X,]]))
        return reportObj