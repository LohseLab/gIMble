#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""usage: gIMble inference                  [-z FILE] -m FILE [-c FILE] (-b|-w)
                                            [-t INT] [-P INT] [-h|--help]
                                            
    Options:
        -h --help                                   show this
        -m, --model_file FILE                       Model file to analyse
        -z, --zarr_file FILE                        ZARR datastore
        -c, --config_file FILE                      Config file with parameters (if not present, empty config file is created)
        -P, --precision INT                         Floating point precision of probabilities [default: 30]
        -b, --blocks                                                               
        -w, --windows
        -t, --threads INT                           Threads [default: 1]
        
"""
import pathlib
import oyaml as yaml
import collections
from timeit import default_timer as timer
from docopt import docopt
import sys
from lib.gimble import RunObj
import lib.math

'''
- Add either in config or print statement: mutype => fmiC
- Config files as function of analysis type
    - ~/git/gIMble/gIMble inference -m models/graph.s_A_B.p2.n_1_1.J_A_B.M_A_B.TEST4.model.tsv --numeric
    - ~/git/gIMble/gIMble inference -m models/graph.s_A_B.p2.n_1_1.J_A_B.M_A_B.TEST4.model.tsv --grid

- In model: Fix C's as identical (thereby decreasing variables)
    ~/git/gIMble/gIMble model -s A,B -p 2 -n 1,1 -j '(A,B)' -m 'A>B' -C 'A=B' -o test

- what happens if we don't have a bed file? only intergenic/no coverage

Grid:
- which gridpoint gives highes Likelihood
    - overall
    - conditional, when Me is set to 0

# bootstrapping
after running simulations
- for each grid points, do simulation replicates
- somehow use rembination rate

1. reading in data
    - best fitting parameter combination
   - IM models:
     Div1, -3.587 × 10**6,       theta=1.119,                                               T=0.3767
    Div2, -3.587 × 10**6,        theta=1.086, λ[{x}] → 0.8375,                              T=0.4049
    Div2B, -3.585 × 10**6,       theta=1.174, λ[{y}] → 1.315,                               T=0.3295
    Div3, -3.585 × 10**6,        theta=1.161, λ[{x}] → 0.9507, λ[{y}] → 1.284,              T=0.3397
    MigChiToRos, -3.576 × 106,   theta=0.4366, λ[{y}] → 0.4147, M → 1.652       
    MigRostoChi, -3.594 × 106,   theta=1.020, λ[{x}] → 1.749, M → 3.895     
    IMChiToRosChi, -3.577 × 106, theta=0.5000, λ[{y}] → 0.4808,                             T=8.000
    IMChiToRosRos, -3.574 × 106, theta=1.055, λ[{x}] → 2.400,                               T=3.708
  
2. read in mathematic grid
    - plot Me/Ne
3. simulations  
    - checking probabilities
    - parametric bootstraps
        - recombination map or just single value
        
gimble inference (based on blocks)
    => makes config for inference
    => writes a new config file to be used in grid/simulation

gimble grid (just likelihood calculation)
    => uses 'centre' of grid in config file 

gimble simulate
    => replicates 
    => simulates windows
    (=> simulate blocks for internal benchmarking)
        - simulate replicates for each of the points in the grid 
    * parsers need checks based on shape

gimble scan (combininig windows data with likelihoods)






'''

class ParameterObj(RunObj):
    '''Sanitises command line arguments and stores parameters.'''

    def __init__(self, params, args):
        super().__init__(params)
        self.zstore = self._get_path(args['--zarr_file'])
        self.model_file = self._get_path(args['--model_file'])
        self.config_file = self._get_path(args['--config_file'])
        self.threads = self._get_int(args['--threads'])
        self._config = self._get_or_write_config()
        self.data_type = self._get_datatype([args['--blocks'], args['--windows']])

    def _get_datatype(self, args):
        if not any(args):
            return None
        elif args[0]:
            return 'blocks'
        elif args[1]:
            return 'windows'
        else:
            sys.exit("[X1] This should not have happened.")

    def _get_or_write_config(self):
        '''
        [To Do] 
            - numerical params: equal values mean equality (simplifies equation).
            - boundary params: need equality list, parameters in equality list MUST have equal values and are then collapsed.
        '''
        if self.config_file is None:
            print("[-] No config file found.")
            print("[+] Generating config file for model %r" % self.model_file)
            '''for now we use the following dict until columns are fixed in gimble model''' 
            config = {
                'version': self._VERSION,
                'random_seed' : 12345,
                'precision': 25,
                'model' : self.model_file,
                'population_ids': collections.defaultdict(dict),
                'k_max': collections.defaultdict(dict),
                'parameters': collections.defaultdict(dict), 
                'boundaries': collections.defaultdict(list),
                }
            config['parameters']['theta'] = 'FLOAT'
            for column in self._parse_model_file(target='header'):
                if column.startswith('C_'): 
                    config['parameters'][column] = 'FLOAT'
                    population_id = column.replace('C_', '') 
                    if len(population_id) == 1:
                        config['population_ids'][population_id] = 'STRING'
                if column.startswith('M_'):
                    config['parameters'][column] = 'FLOAT'
                elif column.startswith('m_'):
                    config['k_max'][column] = 'INT'
            config['parameters']['T'] = 'FLOAT'
            for parameter in config['parameters']:
                config['boundaries'][parameter] = ['MIN', 'MAX']
            config_file = pathlib.Path(self.model_file).with_suffix('.config.yaml')
            yaml.add_representer(collections.defaultdict, yaml.representer.Representer.represent_dict)
            with open(config_file, 'w') as fh:
                yaml.dump(config, fh)
            print("[+] Wrote file %r" % str(config_file))
            sys.exit("[X] Please specify parameters in config file %r" % str(config_file))
        else:
            print("[+] Reading config %r" % self.config_file)
            config_raw = yaml.safe_load(open(self.config_file, 'r'))
            config = {}
            for k, v in config_raw.items():
                if k == 'version':
                    config[k] = v
                elif k == 'population_ids':
                    config[k] = v
                elif k == 'model':
                    config[k] = v
                elif isinstance(v, str):
                    sys.exit("[X] Config file error: %r should be a number (not %r)." % (k, v))
                elif k == 'parameters':
                    config['parameters'], config['boundaries'] = {}, {}
                    for v_k, v_v in config_raw[k].items():
                        if isinstance(v_v, str): # parameter not set
                            if any([isinstance(bound, str) for bound in config_raw['boundaries'][v_k]]):
                                sys.exit("[X] Config file error: set parameter or boundaries for %r (not %r)." % (v_k, v_v))
                            else:
                                config['boundaries'][v_k] = config_raw['boundaries'][v_k]
                        else:
                            config[k][v_k] = v_v
                elif k == 'boundaries':
                    pass
                elif k == 'k_max':
                    config['k_max'] = {}
                    for v_k, v_v in config_raw[k].items():
                        if isinstance(v_v, int): # k_max not set
                            config[k][v_k] = v_v
                        else:
                            sys.exit("[X] Config file error: set value for k_max %r (not %r)." % (v_k, v_v))
                else:
                    config[k] = v
            return config

    def _parse_model_file(self, target='model_name'):
        '''# model = s_A_B.p2.n_1_1.m_AtoB.j_A_B'''
        with open(self.model_file) as fh:
            for l in fh:
                if target == 'model_name':
                    if l.startswith("# model ="):
                        return l.split(" = ")[1].rstrip("\n")
                if target == 'header':
                    if l.startswith("path_idx"):
                        return l.split()

def main(params):
    try:
        '''
        hetB, hetA, hetAB, fixed
        Ne_B, Ne_A
        B,    A
        
        grid:=
            center of grid : 
                - "global" model estimated from overal block_counts across genome in canonical workflow
                - Ne_A, Ne_B, M_e
            boundaries: 
                - all required! 
                - have to be checked for consistency! (logarithmic gridding?)
                Ne_A_min, Ne_A_max, Ne_A_steps (min=0 makes no sense, though)
                Ne_B_min, Ne_B_max, Ne_B_steps (min=0 makes no sense, though)
                M_e_min, M_e_Max, M_e_steps (def want to include min=0)
                     
                                  centre
        step                        *
        |   |   |   |   |   |   |   |   |   |   |   |   |   |   |        
        min                                                     max
        
        l = range(min, max+step, step)  # should allow for linear/log spacing                                                   
        if not centre in l:  # identity check has to allow for precision-offness
            error
        if M_e and not 0 in l: # identity check has to allow for precision-offness
            error   

        '''
        start_time = timer()
        args = docopt(__doc__)
        #log = lib.log.get_logger(params)
        parameterObj = ParameterObj(params, args)
        #data = lib.math.get_data_array(parameterObj)
        #print(data)
        equationSystem = lib.math.EquationSystemObj(parameterObj)
        equationSystem.info()
        equationSystem.initiate_model()
        PODs = equationSystem.calculate_PODs()
        print(PODs)

        
        print("[*] Total runtime: %.3fs" % (timer() - start_time))
    except KeyboardInterrupt:
        print("\n[X] Interrupted by user after %s seconds!\n" % (timer() - start_time))
        exit(-1)