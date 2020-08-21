#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""usage: gimble gridsearch                  -z FILE -m FILE [-c FILE]
                                            (-b|-w) [-t INT] [-P INT] [-h|--help]
                                            
                                            
    Options:
        -h --help                                   show this
        -m, --model_file FILE                       Model file to analyse
        -z, --zarr_file FILE                        ZARR datastore
        -c, --config_file FILE                      Config file with parameters (if not present, empty config file is created)
        -g, --grid_file FILE                        Grid file
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
import lib.gimble
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

class GridsearchParameterObj(lib.gimble.ParameterObj):
    '''Sanitises command line arguments and stores parameters.'''

    def __init__(self, params, args):
        super().__init__(params)
        self.zstore = self._get_path(args['--zarr_file'])
        self.model_file = self._get_path(args['--model_file'])
        self.config_file = self._get_path(args['--config_file'])
        self.threads = self._get_int(args['--threads'])
        self.config = self._parse_config(self.config_file)
        self.data_type = self._get_datatype([args['--blocks'], args['--windows']])
        #self.probcheck_file = self._get_path(args['--probcheck']) if args['--probcheck'] is not None else None

    def _get_datatype(self, args):
        if not any(args):
            return None
        elif args[0]:
            return 'blocks'
        elif args[1]:
            return 'windows'
        else:
            sys.exit("[X1] This should not have happened.")

    # def _parse_config(self):
    #     '''
    #     [To Do] 
    #         - numerical params: equal values mean equality (simplifies equation).
    #         - boundary params: need equality list, parameters in equality list MUST have equal values and are then collapsed.
    #     - Cerberus docs here: https://docs.python-cerberus.org/en/stable/validation-rules.html#valuesrules-rule
    #     '''
    #     print("[+] Reading config %r" % self.config_file)
    #     import oyaml
    #     import cerberus
    #     import sys
    #     validator = cerberus.Validator()
    #     schema = {
    #         'version': {'type': 'string'},
    #         'precision': {'type': 'integer'},
    #         'random_seed': {'type': 'integer'},
    #         'population_ids': {'type': 'dict', 'valuesrules': {'type': 'string'}},
    #         'k_max': {'type': 'dict', 'valuesrules': {'type': 'integer', 'min': 1}},
    #         'parameters': {'type': 'dict', 'valuesrules': {'type': 'float', 'nullable': True}},
    #         'boundaries': {'type': 'dict', 'valuesrules': {'type': 'list', 'nullable': True}}} # does not check for type IN list, has to be done afterwards
    #     config = oyaml.safe_load(open(self.config_file, 'r'))
    #     output = ["[X] YAML Config file format error(s) ..."]
    #     if not validator.validate(config, schema):
    #         for level in validator.errors:
    #             output.append("[X] %r ..." % level)
    #             for error_dict in validator.errors[level]:
    #                 for key, error_list in error_dict.items():
    #                     output.append("[X] \t %r : %s" % (key, "; ".join(error_list)))
    #     # boundary validation
    #     for key, value_list in config['boundaries'].items():
    #         if value_list:
    #             print(key, value_list)
    #             floats, string = value_list[0:4], value_list[4]
    #             if not all([isinstance(value, float) for value in floats]):
    #                 output.append("[X] \t %r : first 4 values should be floats" % (key))
    #             if not isinstance(string, str):
    #                 output.append("[X] \t %r : last value should be 'linear' or 'log'" % (key))
    #     if len(output) > 1:
    #         sys.exit("\n".join(output))
    #     return config

    #def _make_grid(self):
    #    '''
    #    - no need to be an "actual" grid
    #    '''
    #    if len(self._config["parameters"])>0:
    #        for key, value in self._config["parameters"].items():
    #            if len(value) > 1 and key!="recombination":
    #                assert len(value) >= 3, "MIN, MAX and STEPSIZE need to be specified"
    #                sim_range = np.arange(value[0], value[1]+value[2], value[2], dtype=float)
    #                if len(value)==4:
    #                    if not any(np.isin(sim_range, value[3])):
    #                        print(f"[-] Specified range for {key} does not contain specified grid center value")  
    #                self._config["parameters"][key] = sim_range
        
def main(params):
    try:
        start_time = timer()
        args = docopt(__doc__)
        parameterObj = GridsearchParameterObj(params, args)
        print(parameterObj.config)
        # grid
        sys.exit("done.")
        grid_points = lib.math.get_grid(parameterObj) # LoD
        # data
        gimbleStore = lib.gimble.Store(path=parameterObj.zstore, create=False)
        data = gimbleStore.get_bsfs_matrix(
            data='blocks', 
            population_by_letter=parameterObj._config['population_ids'], 
            cartesian_only=True, 
            kmax_by_mutype=parameterObj._config['k_max'])


        #data = lib.math.get_data_array(parameterObj)
        #print(data)
        equationSystem = lib.math.EquationSystemObj(parameterObj)
        equationSystem.info()
        equationSystem.initiate_model()
        equationSystem.calculate_ETPs()
        if parameterObj.probcheck_file is not None:
            equationSystem.check_ETPs()
        from scipy.special import xlogy
        import numpy as np
        composite_likelihood = -np.sum((xlogy(np.sign(equationSystem.ETPs), equationSystem.ETPs) * data))
        print('[+] L=-%s' % (composite_likelihood))
        print("[*] Total runtime: %.3fs" % (timer() - start_time))
    except KeyboardInterrupt:
        print("\n[X] Interrupted by user after %s seconds!\n" % (timer() - start_time))
        exit(-1)