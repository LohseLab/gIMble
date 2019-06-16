#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""usage: gIMble gridsearch         -s FILE -g FILE -l FILE -A STR -w FILE [-k INT --block_size INT]
                                            [--time_MLE FLOAT --derived_MLE FLOAT --mu FLOAT] [--grid FILE]
                                            [--theta_low FLOAT --theta_high FLOAT --migration_low FLOAT --migration_high FLOAT]
                                            [-t INT -o STR -h]

    Options:
        -h, --help                                   show this

        -s, --sample_file <FILE>                    CSV file ("sample_id,population_id")
        -g, --genome_file <FILE>                    Genome file (as used in BedTools)
        -w, --windows_hd5 FILE                      windows HDF5 file
        -k, --kmax INT                              maximum count per mutation type per block [default: 2]
                                                     (beyond which counts are combined)
        
        --grid FILE                                 Grid File

        -l, --model FILE                            Model file
        
        -A, --ancestor_population STR               Ancestor population id
        --block_size INT                            Block size in bases
        
        --time_MLE FLOAT                            Split time (in 2*N_e*generations) 
        --derived_MLE FLOAT                         Coalescence rate scalar for derived population (inverse of Ne)
        --mu FLOAT                                  Scaled mutation rate per block (theta = 4 * N_e * µ)

        --theta_low FLOAT                           Lower bound scaled mutation rate (will be estimated)
        --theta_high FLOAT                          Upper bound scaled mutation rate (will be estimated)
        --migration_low FLOAT                       Lower bound migration rate per generation (M = 4 * N_e * m) 
        --migration_high FLOAT                      Upper bound migration rate per generation (M = 4 * N_e * m) 

        -t, --threads <INT>                         Number of threads to use [default: 1]
        -o, --prefix <STR>                          Folder/prefix for output
        
"""

from timeit import default_timer as timer
from docopt import docopt
import lib.gridsearch 

def main():
    main_time = timer()
    args = docopt(__doc__)
    print("[#] ### gIMble GRIDSEARCH ###")
    parameterObj = lib.gridsearch.ParameterObj(args)
    entityCollection = lib.gridsearch.task_generate_entityCollection(parameterObj)
    parameterObj.parse_window_positions(entityCollection)
    parameterObj.parse_grid_file(entityCollection)
    parameterObj.generate_mutuple_space()
    #if parameterObj.grid_raw is None:
    #    parameterObj.setup_grid()
    #print(parameterObj.boundaries)
    mutuple_count_matrix_by_window_id = parameterObj.get_mutuple_count_matrix_by_window_id(entityCollection)
    
    # needs symbolic 'C_ancestor', 'C_derived', 'Migration', 'BigL', 'hetA', 'hetB', 'hetAB', 'fixed'
    # pathObj_by_path_id = lib.gridsearch.prepare_paths(parameterObj) 
    # needs symbolic 'hetA', 'hetB', 'hetAB', 'fixed'
    # symbolic_equations_by_mutuple = lib.gridsearch.generate_equations(pathObj_by_path_id, parameterObj)
    # needs numeric grid_gimble, references grid_raw using idx
    # probability_matrix_by_grid_idx = lib.gridsearch.score_grid(symbolic_equations_by_mutuple, mutuple_count_matrix_by_window_id, parameterObj)
    # composite_likelihood_by_grid_idx_by_window_id = lib.gridsearch.infer_composite_likelihoods(probability_matrix_by_grid_idx, mutuple_count_matrix_by_window_id)
    composite_likelihood_by_parameter_tuple_by_window_id = lib.gridsearch.compute_composite_likelihoods(mutuple_count_matrix_by_window_id, parameterObj)
    lib.gridsearch.generate_output(composite_likelihood_by_parameter_tuple_by_window_id, parameterObj)
    print("[+] Total runtime: %s seconds" % (timer() - main_time))
        
###############################################################################

if __name__ == '__main__':
    main()