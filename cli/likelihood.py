#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""usage: gIMble likelihood         -l FILE -A STR (-w FILE | -v FILE) 
                                            (-m FLOAT | --m_low FLOAT --m_high FLOAT)
                                            (--theta FLOAT | --theta_low FLOAT --theta_high FLOAT)
                                            (--time FLOAT | --time_low FLOAT --time_high FLOAT)
                                            [-t INT -o STR -h]

    Options:
        -h, --help                                   show this

        -w, --windows_hd5 FILE                      windows HDF5 file
        -v, --variants_hd5 FILE                     variants HDF5 file
        -k, --kmax INT                              maximum mutation count beyond which probabilities are calculated as marginals
        
        -l, --model FILE                            Model file
        
        -A, --ancestor_population STR               Ancestor population id
        -d, --derived_population_size FLOAT         Coalescence rate multiplier for derived population [default: 1.0]
        
        -m, --migration_rate FLOAT                  Migration rate per generation 
        --m_low FLOAT                               Lower bound for migration rate per generation (will be estimated)
        --m_high FLOAT                              Upper bound for migration rate per generation (will be estimated)
        
        --theta FLOAT                               Mutation rate/lineage/base
        --theta_low FLOAT                           Lower bound for mutation rate/lineage/base (will be estimated)
        --theta_high FLOAT                          Upper bound for mutation rate/lineage/base (will be estimated)
    
        --time FLOAT                                Split time in (unit?) 
        --time_low FLOAT                            Lower bound for split time ( will be estimated)
        --time_high FLOAT                           Upper bound for split time (will be estimated)

        -t, --threads <INT>                         Number of threads to use [default: 1]
        -o, --prefix <STR>                          Folder/prefix for output
        
"""

from timeit import default_timer as timer
from docopt import docopt
import lib.likelihood
import collections 

boundaries = collections.OrderedDict({
            'Time' : (0.0, 4.0),
            'theta' : (0.0, 1.0)
        })

def main():
    main_time = timer()
    args = docopt(__doc__)
    print(args)
    parameterObj = lib.likelihood.ParameterObj(args)
    print(parameterObj.base_rate)
    parameterObj.generate_data_space()
    pathObj_by_path_id = lib.likelihood.prepare_paths(parameterObj)
    symbolic_equations_by_data_tuple = lib.likelihood.generate_equations(pathObj_by_path_id, parameterObj)
    pod = lib.likelihood.calculate_pods(symbolic_equations_by_data_tuple, parameterObj)    
    lib.likelihood.estimate_parameters(symbolic_equations_by_data_tuple, boundaries, pod, parameterObj, 12345)
    #parameterObj.write_probs(prob_by_data_string, data)
    #parameterObj.write_path_equations(pathObj_by_path_id)
    print("[+] Total runtime: %s seconds" % (timer() - main_time))
        
###############################################################################

if __name__ == '__main__':
    main()