#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""usage: gIMble simulate -s <STR> -b <INT> -r <FLOAT> -S <STR> -T <FLOAT> [-w <INT> -n <STR> -j <STR> -J <FLOAT> -m <STR> -M <FLOAT> -p <INT> -l <INT> -z <INT>] [-h|--help]

    Options:
        -h --help                         show this
        -s, --pop_ids <STR>               User defined samples, e.g. for populations 'A' and 'B'
        -S, --pop_sizes <STR>             List of Ne for all populations
        -T, --split_time <FLOAT>          Split time
        -b, --blocks <INT>                Number of blocks [default: 1]
        -w, --windows <INT>               Number of windows [default: 1]
        -l, --block_length <INT>          Length of a block [default: 256]
        -z, --window_size <INT>           Number of blocks in a window [default: 1]
        -p, --ploidy <INT>                Ploidy of samples [default: 2]
        -r, --mutation_rate <FLOAT>       Mutation rate
        -n, --samples <STR>               Number of samples by populations (same order) [default: 1]
        -j, --join_string <STR>           Newick string of population-join history (Ne changes)
        -J, --ancestral_Ne <FLOAT>        Ne of ancestral population [default: 0]
        -m, --migration_string <STR>      Migration string
        -M, --migration_rate <FLOAT>      Migration rate [default: 0]
        
"""

from docopt import docopt
from timeit import default_timer as timer
#from sys import stderr, exit
import lib.simulate
import lib.gimble
import lib.log


'''

conda install -c conda-forge zarr scikit-allel pandas numpy tqdm docopt parallel more-itertools networkx scipy sagelib msprime sparse
./gIMble simulate -s 'A,B' -b 10 -r 1e-9 -n 2,2 -T 1e5 -S 1000,1000 -J 1000

'''

def main(run_params):
    try:
        start_time = timer()
        args = docopt(__doc__)
        log = lib.log.get_logger(run_params)
        parameterObj = lib.simulate.ParameterObj(args)
        parameterObj.simulate()

        log.info("[*] Total runtime: %.3fs" % (timer() - start_time))
    except KeyboardInterrupt:
        log.info("\n[X] Interrupted by user after %s seconds!\n" % (timer() - start_time))
        exit(-1)

###############################################################################

if __name__ == '__main__':
    main()