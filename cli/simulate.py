#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""usage: gIMble simulate -o <STR> -s <STR> [-n <STR> -j <STR> -m <STR> -p <INT> -c --nograph --nomodel] [-h|--help]

    Options:
        -h --help                         show this
        -p, --ploidy <INT>                Ploidy of samples [default: 2]
        -s, --pop_ids <STR>               User defined samples, e.g. for populations 'A' and 'B'
        -n, --samples <STR>               Number of samples by populations (same order) [default: 1]
        -j, --join_string <STR>           Newick string of population-join history (Ne changes)
        -m, --migration_string <STR>      Migration string
        -c, --complete_labelling          Complete labelling of StateGraph (as opposed to partial)
        -o, --out_prefix <STR>            Prefix for output files
        --nomodel                         No model output
        --nograph                         No graph output 

"""

from docopt import docopt
from timeit import default_timer as timer
#from sys import stderr, exit
import lib.setup
import lib.gimble
import lib.log


'''

conda install -c conda-forge zarr scikit-allel pandas numpy tqdm docopt parallel more-itertools networkx scipy sagelib msprime sparse

'''

def main(run_params):
    try:
        start_time = timer()
        args = docopt(__doc__)
        print(args)
        log = lib.log.get_logger(run_params)
        parameterObj = lib.setup.ParameterObj(args)
        store = lib.classes.Store(parameterObj)
        print(store.tree())
        log.info("[*] Total runtime: %.3fs" % (timer() - start_time))
    except KeyboardInterrupt:
        log.info("\n[X] Interrupted by user after %s seconds!\n" % (timer() - start_time))
        exit(-1)

###############################################################################

if __name__ == '__main__':
    main()