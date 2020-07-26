#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""usage: gIMble simulate -v <FILE> -b <FILE> -g <FILE> -s <FILE> -o <STR> [-D -h]

    Options:
        -v, --vcf <FILE>                            VCF file of variants
        -b, --bed <FILE>                            BED file of intervals (only defined regions are used)
        -s, --sample <FILE>                         Sample CSV file
        -g, --genome <FILE>                         Genome file
        -o, --outprefix <STR>                       Outprefix
        -D, --debug                                 Print debug information
        -h --help                                   show this

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