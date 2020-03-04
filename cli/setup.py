#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""usage: gimble setup          -v <FILE> -b <FILE> -g <FILE> -s <FILE> -o <STR> 
                                        [ -p <INT> -D -h]

    Options:
        -v, --vcf <FILE>                            VCF file of variants
        -b, --bed <FILE>                            BED file of intervals (only defined regions are used)
        -s, --sample <FILE>                         Sample CSV file
        -g, --genome <FILE>                         Genome file
        -o, --outprefix <STR>                       Outprefix
        -p, --pairedness <INT>                      Number of samples in sample-sets for intervals/GTs [default: 2]

        -D, --debug                                 Print debug information

        -h --help                                   show this

"""

from docopt import docopt
from timeit import default_timer as timer
#from sys import stderr, exit
import lib.gimble


'''

conda install -c conda-forge zarr scikit-allel pandas numpy tqdm docopt parallel more-itertools networkx scipy sagelib networkx pygraphviz sparse

- GOAL is to store BED/VCF file in minimal arrays by CHROMs
- should allow creation of diagnostic plots

Generates ZARR datastore, with
- Sample/Population information
- Path to VCF file
- Path to BED file

'''
class ParameterObj(object):
    def __init__(self, args):
        print(args)
        self.vcf_file = args['--vcf']
        self.bed_file = args['--bed']
        self.genome_file = args['--genome']
        self.sample_file = args['--sample'] 
        self.outprefix = args['--outprefix']
        self.pairedness = int(args['--pairedness'])


def main(run_params):
    try:
        start_time = timer()
        args = docopt(__doc__)
        #log = lib.log.get_logger(run_params)
        parameterObj = ParameterObj(args)
        store = lib.gimble.create_store(parameterObj)
        print(store.tree())
        print(store.attrs())
        #log.info("[*] Total runtime: %.3fs" % (timer() - start_time))
    except KeyboardInterrupt:
        #log.info("\n[X] Interrupted by user after %s seconds!\n" % (timer() - start_time))
        exit(-1)

###############################################################################

if __name__ == '__main__':
    main()