#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""usage: gimble setup          [-v <FILE> -b <FILE> -g <FILE> -s <FILE> -o <STR> -p <INT> -D -h]

    [Options]
        -v, --vcf <FILE>                            VCF file of variants
        -b, --bed <FILE>                            BED file of intervals (only defined regions are used)
        -s, --sample <FILE>                         Sample CSV file
        -g, --genome <FILE>                         Genome file
        -o, --outprefix <STR>                       Outprefix
        -d, --simulated_data <FILE>                 simulated_data
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
        self.vcf_file = args['--vcf']
        self.bed_file = args['--bed']
        self.genome_file = args['--genome']
        self.sample_file = args['--sample']
        self.outprefix = args['--outprefix']
        self.pairedness = 2 # once everything is adapted for any size of sample-sets => int(args['--pairedness'])
        self.simulated_data = args['--simulated_data']
        self.check_valid()

    def check(self):
        if self.simulated_data and self.sample_file:
            pass

def main(run_params):
    try:
        start_time = timer()
        args = docopt(__doc__)
        #log = lib.log.get_logger(run_params)
        parameterObj = ParameterObj(args)
        gimbleStore = lib.gimble.create_store(parameterObj)
        # print(gimbleStore.tree())
        # print(gimbleStore.attrs())
        print("[*] Total runtime: %.3fs" % (timer() - start_time))
    except KeyboardInterrupt:
        print("\n[X] Interrupted by user after %s seconds!\n" % (timer() - start_time))
        exit(-1)

###############################################################################

if __name__ == '__main__':
    main()