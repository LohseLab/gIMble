#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""usage: gimble setup          [-v <FILE> -b <FILE> -g <FILE> -s <FILE> -o <STR> -f -D -h]

    [Input]
        -g, --genome_f <FILE>        Gimble genome file (TSV) of sequence IDs/lengths for filtering BED file.
        -b, --bed_f <FILE>           Gimble BED file of regions for filtering VCF file (horizontally).
        -s, --sample_f <FILE>        Gimble sample file (CSV) for filtering VCF file (vertically). Only two populations are supported.
        -v, --vcf_f <FILE>           VCF file of variants. bgzip'ed. Indexed.
    
    [Options]
        -o, --outprefix <STR>           Prefix to use in output files [default: gimble]
        -f, --force                     Force overwrite if GimbleStore already exists.
        -D, --debug                     Show debugging information

        -h --help                       Show this
    
"""
import sys

from timeit import default_timer as timer
from docopt import docopt
import lib.gimble 

'''
[To Do]
- rename to: 'init', 'load', ... ?
'''

class SetupParameterObj(lib.gimble.ParameterObj):
    '''Sanitises command line arguments and stores parameters.'''

    def __init__(self, params, args):
        super().__init__(params)
        self.vcf_f = self._get_path(args['--vcf_f'])
        self.bed_f = self._get_path(args['--bed_f'])
        self.genome_f = self._get_path(args['--genome_f'])
        self.sample_f = self._get_path(args['--sample_f'])
        self.outprefix = args['--outprefix']
        self.overwrite = args['--force']
        self._pairedness = 2
        self._check()

    def _check(self):
        required_values_by_arg = {
            '--vcf_f': self.vcf_f,
            '--bed_f': self.bed_f,
            '--genome_f': self.genome_f,
            '--sample_f': self.sample_f
        }
        missing_args = [k for k,v in required_values_by_arg.items() if v is None]
        if missing_args:
            sys.exit("[X] Please provide arguments for %s" % (", ".join(missing_args)))

def main(params):
    try:
        start_time = timer()
        print("[+] Running 'gimble setup'")
        args = docopt(__doc__)
        parameterObj = SetupParameterObj(params, args)
        gimbleStore = lib.gimble.Store(prefix=parameterObj.outprefix, create=True, overwrite=parameterObj.overwrite)
        gimbleStore.setup_seq(parameterObj)
        gimbleStore.info()
        print("[*] Total runtime: %.3fs" % (timer() - start_time))
    except KeyboardInterrupt:
        print("\n[X] Interrupted by user after %s seconds!\n" % (timer() - start_time))
        exit(-1)