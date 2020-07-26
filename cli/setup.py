#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""usage: gimble setup          [-v <FILE> -b <FILE> -g <FILE> -s <FILE> -o <STR> -D -h]

    [Input]
        -g, --genome_file <FILE>        Genome file (TSV) of sequence IDs for filtering BED file. 
                                          Lists the names of the genomic sequences and their length.
        -b, --bed_file <FILE>           BED file for filtering VCF file.  
                                          Lists BedTools multiinter intervals from mosdepth 
                                          'Callable regions' output based on BAM files used in 
                                          VCF file. Fifth column must contain sample 
                                          names as in VCF file.  
        -s, --sample_file <FILE>        Sample file (CSV) for filtering VCF file. 
                                          Lists sample names as in VCF file and their populations.
                                          Only two populations are supported.                                              
        -v, --vcf_file <FILE>           VCF file of variants. Indexed.
    
    [Options]
        -o, --outprefix <STR>           Prefix to use in output files [default: gimble]
        -h --help                       show this
    
"""
from timeit import default_timer as timer
import sys
from docopt import docopt
from lib.gimble import RunObj
import lib.gimble 

class ParameterObj(RunObj):
    '''Sanitises command line arguments and stores parameters.'''

    def __init__(self, params, args):
        super().__init__(params)
        self.vcf_file = self._get_path(args['--vcf_file'])
        self.bed_file = self._get_path(args['--bed_file'])
        self.genome_file = self._get_path(args['--genome_file'])
        self.sample_file = self._get_path(args['--sample_file'])
        self.outprefix = args['--outprefix']
        self._pairedness = 2
        self._check()

    def _check(self):
        required_values_by_arg = {
            '--vcf_file': self.vcf_file,
            '--bed_file': self.bed_file,
            '--genome_file': self.genome_file,
            '--sample_file': self.sample_file
        }
        missing_args = [k for k,v in required_values_by_arg.items() if v is None]
        if missing_args:
            sys.exit("[X] Please provide arguments for %s" % (", ".join(missing_args)))


def main(params):
    try:
        start_time = timer()
        print("[+] Running 'gimble setup'")
        args = docopt(__doc__)
        parameterObj = ParameterObj(params, args)
        gimbleStore = lib.gimble.create_store(parameterObj)
        gimbleStore.info()
        print("[*] Total runtime: %.3fs" % (timer() - start_time))
    except KeyboardInterrupt:
        print("\n[X] Interrupted by user after %s seconds!\n" % (timer() - start_time))
        exit(-1)