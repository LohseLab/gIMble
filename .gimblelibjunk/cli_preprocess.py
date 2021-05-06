#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""usage: gimbl preprocess                 -v FILE -f FILE -b DIR [-g INT -m INT -M INT -t INT -o STR -k] [-h|--help]
                                            
    [Options]
        -v, --vcf_file FILE                         VCF file (raw, unfiltered)
        -f, --fasta_file FILE                       FASTA file (reference used in VCF)
        -b, --bam_dir FILE                          Directory containing BAM files (only those used in VCF)
        -g, --snpgap INT                            SnpGap [default: 2]
        -m, --min_depth INT                         Min read depth [default: 8]
        -M, --max_depth INT                         Max read depth (as multiple of SD from mean) [default: 2]
        -t, --threads INT                           Threads [default: 1]
        -o, --outprefix STR                         Outprefix [default: gimble]
        -k, --keep_tmp                              Do not delete temporary files [default: False]
        -h, --help                                  Show this
"""

from timeit import default_timer as timer
from docopt import docopt

import tempfile
import pathlib
import shutil

import gimble.parameters
import gimble.preprocess
'''
[To Do]
- diagnostic plots 
    - plotting gimble.bed
    - plotting gimble.vcf
        - missingness (across non-overlapping 10kb windows)
        - multiallelicity (across non-overlapping 10kb windows)
        - snp-rate (across non-overlapping 10kb windows)
        - pca of snps
        - pairwise distance (diagonal heatmap: https://seaborn.pydata.org/examples/many_pairwise_correlations.html)
            - how to calculate?

'''

class PreprocessParameterObj(gimble.parameters.ParameterObj):
    '''Sanitises command line arguments and stores parameters'''
    def __init__(self, params, args):
        super().__init__(params)
        self.fasta_file = self._get_path(args['--fasta_file'], path=True)
        self.vcf_file = self._get_path(args['--vcf_file'], path=True)
        self.bam_dir = self._get_path(args['--bam_dir'], path=True)
        self.bam_files = list(self.bam_dir.glob('*.bam'))
        self.snpgap = self._get_int(args['--snpgap'])
        self.min_depth = self._get_int(args['--min_depth'])
        self.max_depth = self._get_float(args['--max_depth'])
        self.tmp_dir = tempfile.mkdtemp(prefix='.tmp_gimble_', dir=".")
        self.threads = self._get_int(args['--threads'])
        self.outprefix = args['--outprefix']
        self.keep_tmp = args['--keep_tmp']
        self.gimble_genome_file = pathlib.Path('%s.genomefile' % self.outprefix)
        self.gimble_coverage_summary_file = pathlib.Path('%s.coverage_summary.csv' % self.outprefix)
        self.gimble_vcf_file = pathlib.Path('%s.vcf.gz' % self.outprefix)
        self.gimble_sample_file = pathlib.Path("%s.samples.csv" % self.outprefix)
        self.gimble_log_file = pathlib.Path("%s.log.txt" % self.outprefix)
        self.bed_multi_callable_file =  self.tmp_dir / pathlib.Path('%s.multi_callable.bed' % self.outprefix)
        self.coverage_data = []
        self.commands = []

    def clean_up(self):        
        if self.keep_tmp:
            print("[+] Not deleting Temporary folder %r since '--keep_tmp' was specified." % self.tmp_dir)
        else:
            print("[+] Deleting temporary folder %r ..." % self.tmp_dir)
            shutil.rmtree(self.tmp_dir)

    def write_log(self):
        with open(self.gimble_log_file, 'w') as fh:
            fh.write("%s\n" % "\n".join(self.commands))

def main(params):
    try:
        start_time = timer()
        args = docopt(__doc__)
        parameterObj = PreprocessParameterObj(params, args)
        gimble.preprocess.preprocess(parameterObj)
        print("[*] Total runtime: %.3fs" % (timer() - start_time))
    except KeyboardInterrupt:
        print("\n[X] Interrupted by user after %s seconds!\n" % (timer() - start_time))
        exit(-1)
