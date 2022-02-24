"""usage: gimble measure                 -v <FILE> -b <FILE> -g <FILE> -s <FILE> [-z <STR> -f -h]

    [Input]
        -g, --genome_f <FILE>            Gimble genome file (TSV) of sequence IDs/lengths for filtering BED file.
        -b, --bed_f <FILE>               Gimble BED file of regions for filtering VCF file (horizontally).
        -s, --sample_f <FILE>            Gimble sample file (CSV) for filtering VCF file (vertically). Only two populations are supported.
        -v, --vcf_f <FILE>               VCF file of variants. bgzip'ed. Indexed.
    
    [Options]
        -z, --zarr <STR>                 Prefix to use for ZARR store [default: gimble]
        -f, --force                      Force overwrite if GimbleStore already exists.

        -h --help                        Show this
    
"""

'''
User needs to modify sample and genomefile and bed file. A good idea is to make a copy of these files.

The genomefile is a way to control which sequences are actually processed within gimble.
It is a good idea to restrict the analysis to chromosome-length sequences as only these will carry enough signal to analyse the demography
Caution with sex chromosomes if any of the samples is haploid for these

# remove regions through subtractions from gimble bed
bedtools subtract -a ../preprocess/heliconius_20220202.bed -b /data/lohse/modelgenomland/analyses/heliconius/9_bed_files/hmel2_5.chromosomes.annotations.bed > heliconius_20220202.non_genic.bed
bedtools subtract -a heliconius_20220202.non_genic.bed -b /data/lohse/modelgenomland/analyses/heliconius/9_bed_files/hmel2_5.chromosomes.repeatmasked.bed > heliconius_20220202.non_genic_non_repeats.bed


'''
import sys

from timeit import default_timer as timer
from docopt import docopt
import lib.gimble 

class ParseParameterObj(lib.gimble.ParameterObj):
    '''Sanitises command line arguments and stores parameters.'''

    def __init__(self, params, args):
        super().__init__(params)
        self.vcf_f = self._get_path(args['--vcf_f'])
        self.bed_f = self._get_path(args['--bed_f'])
        self.genome_f = self._get_path(args['--genome_f'])
        self.sample_f = self._get_path(args['--sample_f'])
        self.outprefix = args['--zarr']
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
        args = docopt(__doc__)
        parameterObj = ParseParameterObj(params, args)
        gimbleStore = lib.gimble.Store(prefix=parameterObj.outprefix, create=True, overwrite=parameterObj.overwrite)
        gimbleStore.measure(
            genome_f=parameterObj.genome_f, 
            sample_f=parameterObj.sample_f, 
            bed_f=parameterObj.bed_f, 
            vcf_f=parameterObj.vcf_f)
        gimbleStore.log_action(module=parameterObj._MODULE, command=parameterObj._get_cmd())
        print("[*] Total runtime was %s" % (lib.gimble.format_time(timer() - start_time)))
    except KeyboardInterrupt:
        print("\n[X] Interrupted by user after %s !\n" % (lib.gimble.format_time(timer() - start_time)))
        exit(-1)