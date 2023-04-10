"""
usage: gimble parse                      -g <g> -v <v> -b <b> -s <s> [-z <z>] [-f] [-h]

    [Options]
        -g, --genome_f=<g>               Gimble genome file (TSV) of sequence IDs/lengths for filtering BED file.
        -v, --vcf_f=<v>                  VCF file of variants. bgzip'ed. Indexed.
        -b, --bed_f=<b>                  Gimble BED file of regions for filtering VCF file (horizontally).
        -s, --sample_f=<s>               Gimble sample file (CSV) for filtering VCF file (vertically). 
                                            Only two populations are supported.
        -z, --zarr=<z>                   Prefix to use for GimbleStore [default: gimble]
        -f, --force                      Force overwrite of existing data
        -h, --help                       Show this
    
"""

import sys
from timeit import default_timer as timer
from docopt import docopt
import lib.runargs 

class ParseParameterObj(lib.runargs.RunArgs):
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
        print("[+] Running 'gimble parse' ...")
        parameterObj = ParseParameterObj(params, args)
        import lib.gimble
        gimbleStore = lib.gimble.Store(prefix=parameterObj.outprefix, create=True, overwrite=parameterObj.overwrite)
        gimbleStore.measure(
            genome_f=parameterObj.genome_f, 
            sample_f=parameterObj.sample_f, 
            bed_f=parameterObj.bed_f, 
            vcf_f=parameterObj.vcf_f)
        gimbleStore.log_action(module=parameterObj._MODULE, command=parameterObj._get_cmd())
        print("[*] Total runtime was %s" % (lib.runargs.format_time(timer() - start_time)))
    except KeyboardInterrupt:
        print("\n[X] Interrupted by user after %s !\n" % (lib.runargs.format_time(timer() - start_time)))
        exit(-1)