#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""usage: gIMble variants -s <FILE> -g <FILE> -v <FILE> -b <FILE> [-o <STR> -t <INT> -h]
    
    -s, --sample_file <FILE>                    CSV file ("sample_id,population_id")
    -g, --genome_file <FILE>                    Genome file (as used in BedTools)

    -v, --vcf <FILE>                            VCF file (tabix'ed)
    -b, --blocks <FILE>                         blocks HDF5 file 
    
    -t, --threads <INT>                         Number of threads to use [default: 1]
    -o, --prefix <STR>                          Folder/prefix for output
    
    -h, --help

"""

from docopt import docopt
from timeit import default_timer as timer
import lib.variants

def main():
    start_time = timer()
    args = docopt(__doc__)
    print("[#] ### gIMble VARIANTS ###")
    parameterObj = lib.variants.task_parse_parameters(args)
    entityCollection = lib.variants.task_generate_entityCollection(parameterObj)
    lib.variants.task_load_blockObjs(parameterObj, entityCollection)
    lib.variants.task_infer_mutypes(parameterObj, entityCollection)
    lib.variants.task_write_variant_output(parameterObj, entityCollection)
    print("[+] Total runtime: %.3fs" % (timer() - start_time))

if __name__ == "__main__":
    pass