"""usage: gIMble blocks -s <FILE> -g <FILE> -b <FILE> [-n <INT> -x <INT> -l <INT> -m <INT> -o <STR> -t <INT> -h -d]
    
    -s, --sample_file <FILE>                    CSV file ("sample_id,population_id")
    -g, --genome_file <FILE>                    Genome file (as used in BedTools)

    -b, --bed_file <FILE>                       multiintersectBED file of 

    -x, --block_length <INT>                    Block length in bases [default: 64]
    -l, --max_block_span <INT>                  Maximum span of block (from start to end) [default: 80]
    -m, --min_interval_length <INT>             Minimum length in bases of a BED interval to be considered [default: 0]
    -n, --min_samples <INT>                     Minimum number of samples per population in interval [default: 1]
    -t, --threads <INT>                         Number of threads to use [default: 1]
    -o, --prefix <STR>                          Folder/prefix for output
    
    -h, --help
    -d, --debug                                 Write debugging logs
"""

from timeit import default_timer as timer
from docopt import docopt
import lib.blocks
import lib.gimblelog

def main(run_params):
    start_time = timer()
    args = docopt(__doc__)
    print(args)
    log = lib.gimblelog.get_logger(run_params, args['--debug'])
    blockParameterObj = lib.blocks.task_generate_parameterObj(args)
    entityCollection = lib.blocks.task_generate_entityCollection(blockParameterObj)
    region_dfs = lib.blocks.task_generate_region_dfs(blockParameterObj, entityCollection)
    lib.blocks.task_make_blocks(blockParameterObj, region_dfs, entityCollection)
    lib.blocks.task_write_block_output(blockParameterObj, entityCollection)
    log.info("[*] Total runtime: %.3fs" % (timer() - start_time))

if __name__ == "__main__":
    pass