"""usage: blocktools windows -s <FILE> -g <FILE> -b <FILE> -v <FILE> [-w <INT> -l <INT> -m <INT> -o <STR> -t <INT> -h]
    
    -s, --sample_file <FILE>                    CSV file ("sample_id,population_id")
    -g, --genome_file <FILE>                    Genome file (as used in BedTools)

    -b, --blocks_h5 <FILE>                         blocks HDF5 file
    -v, --variants_h5 <FILE>                    variants HDF5 file

    -w, --window_size <INT>                     Window size in blocks [default: 500]
    -l, --overlap <INT>                         Window overlap in blocks [default: 50]
    -m, --max_block_distance <INT>              Maximum distance in bases between blocks in same window [default: 10000]

    -t, --threads <INT>                         Number of threads to use [default: 1]
    -o, --prefix <STR>                          Folder/prefix for output
    
    -h, --help

"""

from docopt import docopt
from timeit import default_timer as timer
import lib.windows

def main():
    start_time = timer()
    args = docopt(__doc__)
    parameterObj = lib.windows.task_parse_parameters(args)
    entityCollection = lib.windows.task_generate_entityCollection(parameterObj)
    block_id_batches = lib.windows.task_get_block_regions(parameterObj, entityCollection)
    mutype_counters_by_block_id = lib.windows.task_get_variants_df(parameterObj, entityCollection)
    lib.windows.task_make_windows(parameterObj, entityCollection, block_id_batches)
    lib.windows.task_write_window_output(parameterObj, entityCollection, mutype_counters_by_block_id)
    print("[+] Total runtime: %.3fs" % (timer() - start_time))

if __name__ == "__main__":
    pass