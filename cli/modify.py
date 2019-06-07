"""usage: 
            gIMble modify blocks    -s <FILE> -g <FILE> -b <FILE> -c <FILE> [-t <INT> -o <STR> -h]
            gIMble modify variants  -v <FILE>             [-m <INT> -M <INT> -t <INT> -o <STR> -h]
    

    -s, --sample_file <FILE>                    CSV file ("sample_id,population_id")
    -g, --genome_file <FILE>                    Genome file (as used in BedTools)
    -b, --blocks FILE                           blocks HDF5 file
    -c, --coordinates FILE                      Coordinates file
    
    -v, --variants_h5 FILE                      variants HDF5 file
    -m, --max_missing INT                       max missing genotypes per pair [default: 4]
    -M, --max_multiallelic INT                  max multiallelic genotypes per pair [default: 4]

    -t, --threads INT                           Number of threads to use [default: 1]
    -o, --prefix STR                            Folder for output
    
    -h, --help


"""

from docopt import docopt
from timeit import default_timer as timer
import lib.modify

def main():
    start_time = timer()
    args = docopt(__doc__)
    if args['blocks']:
        parameterObj = lib.modify.task_parse_parameters(args, mode='blocks')
        entityCollection = lib.modify.task_generate_entityCollection(parameterObj)
        coordinateTransformObj = lib.modify.task_parse_coordinates(parameterObj)
        lib.modify.task_load_blockObjs(parameterObj, entityCollection)
        lib.modify.task_transform_coordinates(parameterObj, entityCollection, coordinateTransformObj)
        lib.modify.task_write_modify_output(parameterObj, entityCollection)
    elif args['variants']:
        parameterObj = lib.modify.task_parse_parameters(args, mode='variants')
        lib.modify.task_filter_variants(parameterObj)
    
    print("[+] Total runtime: %.3fs" % (timer() - start_time))

if __name__ == "__main__":
    pass