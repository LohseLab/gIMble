"""
usage: gimble gridsearch               -z <z> -g <g> -d <d> [-w] [-p <p> -c <c> -f] [-h|--help]
                                                                                                                   
    [Options]
        -z, --zarr_file=<z>            Path to existing GimbleStore
        -g, --grid_key=<g>             Makegrid key run in GimbleStore
        -d, --data_key=<d>             Dataset key ('tally/...' or 'simulate/...')
        -w, --windowsum                Sum windows in dataset [default: False]
        -p, --processes=<p>            Number of processes [default: 1]
        -c, --chunksize=<c>            Size of chunks to use in parallelisation 
                                           (greater chunksize => greater RAM requirements)
                                           [default: 500]
        -f, --overwrite                Overwrite results in GimbleStore
        -h, --help                     Show this
"""
from timeit import default_timer as timer
from docopt import docopt
import lib.runargs

class GridsearchParameterObj(lib.runargs.RunArgs):
    def __init__(self, params, args):
        super().__init__(params)
        self.zstore = self._get_path(args['--zarr_file'])
        self.data_key = args['--data_key']
        self.windowsum = args['--windowsum']
        self.grid_key = args['--grid_key']
        self.overwrite = args['--overwrite']
        self.num_cores = self._get_int(args['--processes'])    # number of cores for independent processes
        self.chunksize = self._get_int(args['--chunksize'])    # size of chunks in first dimension of tally/grid dask array 

def main(params):
    try:
        start_time = timer()
        args = docopt(__doc__)
        print("[+] Running 'gimble gridsearch' ...")
        parameterObj = GridsearchParameterObj(params, args)
        import lib.gimble
        gimbleStore = lib.gimble.Store(path=parameterObj.zstore, create=False)
        gimbleStore.gridsearch(
            data_key=parameterObj.data_key,
            grid_key=parameterObj.grid_key,
            windowsum=parameterObj.windowsum,
            num_cores=parameterObj.num_cores,
            chunksize=parameterObj.chunksize,
            overwrite=parameterObj.overwrite,
            )
        print("[*] Total runtime was %s" % lib.runargs.format_time(timer() - start_time))
    except KeyboardInterrupt:
        print("\n[X] Interrupted by user after %s !\n" % lib.runargs.format_time(timer() - start_time))
        exit(-1)