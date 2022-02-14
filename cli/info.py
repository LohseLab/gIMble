"""usage: gimbl info                  -z FILE [-h|--help]
                                            
    Options:
        -z, --zarr_file FILE                        ZARR datastore

        -h --help                                   show this
"""

from timeit import default_timer as timer
from docopt import docopt
import lib.gimble

class InfoParameterObj(lib.gimble.ParameterObj):
    '''Sanitises command line arguments and stores parameters.'''

    def __init__(self, params, args):
        super().__init__(params)
        self.zstore = self._get_path(args['--zarr_file'])
        self.tree = False

'''[To Do]
- MUST write info report to file

info should say which parts have no data! (as opposed to just not listing anything)
- BED intervals (has to say BED otherwise people get confused)
'''
def main(params):
    try:
        start_time = timer()
        args = docopt(__doc__)
        parameterObj = InfoParameterObj(params, args)
        gimbleStore = lib.gimble.Store(path=parameterObj.zstore)
        print("[+] Getting report. This might take a while ...")
        info = gimbleStore.info(
            version=parameterObj._VERSION,
            tree=parameterObj.tree)
        print(info)
        print("[*] Total runtime was %s" % (lib.gimble.format_time(timer() - start_time)))
    except KeyboardInterrupt:
        print("\n[X] Interrupted by user after %s !\n" % (lib.gimble.format_time(timer() - start_time)))
        exit(-1)