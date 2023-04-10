"""
usage: gimble info                 -z <z> [-h|--help]
                                            
    [Options]
        -z, --zarr_file=<z>          Path to existing GimbleStore 
        -h --help                    Show this
"""

from timeit import default_timer as timer
from docopt import docopt
import lib.runargs

class InfoParameterObj(lib.runargs.RunArgs):
    '''Sanitises command line arguments and stores parameters.'''

    def __init__(self, params, args):
        super().__init__(params)
        self.zstore = self._get_path(args['--zarr_file'])
        self.tree = False

'''[To Do]
- MUST write info report to file

should not print shape, but kmax (and explanation what it means)
                                                                                       kmax*
    └── 'tally/blocks_kmax2' .....................................................  [2,2,2,2]
            Het_A=0.01583 Het_B=0.01553 D_xy=0.02204 F_st=0.16862 Marginality=23.46%
[+] └── 'tally/blocks_kmax3' .....................................................  [2,3,2,3]
            Het_A=0.01583 Het_B=0.01553 D_xy=0.02204 F_st=0.16862 Marginality=11.13%

    * kmax is bla ... 

info should say which parts have no data! (as opposed to just not listing anything)
- BED intervals (has to say BED otherwise people get confused)
'''
def main(params):
    try:
        start_time = timer()
        args = docopt(__doc__)
        print("[+] Running 'gimble info' ...")
        parameterObj = InfoParameterObj(params, args)
        import lib.gimble
        gimbleStore = lib.gimble.Store(path=parameterObj.zstore)
        print("[+] Getting report. This might take a while ...")
        info = gimbleStore.info(
            version=parameterObj._VERSION,
            tree=parameterObj.tree)
        print(info)
        print("[*] Total runtime was %s" % (lib.runargs.format_time(timer() - start_time)))
    except KeyboardInterrupt:
        print("\n[X] Interrupted by user after %s !\n" % (lib.runargs.format_time(timer() - start_time)))
        exit(-1)