"""
usage: gimble windows                    -z <z> [-w <w> -s <s>] [-f] [-h]

    [Options]
        -z, --zarr_file=<z>              Path to existing GimbleStore
        -w, --blocks=<w>                 Number of blocks in windows [default: 500]
        -s, --steps=<s>                  Number of steps (blocks) by which windows are shifted [default: 50]
        -f, --force                      Force overwrite of existing data
        -h, --help                       show this

"""

from docopt import docopt
from timeit import default_timer as timer
import sys
import lib.runargs

class WindowsParameterObj(lib.runargs.RunArgs):
    '''Sanitises command line arguments and stores parameters'''
    def __init__(self, params, args):
        super().__init__(params)
        self.zstore = self._get_path(args['--zarr_file'])
        self.window_size = self._get_int((args['--blocks']))
        self.window_step = self._get_int((args['--steps']))
        self.overwrite = args['--force']
        self.check_block_steps()

    def check_block_steps(self):
        if not self.window_size % self.window_step == 0:
            sys.exit("[X] Quotient of '--blocks' and '--steps' must be an integer. It is %s/%s=%s" % (self.window_size, self.window_step, self.window_size/self.window_step))


def main(params):
    try:
        start_time = timer()
        args = docopt(__doc__)
        print("[+] Running 'gimble windows' ...")
        parameterObj = WindowsParameterObj(params, args)
        import lib.gimble
        gimbleStore = lib.gimble.Store(path=parameterObj.zstore, create=False)
        gimbleStore.windows(
            window_size=parameterObj.window_size, 
            window_step=parameterObj.window_step, 
            overwrite=parameterObj.overwrite)
        gimbleStore.log_action(module=parameterObj._MODULE, command=parameterObj._get_cmd())
        print("[*] Total runtime was %s" % (lib.runargs.format_time(timer() - start_time)))
    except KeyboardInterrupt:
        print("\n[X] Interrupted by user after %s !\n" % (lib.runargs.format_time(timer() - start_time)))
        exit(-1)