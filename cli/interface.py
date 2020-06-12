"""
Usage: gimble <module> [<args>...] [-D -V -h]

  [Modules]
    setup                 Setup DataStore
    info                  Print information about DataStore
    blocks                Generate blocks from data in DataStore 
    windows               Generate windows from blocks in DataStore (requires blocks)
    model                 Build demographic model
    simulate              Simulate data [TBI] 
    inference             Make inference [TBI] (requires blocks)
    grid                  Make grid [TBI]
    scan                  Scan using grid [TBI] (requires windows)
    
    
  [Options]
    -h, --help                         Show this screen.
    -D, --debug                        Print debug information.
    -V, --version                      Show version.

  [Dependencies] 
    -------------------------------------------------------------------------------------------------------------------------------------------------------
    | $ conda install -c conda-forge oyaml zarr scikit-allel pandas numpy tqdm docopt parallel more-itertools networkx scipy sagelib networkx pygraphviz  |
    -------------------------------------------------------------------------------------------------------------------------------------------------------

"""

import sys
import os
from docopt import docopt
from timeit import default_timer as timer

def main(gimble_dir):
    try:
        __version__ = '0.5.0'
        version = "gimble v%s" % __version__
        start_time = timer()
        args = docopt(__doc__, version=version, options_first=True)
        if '--version' in args['<args>'] or '-V' in args['<args>']:
            sys.exit(version)
        params = {
            'module': args['<module>'],
            'path': gimble_dir,
            'cwd': os.getcwd(),
            'debug': True if '--debug' in args['<args>'] or '-D' in args['<args>'] else False,
            'version': version
        }
        if args['<module>'] == 'setup':
            import cli.setup as setup
            setup.main(params)
        elif args['<module>'] == 'blocks':
            import cli.blocks as blocks
            blocks.main(params)
        elif args['<module>'] == 'windows':
            import cli.windows as windows
            windows.main(params)
        elif args['<module>'] == 'model':
            import cli.model as model
            model.main(params)
        elif args['<module>'] == 'info':
            import cli.info as info
            info.main(params)
        elif args['<module>'] == 'inference':
            import cli.inference as inference
            inference.main(params)
        elif args['<module>'] == 'simulate':
            import cli.simulate as simulate
            simulate.main(params)
        else:
            sys.exit("%r is not a gimble module. See 'gimble -help'." % args['<module>'])
    except KeyboardInterrupt:
        sys.stderr.write("\n[X] Interrupted by user after %i seconds!\n" % (timer() - start_time))
        sys.exit(-1)