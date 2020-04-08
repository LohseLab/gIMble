"""
Usage: gIMble <module> [<args>...] [-D -V -h]

Modules:

    setup                 Setup DataStore
    blocks                Make blocks
    windows               Make windows
    model                 Build new model
    inference             Make inference
    
Options:
    -h, --help                         Show this screen.
    -D, --debug                        Print debug information
    -V, --version                      Show version

"""
import sys
from docopt import docopt
from timeit import default_timer as timer

def main():
    try:
        __version__ = '0.5.0'
        start_time = timer()
        args = docopt(__doc__, version=__version__, options_first=True)
        # exit if --version
        if '--version' in args['<args>'] or '-V' in args['<args>']:
            sys.exit("gIMble v%s" % __version__)
        # setup run_params
        run_params = {
            'module': args['<module>'],
            'debug': True if '--debug' in args['<args>'] or '-D' in args['<args>'] else False,
            'version': __version__
        }
        if args['<module>'] == 'setup':
            import cli.setup as setup
            setup.main(run_params)
        elif args['<module>'] == 'blocks':
            import cli.blocks as blocks
            blocks.main(run_params)
        elif args['<module>'] == 'windows':
            import cli.windows as windows
            windows.main(run_params)
        elif args['<module>'] == 'model':
            import cli.model as model
            model.main(run_params)
        elif args['<module>'] == 'inference':
            import cli.inference as inference
            inference.main(run_params)
        elif args['<module>'] == 'simulate':
            import cli.simulate as simulate
            simulate.main(run_params)
        
        elif args['<module>'] == 'varfilter':
            import cli.varfilter as varfilter
            varfilter.main(run_params)
        elif args['<module>'] == 'windows':
            import cli.windows as windows
            windows.main(run_params)
        elif args['<module>'] == 'likelihood':
            import cli.likelihood as likelihood
            likelihood.main(run_params)
        elif args['<module>'] == 'estimate':
            import cli.estimate as estimate
            estimate.main(run_params)
        elif args['<module>'] == 'gridsearch':
            import cli.gridsearch as gridsearch
            gridsearch.main(run_params)
        elif args['<module>'] == 'makegrid':
            import cli.makegrid as makegrid
            makegrid.main(run_params)
        else:
            sys.exit("%r is not a gIMble module. See 'gIMble -help'." % args['<module>'])
    except KeyboardInterrupt:
        sys.stderr.write("\n[X] Interrupted by user after %i seconds!\n" % (timer() - start_time))
        sys.exit(-1)
