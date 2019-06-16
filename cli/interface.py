"""
Usage: gIMble <module> [<args>...]

Modules:
    blocks                Makes blocks
    variants              Fetches and analyses variants for blocks 
    modify                Modifies/filters blocks/variants
    windows               Constructs windows of blocks
    likelihood            Infer likelihood for data given model 
    gridsearch            TBE

Options:
    -h, --help                         Show this screen.
    -v, --version                      Show version.

"""
import sys
from docopt import docopt
from timeit import default_timer as timer

def main():
    try:
        __version__ = '0.4.0'
        start_time = timer()
        args = docopt(__doc__, version=__version__, options_first=True)
        if args['<module>'] == 'blocks':
            import cli.blocks as blocks
            blocks.main()
        elif args['<module>'] == 'variants':
            import cli.variants as variants
            variants.main()
        elif args['<module>'] == 'modify':
            import cli.modify as modify
            modify.main()
        elif args['<module>'] == 'windows':
            import cli.windows as windows
            windows.main()
        elif args['<module>'] == 'likelihood':
            import cli.likelihood as likelihood
            likelihood.main()
        elif args['<module>'] == 'gridsearch':
            import cli.gridsearch as gridsearch
            gridsearch.main()
        else:
            sys.exit("%r is not a gIMble module. See 'gIMble -help'." % args['<module>'])
    except KeyboardInterrupt:
        sys.stderr.write("\n[X] Interrupted by user after %i seconds!\n" % (timer() - start_time))
        sys.exit(-1)
