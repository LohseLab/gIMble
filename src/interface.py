"""
Usage: ./gIMble <module> [<args>...]

Modules:
    graph           Generate a state graph
    probs           Analyse state graph

Options:
    -h, --help                         Show this screen.
    -v, --version                      Show version.

Help:
    https://gimble.readme.io/ (TBD)
"""

import sys
from docopt import docopt
from timeit import default_timer as timer

def main():
    try:
        __version__ = '0.1.0'
        start_time = timer()
        args = docopt(__doc__, version=__version__, options_first=True)
        if args['<module>'] == 'graph':
            import src.graph as graph
            graph.main()
        elif args['<module>'] == 'probs':
            import src.probs as probs
            probs.main()
        else:
            sys.exit("%r is not a gIMble module. See 'gIMble -help'." % args['<module>'])
    except KeyboardInterrupt:
        sys.stderr.write("\n[X] Interrupted by user after %i seconds!\n" % (timer() - start_time))
        sys.exit(-1)
