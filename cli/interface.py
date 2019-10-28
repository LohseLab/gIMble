"""
Usage: gIMble <module> [<args>...]

Modules:

    blocks                Makes blocks
    variants              Fetches and analyses variants for blocks 
    varfilter             Filter variants
    portblocks            Port blocks to new coordinate system
    
    windows               Constructs windows of blocks
    winfilter             Filter windows
    
    likelihood            Infer likelihood for data given model and parameters
    lsearch               Estimate parameters for data given a model 
    
    makegrid              Precompute grid
    gridsearch            Perform gridsearch on precomputed grid
    

Options:
    -h, --help                         Show this screen.
    -v, --version                      Show version.

"""
import sys
from docopt import docopt
from timeit import default_timer as timer

def main():
    try:
        __version__ = '0.5.0'
        start_time = timer()
        args = docopt(__doc__, version=__version__, options_first=True)
        run_params = {
            'module': args['<module>'],
            'version': __version__
        }
        if args['<module>'] == 'blocks':
            import cli.blocks as blocks
            blocks.main(run_params)
        elif args['<module>'] == 'variants':
            import cli.variants as variants
            variants.main(run_params)
        elif args['<module>'] == 'modify':
            import cli.modify as modify
            modify.main(run_params)
        elif args['<module>'] == 'portblocks':
            import cli.portblocks as portblocks
            portblocks.main(run_params)
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
        elif args['<module>'] == 'model':
            import cli.model as model
            model.main(run_params)
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
