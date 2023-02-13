"""usage: gimble <module> [<args>...] [-V -h]

  [Input]
    preprocess            Preprocess input files
    parse                 Parse files into GimbleStore
    blocks                Generate blocks from parsed data in GimbleStore (requires 'parse')
    windows               Generate windows from blocks in GimbleStore (requires 'blocks')
    tally                 Tally variation for inference (requires 'blocks' or 'windows')

  [Simulation]
    simulate              Simulate data  
    
  [Inference]
    optimize              Perform global parameter optimisation on tally/simulation
    makegrid              Precalculate grid of parameters
    gridsearch            Evaluate tally/simulation against a precomputed grid (requires 'makegrid')

  [Info]
    info                  Print information about GimbleStore
    query                 Extract information from GimbleStore

  [Misc]
    partitioncds          Partition CDS sites in BED file by degeneracy in sample GTs

  [Options]
    -h, --help            Show this screen
    -V, --version         Show version

  [Dependencies] 
    
    ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    | $ conda install bedtools bcftools samtools vcflib mosdepth pysam numpy docopt tqdm pandas tabulate zarr scikit-allel parallel more-itertools networkx giac sagelib matplotlib msprime networkx pygraphviz sympy cerberus maxima -c conda-forge -c bioconda |
    ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
"""

import sys
import os
from docopt import docopt
from timeit import default_timer as timer

def main(gimble_dir):
    try:
        __version__ = '0.7.1'
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
        if args['<module>'] == 'measure':
            import cli.measure as measure
            measure.main(params)
        elif args['<module>'] == 'preprocess':
            import cli.preprocess as preprocess
            preprocess.main(params)
        elif args['<module>'] == 'partitioncds':
            import cli.partitioncds as partitioncds
            partitioncds.main(params)
        elif args['<module>'] == 'plotbed':
            import cli.plotbed as plotbed
            plotbed.main(params)
        elif args['<module>'] == 'query':
            import cli.query as query
            query.main(params)
        elif args['<module>'] == 'tally':
            import cli.tally as tally
            tally.main(params)
        elif args['<module>'] == 'blocks':
            import cli.blocks as blocks
            blocks.main(params)
        elif args['<module>'] == 'windows':
            import cli.windows as windows
            windows.main(params)
        elif args['<module>'] == 'makeconfig':
            import cli.makeconfig as makeconfig
            makeconfig.main(params)
        elif args['<module>'] == 'info':
            import cli.info as info
            info.main(params)
        elif args['<module>'] == 'simulate':
            import cli.simulate as simulate
            simulate.main(params)
        elif args['<module>'] == 'gridsearch':
            import cli.gridsearch as gridsearch
            gridsearch.main(params)
        elif args['<module>'] == 'makegrid':
            import cli.makegrid as makegrid
            makegrid.main(params)
        elif args['<module>'] == 'optimize' or args['<module>'] == 'optimise':
            import cli.optimize as optimize
            optimize.main(params)
        else:
            sys.exit("%r is not a gimble module. See 'gimble -help'." % args['<module>'])
    except KeyboardInterrupt:
        sys.stderr.write("\n[X] Interrupted by user after %i seconds!\n" % (timer() - start_time))
        sys.exit(-1)