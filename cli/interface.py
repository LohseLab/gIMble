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
    makegrid_legacy       Precalculate grid of parameters (legacy)
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
import importlib
from docopt import docopt
from timeit import default_timer as timer

RUNNER_BY_MODULE = {
    'preprocess': 'cli.preprocess',   
    'parse': 'cli.parse',        
    'blocks': 'cli.blocks',       
    'windows': 'cli.windows',      
    'tally': 'cli.tally',
    'simulate': 'cli.simulate',
    'optimize': 'cli.optimize',
    'optimize_legacy': 'cli.optimize_legacy',
    'makegrid': 'cli.makegrid',
    'makegrid_legacy': 'cli.makegrid_legacy',
    'gridsearch': 'cli.gridsearch',
    'info': 'cli.info',
    'query': 'cli.query',
    'partitioncds': 'cli.partitioncds',
}
MODULES = RUNNER_BY_MODULE.keys()

def main(gimble_dir):
    try:
        start_time = timer()
        __version__ = '0.8.1'
        version = "gimble v%s" % __version__
        args = docopt(__doc__, version=version, options_first=True)
        if '--version' in args['<args>'] or '-V' in args['<args>']:
            sys.exit("gimble v%s" % __version__)
        params = {
            'module': args['<module>'],
            'path': gimble_dir,
            'cwd': os.getcwd(),
            'version': version
        }
        if not params['module'] in MODULES:
            print("[X] %r is not a gimble module.\n" % params['module'])
            sys.exit(__doc__)
        else:
            runner = importlib.import_module(RUNNER_BY_MODULE[params['module']])
            runner.main(params)
    except KeyboardInterrupt:
        sys.stderr.write("\n[X] Interrupted by user after %i seconds!\n" % (timer() - start_time))
        sys.exit(-1)