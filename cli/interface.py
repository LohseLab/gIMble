"""usage: gimble <module> [<args>...] [-V -h]

  [Input]
    preprocess            Install gimbleprep instead
    parse                 Parse files into GimbleStore
    blocks                Generate blocks from parsed data in GimbleStore (requires 'parse')
    windows               Generate windows from blocks in GimbleStore (requires 'blocks')
    tally                 Tally variation for inference (requires 'blocks' or 'windows')

  [Simulation]
    simulate              Simulate data based on specific parameters or gridsearch results  
    
  [Inference]
    optimize              Perform global parameter optimisation on tally/simulation
    makegrid              Precalculate grid of parameters
    gridsearch            Evaluate tally/simulation against a precomputed grid (requires 'makegrid')

  [Info]
    info                  Print metrics about data in GimbleStore
    list                  List information saved in GimbleStore
    query                 Extract information from GimbleStore
    delete                Delete information in GimbleStore

  [Experimental]
    partitioncds          Partition CDS sites in BED file by degeneracy in sample GTs

  [Options]
    -h, --help            Show this screen
    -V, --version         Show version
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
    'list': 'cli.list',
    'delete': 'cli.delete',
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

installation_steps = """[========================= Missing dependencies =========================]
1. Get conda from https://conda.io/miniconda.html

2. Create the following conda environment 
>>> conda create -n gimble python=3.7.12 numpy docopt tqdm pandas tabulate zarr scikit-allel parallel matplotlib msprime demes dask numcodecs python-newick nlopt -c conda-forge -c bioconda -y

3. Load the environment (needs to be activated when using gimble)
>>> conda activate gimble

4. Install agemo (make sure you have the conda environment activated)
>>> (gimble) pip install agemo

5. Rock'n'roll ...
[========================================================================]
"""
def main(gimble_dir=None):
    if gimble_dir is None:
        gimble_dir = os.path.dirname(os.path.join(os.path.realpath(__file__), '..'))
    try:
        start_time = timer()
        __version__ = '1.0.3'
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
            try:
                runner = importlib.import_module(RUNNER_BY_MODULE[params['module']])
                runner.main(params)
            except ImportError as error:
                print("[X] ImportError: %s" % error)
                print(installation_steps)
    except KeyboardInterrupt:
        sys.stderr.write("\n[X] Interrupted by user after %i seconds!\n" % (timer() - start_time))
        sys.exit(-1)
