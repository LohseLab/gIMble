"""
usage: gimble simulate                          (-z <z> | -o <o>) (-s <s>) (-a <a> -b <b>)
                                                [-r <r> -w <w> -n <n> -l <l>] -u <u> [-k <k>] [--continuous_genome]
                                                ((-m <m> -A <A> -B <B> [-C <C> -T <T> -M <M>]) | (-g <g> [-t <t>]))
                                                [--rec_rate <rate> | --rec_map <bed>]
                                                [-e <e> -p <p> -f] [-h|--help]
                                            
        -z, --zarr=<z>                          Path to zarr store
        -o, --outprefix=<o>                     Prefix to make new zarr store
        -s, --simulate_label=<s>                Label used to store simulation for later access
        -p, --processes=<p>                     Processes [default: 1]
        -e, --seed=<e>                          Seed used for randomness [default: 19]
        -f, --force                             Force overwrite of existing analysis
        -h --help                               Show this

    [Simulation]                                
        -a, --samples_A=<a>                     Number of diploid individuals in population A 
        -b, --samples_B=<b>                     Number of diploid individuals in population B 
        -r, --replicates=<r>                    Number of replicates [default: 1]
        -w, --windows=<w>                       Number of windows per replicate [default: 1]
        -n, --blocks=<n>                        Number of blocks per window 
                                                    n blocks will yield n*(a*b) pair-blocks
                                                    [default: 1]
        -l, --block_length=<l>                  Number of sites per block [default: 64]
        --continuous_genome                     By default, mutations are simulated under a finite-site 
                                                mutation model. Use this option to use an infinite-site 
                                                mutation model [default: False]
        -u, --mu=<u>                            Mutation rate (in mutations/site/generation)
        -k, --kmax=<k>                          Max count per mutation type beyond which counts 
                                                are treated as marginals. Order of mutation 
                                                types is (hetB, hetA, hetAB, fixed)
                                                [default: 2,2,2,2]

    [Model based simulation]
        -m, --model=<m>                         Model name: DIV, MIG_AB, MIG_BA, IM_AB or IM_BA
        -A, --Ne_A=<A>                          Effective population size of population A (in years)
        -B, --Ne_B=<B>                          Effective population size of population B (in years)
        -C, --Ne_AB=<C>                         Effective population size of ancestral population A_B (in years)
        -T, --T=<T>                             Split time (in generations)                     
        -M, --me=<M>                            Migration rate (per lineage probability of migrating) 
                                                **backwards** in time with direction determined by model name: 
                                                - MIG_AB and IM_AB: A->B 
                                                - MIG_BA and IM_BA: B->A

    [Gridsearch based simulation]
        -g, --gridsearch_key=<g>                Key of gridsearch stored in gimble store ('gridsearch/...'). 
                                                If specified, simulation will parameterize each window with 
                                                the parameter combination which obtained the best composite 
                                                log-likelihood (lnCL) during gridsearch. 
                                                This is relevant for bootstrapping approaches. Number of 
                                                windows in gridsearch result must match '--windows' 
                                                (and '--rec_map', if specified).
        -t, --constraint=<t>                    This option sets constraints under which lnCL are evaluated. 
                                                    example: '--constraint me=2.21e-06'
        
    [Recombination]
        --rec_rate=<rate>                       Constant recombination rate per site (in cM/Mb) [default: 0.0]
        --rec_map=<bed>                         Recombination map in BED format w/ recombination rate 
                                                in 4th column. Number of rows needs to match '--windows'
"""

from timeit import default_timer as timer
from docopt import docopt
import lib.runargs 
import sys              

class SimulateParameterObj(lib.runargs.RunArgs):
    """Sanitises command line arguments and stores parameters."""

    def __init__(self, params, args):
        super().__init__(params)
        self.zstore = self._get_path(args["--zarr"])
        self.prefix = self._get_prefix(args["--outprefix"])
        self.processes = self._get_int(args["--processes"])
        self.simulate_label = self._check_label(args['--simulate_label'])
        self.overwrite = args['--force']
        self.seed = self._get_int(args['--seed'])
        self.samples_A = self._get_int(args['--samples_A'])
        self.samples_B = self._get_int(args['--samples_B'])
        self.replicates = self._get_int(args['--replicates'])
        self.windows = self._get_int(args['--windows'])
        self.blocks = self._get_int(args['--blocks'])
        self.block_length = self._get_int(args['--block_length'])
        self.continuous_genome = args['--continuous_genome']
        self.kmax = self._check_kmax(args['--kmax'])
        self.Ne_A = self._get_float(args['--Ne_A'], ret_none=True)
        self.Ne_B = self._get_float(args['--Ne_B'], ret_none=True)
        self.Ne_A_B = self._get_float(args['--Ne_AB'], ret_none=True)
        self.T = self._get_float(args['--T'], ret_none=True)
        self.me = self._get_float(args['--me'], ret_none=True)
        self.mu = self._get_float(args['--mu'], ret_none=True)
        self.model = self._check_model(args['--model'], ret_none=True)
        self.gridsearch_key = args['--gridsearch_key']
        self.constraint = self._get_constraint(args['--constraint'])
        self.rec_rate = self._get_float(args['--rec_rate'])
        self.rec_map = self._get_path(args['--rec_map'])

    def _check_label(self, label):
        invalid_chars = set([c for c in label if not c.isalnum() and not c in set([".", "-", "_"])])
        if invalid_chars:
            sys.exit("[X] --simulate_label contains invalid characters (%r). Should only contain alphanumericals and -_." % "".join(invalid_chars))
        return label

    def _get_constraint(self, constraint):
        if not constraint:
            return {}
        try:
            return {element.split("=")[0]: float(element.split("=")[1]) for element in constraint.split(",") if element}
        except (ValueError, IndexError):
            sys.exit("[X] '--constraint %s' is not in the right format." % constraint)

def main(params):
    try:
        start_time = timer()
        args = docopt(__doc__)
        print("[+] Running 'gimble simulate' ...")
        parameterObj = SimulateParameterObj(params, args)
        import lib.gimble
        gimble_store = lib.gimble.Store(
            path=parameterObj.zstore, 
            prefix=parameterObj.prefix, 
            create=(True if parameterObj.prefix else False))
        gimble_store.simulate(
            zstore=parameterObj.zstore,
            prefix=parameterObj.prefix,
            processes=parameterObj.processes,
            simulate_label=parameterObj.simulate_label,
            overwrite=parameterObj.overwrite,
            seed=parameterObj.seed,
            samples_A=parameterObj.samples_A,
            samples_B=parameterObj.samples_B,
            replicates=parameterObj.replicates,
            windows=parameterObj.windows,
            blocks=parameterObj.blocks,
            block_length=parameterObj.block_length,
            continuous_genome=parameterObj.continuous_genome,
            kmax=parameterObj.kmax,
            Ne_A=parameterObj.Ne_A,
            Ne_B=parameterObj.Ne_B,
            Ne_A_B=parameterObj.Ne_A_B,
            T=parameterObj.T,
            me=parameterObj.me,
            mu=parameterObj.mu,
            model=parameterObj.model,
            gridsearch_key=parameterObj.gridsearch_key,
            constraint=parameterObj.constraint,
            rec_rate=parameterObj.rec_rate,
            rec_map=parameterObj.rec_map
            )
        print("[*] Total runtime: %s" % lib.runargs.format_time(timer() - start_time))
    except KeyboardInterrupt:
        print("\n[X] Interrupted by user after %s !\n" % lib.runargs.format_time(timer() - start_time))
        exit(-1)
