"""
usage: gimble makegrid                  (-z <FILE> | -o <STR>) -l <STR> -m <STR> -r <STR> -u <FLOAT> \
                                        -b <INT> -A <VALUES> -B <VALUES> [-C <VALUES>] [-T <VALUES>] [-M <VALUES>] 
                                        [-k <STR>] [-n <INT> -s <INT>] [-f] [-h|--help]
                                                 
        -z, --zarr_file <FILE>              Path to existing GimbleStore 
        -o, --outprefix <STR>               Prefix to use for GimbleStore [default: gimble]
        -l, --makegrid_label <STR>          Label used to store grid for later access
        
    [Model]
        -m, --model <STR>                   Model name: DIV, MIG_AB, MIG_BA, IM_AB or IM_BA
        -b, --block_length <INT>            Block length MUST match block length of data to be searched
        -r, --ref_pop <STR>                 Population ID of reference population used for scaling
                                            - A or B or A_B (for models DIV, IM_AB, IM_BA)
                                            - A or B (for models MIG_AB, MIG_BA)
        -u, --mu <FLOAT>                    Mutation rate (in mutations/site/generation)
        -k, --kmax <STR>                    Max count per mutation type beyond which counts 
                                            are treated as marginals. Order of mutation 
                                            types is (hetB, hetA, hetAB, fixed)
                                            [default: 2,2,2,2]

    [Grid parameters]                       single float or distribution in the format [min,max,number_steps,lin|log]
                                                example 1: --T=100000 for T = [100000]
                                                example 2: --Ne_A=10000,20000,3,lin for Ne_A = [10000, 15000, 20000]
                                                example 3: --me=0.0,1e-4,5,log for me = [0.0, 1e-7, 1e-6, 1e-5, 1e-4]

        -A, --Ne_A <VALUES>                 Effective population size of population A (in years)
        -B, --Ne_B <VALUES>                 Effective population size of population B (in years) 
        -C, --Ne_A_B <VALUES>               Effective population size of ancestral population A_B (in years)
        -T, --T <VALUES>                    Split time (in generations) 
        -M, --me <VALUES>                   Migration rate (per lineage probability of migrating) 
                                            **backwards** in time with direction determined by model name: 
                                            - MIG_AB and IM_AB: A->B 
                                            - MIG_BA and IM_BA: B->A
    [Options]
        -n, --processes <INT>               Number of processes [default: 1] 
        -s, --seed <INT>                    Seed used for randomness [default: 19]
        -f, --force                         Force overwrite of existing grid in GimbleStore
        -h --help                           show this
"""

from timeit import default_timer as timer
from docopt import docopt
import lib.gimble
import sys

class MakeGridParameterObj(lib.gimble.ParameterObj):
    '''Sanitises command line arguments and stores parameters.'''

    def __init__(self, params, args):
        super().__init__(params)
        self.zstore = self._get_path(args["--zarr_file"])
        self.prefix = self._get_prefix(args["--outprefix"])
        self.Ne_A = self._get_makegrid_parameters("Ne_A", args['--Ne_A']) 
        self.Ne_B = self._get_makegrid_parameters("Ne_B", args['--Ne_B']) 
        self.Ne_A_B = self._get_makegrid_parameters("Ne_A_B", args.get('--Ne_A_B', None)) 
        self.T = self._get_makegrid_parameters("T", args.get('--T', None))
        self.me = self._get_makegrid_parameters("me", args.get('--me', None))
        self.makegrid_label = args['--makegrid_label']
        self.model = self._check_model(args['--model']) 
        self.block_length = self._get_int(args['--block_length']) 
        self.ref_pop = self._get_ref_pop(args['--ref_pop'])
        self.mu = self._get_float(args['--mu']) 
        self.kmax = self._check_kmax(args['--kmax']) 
        self.processes = self._get_int(args['--processes']) 
        self.seed = self._get_int(args['--seed']) 
        self.overwrite = args['--force']
    
    def _get_ref_pop(self, ref_pop):
        pops_by_model = {
            'DIV': set(['A', 'B', 'A_B']),
            'MIG_AB': set(['A', 'B']),
            'MIG_BA': set(['A', 'B']),
            'IM_AB': set(['A', 'B', 'A_B']),
            'IM_BA': set(['A', 'B', 'A_B'])}
        if not ref_pop in pops_by_model[self.model]:
            sys.exit("[X] --ref_pop for model %r must be one of the following: %s" % (self.model, ", ".join(pops_by_model[self.model])))
        return "Ne_%s" % ref_pop

    def _get_makegrid_parameters(self, parameter, arg):
        if arg is None:
            return arg
        l = arg.split(",")
        if len(l) == 1:
            try:
                return [float(l[0])]
            except ValueError:
                pass
        elif len(l) == 4:
            _min, _max, _steps, _distr = l
            try:
                _min_float, _max_float, _steps_int, _distr = float(_min), float(_max), int(_steps), _distr
                if _steps_int < 2:
                    sys.exit("[X] Parameter %r must have 'number_steps' >= 2." % (parameter))
                if _distr in set(["lin", "log"]) and _min_float < _max_float:
                    return [_min_float, _max_float, _steps_int, _distr]
            except ValueError:
                pass
        else:
            pass
        sys.exit("[X] Parameter %r must be a single float or distribution in the format [min,max,number_steps,lin|log], not: %s" (parameter, ",".join(l)))

def main(params):
    try:
        start_time = timer()
        args = docopt(__doc__)
        parameterObj = MakeGridParameterObj(params, args)
        gimbleStore = lib.gimble.Store(
            path=parameterObj.zstore, 
            prefix=parameterObj.prefix, 
            create=(False if parameterObj.zstore else True))
        #gimbleStore.makegrid(
        #    config=parameterObj.config,
        #    num_cores=parameterObj.num_cores,
        #    overwrite=parameterObj.overwrite,
        #    )
        gimbleStore.makegrid(
            Ne_A=parameterObj.Ne_A,
            Ne_B=parameterObj.Ne_B,
            Ne_A_B=parameterObj.Ne_A_B,
            T=parameterObj.T,
            me=parameterObj.me,
            makegrid_label=parameterObj.makegrid_label,
            model=parameterObj.model,
            block_length=parameterObj.block_length,
            ref_pop=parameterObj.ref_pop,
            mu=parameterObj.mu,
            kmax=parameterObj.kmax,
            processes=parameterObj.processes,
            seed=parameterObj.seed,
            overwrite=parameterObj.overwrite,
           )
        print("[*] Total runtime was %s" % (lib.gimble.format_time(timer() - start_time)))
    except KeyboardInterrupt:
        print("\n[X] Interrupted by user after %s !\n" % (lib.gimble.format_time(timer() - start_time)))
        exit(-1)