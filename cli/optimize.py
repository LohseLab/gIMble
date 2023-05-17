"""
usage: gimble optimize                      -z <z> -l <l> -d <d> [-w] 
                                            [-y <y>] -m <m> -r <r> -u <u>
                                            [-T <T>] [-M <M>] -A <A> -B <B> [-C <C>] 
                                            [-g <g>] [-p <p>] [-s <s>] [-i <i>] [--xtol <xtol>] [--ftol <ftol>] 
                                            [-e <e>] [-f] [-h|--help]

        -z, --zarr_file=<z>                 Path to existing GimbleStore 
        -l, --optimize_label=<l>            Label used to store grid for later access

    [Data]
        -d, --data_key=<d>                  Data key ('tally/...' or 'simulate/...')
        -w, --windowsum                     Apply optimize to sum of windows in dataset [default: False]

    [Model]
        -m, --model=<m>                     Model name: DIV, MIG_AB, MIG_BA, IM_AB or IM_BA
        -r, --ref_pop=<r>                   Population ID of reference population used for scaling
                                                - A or B or A_B (for models DIV, IM_AB, IM_BA)
                                                - A or B (for models MIG_AB, MIG_BA)
        -u, --mu=<u>                        Mutation rate (in mutations/site/generation)
        -y, --sync_pops=<y>                 Synchronization of Ne parameters during optimization. Optional.
                                                - A,B or A,A_B or B,A_B or A,B,A_B (for models DIV, IM_AB, IM_BA)
                                                - A,B (for models MIG_AB, MIG_BA)

    [Parameters]                            Single floats OR boundaries for optimization in the format [min,max]. 
                                                example 1: --T=100000 for fixed parameter T = 100000
                                                example 2: --Ne_A=10000,20000 for boundaries for optimization 
                                                                of Ne_A between 10000 and 20000

        -A, --Ne_A <A>                      Effective population size of population A (in years)
        -B, --Ne_B <B>                      Effective population size of population B (in years) 
        -C, --Ne_A_B <C>                    Effective population size of ancestral population A_B (in years)
        -T, --T <T>                         Split time (in generations) 
        -M, --me <M>                        Migration rate (per lineage probability of migrating) 
                                                **backwards** in time with direction determined by model name: 
                                                - MIG_AB and IM_AB: A->B 
                                                - MIG_BA and IM_BA: B->A
    [Optimization]
        -g, --algorithm=<g>                 NLOPT optimization algorithm [default: CRS2]
                                                - CRS2
                                                - sbplx
                                                - neldermead
        -p, --processes=<p>                 Number of processes. Only relevant for optimization of windows [default: 1] 
        -s, --start_point=<s>               Point from which to start optimization [default: midpoint]
                                                - 'midpoint' : midpoint between all boundary values
                                                - 'random': based on random seed in INI file
        -i, --max_iterations=<i>            Maximum number of iterations to perform when 
                                                optimizing. Depending on the algorithm --max_iterations can be 
                                                exceeded slightly. Deactivate with -1 [default: 10000]
        --xtol=<xtol>                       Relative tolerance on norm of vector of optimisation parameters.
                                                Float between 0 and 1, deactivate with -1 [default: -1]
        --ftol=<ftol>                       Relative tolerance on lnCL. 
                                                Float between 0 and 1, deactivate with -1 [default: -1]

    [Options]
        -e, --seed=<e>                      Seed used for randomness [default: 19]
        -f, --force                         Force overwrite of existing analysis.
        -h,--help                           show this
"""

'''
--keep_anomalies                    During NLOPT optimizations certain parameter combinations can yield anomalous 
                                    behaviour of the likelihood function. By default, these likelihood-anomalies
                                    are detected and set to -inf to avoid polluting the optimization. This option
                                    will turn this behaviour off.  
'''
from timeit import default_timer as timer
from docopt import docopt
import lib.runargs
import sys
import itertools

class OptimizeParameterObj(lib.runargs.RunArgs):
    '''Sanitises command line arguments and stores parameters.'''

    def __init__(self, params, args):
        super().__init__(params)
        self.zstore = self._get_path(args['--zarr_file'])
        self.optimize_label = args['--optimize_label']
        self.data_key = args['--data_key']
        self.windowsum = args['--windowsum']
        self.Ne_A = self._get_optimize_parameter("Ne_A", args['--Ne_A']) 
        self.Ne_B = self._get_optimize_parameter("Ne_B", args['--Ne_B']) 
        self.Ne_A_B = self._get_optimize_parameter("Ne_A_B", args.get('--Ne_A_B', None)) 
        self.T = self._get_optimize_parameter("T", args.get('--T', None))
        self.me = self._get_optimize_parameter("me", args.get('--me', None))
        self.model = self._check_model(args['--model']) 
        self.sync_pops = self._get_sync_pops(args['--sync_pops'])
        self.ref_pop = self._get_ref_pop(args['--ref_pop'])
        self.mu = self._get_float(args['--mu']) 
        self.processes = self._get_int(args['--processes']) 
        self.seed = self._get_int(args['--seed']) 
        self.overwrite = args['--force']
        self.algorithm = self._get_nlopt_algorithm_name(args['--algorithm'])
        self.start_point_method = self._check_start_point(args['--start_point'])
        self.max_iterations = self._get_int(args['--max_iterations'])
        self.xtol_rel = self._get_float(args['--xtol'])
        self.ftol_rel = self._get_float(args['--ftol'])

    def _get_nlopt_algorithm_name(self, algorithm):
        ALGORITHMS = {'neldermead', 'sbplx', 'CRS2'}
        if algorithm not in ALGORITHMS:
            sys.exit("[X] --algorithm must be one of the following: %s" % (", ".join(list(ALGORITHMS))))
        return algorithm

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

    def _get_sync_pops(self, sync_pops):
        pops_by_model = {
            'DIV': set(['A', 'B', 'A_B']),
            'MIG_AB': set(['A', 'B']),
            'MIG_BA': set(['A', 'B']),
            'IM_AB': set(['A', 'B', 'A_B']),
            'IM_BA': set(['A', 'B', 'A_B'])}
        if sync_pops is None:
            return sync_pops 
        sync_sets = set([frozenset(combination) for combination in itertools.combinations(pops_by_model[self.model], 2)] + [frozenset(pops_by_model[self.model]),])
        split_sync_pops = sync_pops.split(',')
        if not len(split_sync_pops) == len(set(split_sync_pops)):
            sys.exit("[X] --sync_pops has repeated population IDs. %s" % (", ".join(split_sync_pops)))
        sync_pops = frozenset(sync_pops.split(','))
        if not sync_pops in sync_sets:
            sync_sets_string = " or ".join([",".join(sync_set) for sync_set in sorted(sync_sets)])
            sys.exit("[X] --sync_pops for model %r must be one of the following: %s" % (self.model, sync_sets_string)) 
        pop_sizes = ["Ne_%s" % pop_id for pop_id in list(sync_pops)]
        parameters_by_pop_size = {pop_size : str(getattr(self, pop_size)) for pop_size in pop_sizes}
        pop_size_by_parameters = {}
        for pop_size, parameter in parameters_by_pop_size.items(): 
            pop_size_by_parameters.setdefault(parameter, []).append(pop_size) 
        if not len(pop_size_by_parameters) == 1:
            sys.exit("[X] Populations in --sync_pops (%s) must have identical values/boundaries: %s" % (", ".join(sync_pops), ", ".join(pop_sizes))) 
        return pop_sizes

    def _get_optimize_parameter(self, parameter, arg):
        if arg is None:
            return arg
        l = arg.replace("=","").split(",")
        if len(l) == 1:
            try:
                return [float(l[0])]
            except ValueError:
                pass
        elif len(l) == 2:
            _min, _max = l
            try:
                _min_float, _max_float = float(_min), float(_max)
                if _min_float >= _max_float:
                    sys.exit("[X] Invalid boundaries for %r: min=%s and max=%s" % (parameter, _min_float, _max_float))
                return [_min_float, _max_float]
            except ValueError:
                pass
        else:
            pass
        sys.exit("[X] Parameter %r must be a single float OR boundaries for optimization in the format [min,max], not: %s" % (parameter, ",".join(l)))

    def _check_start_point(self, start_point_string):
        if start_point_string in set(['random', 'midpoint']):
            return start_point_string
        sys.exit("[X] '--start_point' must be 'midpoint' or 'random'")

def main(params):
    try:
        start_time = timer()
        args = docopt(__doc__)
        print("[+] Running 'gimble optimize' ...")
        parameterObj = OptimizeParameterObj(params, args)
        import lib.gimble
        gimbleStore = lib.gimble.Store(path=parameterObj.zstore, create=False)
        gimbleStore.optimize(
            optimize_label=parameterObj.optimize_label,
            data_key=parameterObj.data_key,
            windowsum=parameterObj.windowsum,
            Ne_A=parameterObj.Ne_A,
            Ne_B=parameterObj.Ne_B,
            Ne_A_B=parameterObj.Ne_A_B,
            T=parameterObj.T,
            me=parameterObj.me,
            model=parameterObj.model,
            sync_pops=parameterObj.sync_pops,
            ref_pop=parameterObj.ref_pop,
            mu=parameterObj.mu,
            processes=parameterObj.processes,
            seed=parameterObj.seed,
            overwrite=parameterObj.overwrite,
            start_point_method=parameterObj.start_point_method,
            nlopt_maxeval=parameterObj.max_iterations,
            nlopt_xtol_rel=parameterObj.xtol_rel,
            nlopt_ftol_rel=parameterObj.ftol_rel,
            nlopt_algorithm=parameterObj.algorithm,
            )
        print("[*] Total runtime was %s" % (lib.runargs.format_time(timer() - start_time)))
    except KeyboardInterrupt:
        print("\n[X] Interrupted by user after %s !\n" % (lib.runargs.format_time(timer() - start_time)))
        exit(-1)