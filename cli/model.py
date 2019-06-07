#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""usage: gIMble model -o STR -s STR [-p INT] [-M] [-e -E] [-P -G] [-h|--help]

    Options:
        -h --help                       show this
        -p, --ploidy INT                Ploidy of samples [default: 2]
        -s, --samples STR               User defined samples, e.g. for populations 'a' and 'b':
                                            - [['b'], {'a'}] : 2 pops, migration ({a} -> [b])
                                            - [{'b'}, ['a']] : 2 pops, migration ({b} -> [a])
        -M, --migration                 Allow migration between populations 
        -e, --mass_migration            Allow exodus between populations from {}->[]
        -E, --mass_migration_reverse    Allow exodus between populations from []->{}
        -o, --out_prefix STR            Prefix for output files
"""

from docopt import docopt
from timeit import default_timer as timer
from sys import stderr, exit
import lib.model

'''
Problems:
    Establishment of population_ids -> C_ancestor/C_derived
    For output: label 'E' to 'bigL'
    Output which pop is ancestor/derived
'''


def main():
    try:
        main_time = timer()
        args = docopt(__doc__)
        parameterObj = lib.model.ParameterObj(args)
        parameterObj.setup()
        state_graph, stateObj_by_idx = lib.model.build_state_graph(parameterObj)
        header, paths = lib.model.get_state_paths(state_graph, parameterObj)
        parameterObj.write_paths(header, paths)
        #parameterObj.plot_graph(state_graph, parameterObj)
        #parameterObj.write_graph(state_graph)
        #state_graph = parameterObj.read_graph(graph_file)
        
    except KeyboardInterrupt:
        stderr.write("\n[X] Interrupted by user after %s seconds!\n" % (timer() - main_time))
        exit(-1)

###############################################################################

if __name__ == '__main__':
    main()