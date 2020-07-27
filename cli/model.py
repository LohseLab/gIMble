#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""usage: gIMble model -o <STR> -s <STR> [-n <STR> -j <STR> -m <STR> -p <INT> -c --nograph --nomodel] [-h|--help]

    Options:
        -h --help                         show this
        -p, --ploidy <INT>                Ploidy of samples [default: 2]
        -s, --pop_ids <STR>               User defined samples, e.g. for populations 'A' and 'B'
        -n, --samples <STR>               Number of samples by populations (same order) [default: 1]
        -j, --join_string <STR>           Newick string of population-join history (Ne changes)
        -m, --migration_string <STR>      Migration string
        -c, --complete_labelling          Complete labelling of StateGraph (as opposed to partial)
        -o, --out_prefix <STR>            Prefix for output files
        --nomodel                         No model output
        --nograph                         No graph output 
"""

from docopt import docopt
from timeit import default_timer as timer
from sys import stderr, exit
import lib.model

'''
conda install -c conda-forge networkx matplotlib docopt tqdm pandas numpy psutil pygraphviz

./gIMble model -s A,B -p 2 -n 1,1 -m 'A>B' -j '(A,B)' -o graph.s_A_B.p2.n_1_1.m_AtoB.j_A_B.post_fix.new ; \
./gIMble model -s A,B -p 2 -n 1,1 -m 'B>A' -j '(A,B)' -o graph.s_A_B.p2.n_1_1.m_BtoA.j_A_B.post_fix.new ; \
./gIMble model -s A,B -p 2 -n 2,0 -m 'A>B' -j '(A,B)' -o graph.s_A_B.p2.n_2_0.m_AtoB.j_A_B.post_fix.new ; \
./gIMble model -s A,B -p 2 -n 2,0 -m 'B>A' -j '(A,B)' -o graph.s_A_B.p2.n_2_0.m_BtoA.j_A_B.post_fix.new ; \
./gIMble model -s A,B -p 2 -n 0,2 -m 'A>B' -j '(A,B)' -o graph.s_A_B.p2.n_0_2.m_AtoB.j_A_B.post_fix.new ; \
./gIMble model -s A,B -p 2 -n 0,2 -m 'B>A' -j '(A,B)' -o graph.s_A_B.p2.n_0_2.m_BtoA.j_A_B.post_fix.new

[ To Do ]
- Output config file for simulation/inference based on model
- add header to model file : -s A,B -p 2 -n 1,1

'''


def main(run_params):
    try:
        main_time = timer()
        args = docopt(__doc__)
        parameterObj = lib.model.ParameterObj(args)
        stateGraph = lib.model.graph_generator(parameterObj)
        if parameterObj.model_output:
            stateGraph.write_model(parameterObj)
        if parameterObj.graph_output:
            stateGraph.plot_dot_graph(parameterObj)
    except KeyboardInterrupt:
        stderr.write("\n[X] Interrupted by user after %s seconds!\n" % (timer() - main_time))
        exit(-1)

###############################################################################

if __name__ == '__main__':
    main()