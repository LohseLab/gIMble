#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""usage: gIMble model -s <STR> [-n <STR> -j <STR> -m <STR> -p <INT> -c --nograph --nomodel --noyaml] [-h|--help]

    Options:
        -h --help                         show this
        -p, --ploidy <INT>                Ploidy of samples [default: 2]
        -s, --pop_ids <STR>               User defined populations (single letters) [default: 'A,B']
        -n, --samples <STR>               Number of samples by population (same order) [default: '1,1']
        -j, --join_string <STR>           Newick string of population-join history (Ne changes) [default: '']
        -m, --migration_string <STR>      Migration string, e.g. 'A>B' or 'B<A' [default: '']
        --nomodel                         No model output
        --nograph                         No graph output 
        --noyaml                          No yaml output 
"""

from docopt import docopt
from timeit import default_timer as timer
from sys import stderr, exit
import lib.model
from lib.gimble import RunObj
import ast
import re
import collections
import sys
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

class ParameterObj(RunObj):
    '''Sanitises command line arguments and stores parameters'''

    def __init__(self, params, args):
        super().__init__(params)
        '''
        2. write rules regarding model strings
        - Population IDs must be single letters, e.g. A and B
        - 
        3. write config YAML
        '''
        self.ploidy = self._get_int(args['--ploidy'])
        self.pop_ids = sorted(args['--pop_ids'].split(','))
        self.samples_by_pop_id = {pop_id: int(samples) for pop_id, samples in zip(args['--pop_ids'].split(','), args['--samples'].split(','))}
        self.join_events = self._parse_join_events(args['--join_string']) if args['--join_string'] else []
        self.migration_events = self._parse_migration_events(args['--migration_string'])
        self.sample_stateObj = self._get_sample_stateObj()
        self.set_of_events = self._get_set_of_events() 
        self.out_prefix = self._get_outprefix()
        self.model_output = False if args['--nomodel'] else True 
        self.graph_output = False if args['--nograph'] else True
        self.yaml_output = False if args['--noyaml'] else True
        print(self.__dict__) 

    def _get_outprefix(self):
        pop_string = "_".join(self.pop_ids)
        sample_string = "_".join([str(self.samples_by_pop_id[pop_id]) for pop_id in self.pop_ids])
        prefix = ["gimble.model.s_%s.n_%s" % (pop_string, sample_string)]
        if self.migration_events:
            prefix.append(".".join(self.migration_events))
        if self.join_events:
            prefix.append("J_%s" % ".".join(self.join_events))
        return ".".join(prefix)
            
        
    def _get_ploidy(ploidy):
        return int(ploidy)
        
    def _get_set_of_events(self):
        return set(["C_%s" % pop_id for pop_id in self.pop_ids + self.join_events] + self.join_events + self.migration_events)

    def _get_sample_stateObj(self):
        pop_dict = {}
        for pop_id in self.pop_ids:
            pop_list = []
            for i in range(self.samples_by_pop_id[pop_id]):
                for j in range(self.ploidy):
                    pop_list.append("%s%s" % (pop_id.lower(), i+1))
            pop_dict[pop_id.upper()] = lib.model.PopObj(pop_list)
        return lib.model.StateObj(pop_dict)

    def _parse_join_events(self, join_string):
        char_counter = collections.Counter(join_string)
        print(char_counter)
        if not all([char_counter[pop_id] <= 1 for pop_id in self.pop_ids]) == True:
            return "[X] Duplicated population IDs"
        #print(lib.model.joins(ast.literal_eval(re.sub(r'([a-zA-Z]+)',r'"\1"', join_string))))
        return lib.model.joins(ast.literal_eval(re.sub(r'([a-zA-Z]+)',r'"\1"', join_string)))

    def _parse_migration_events(self, migration_string):
        events = []
        print(migration_string)
        if not migration_string == "''":
            for m in migration_string.split(','):
               if '>' in m:
                   a, b = m.split('>')
                   events.append("M_%s_%s" % (a, b))
               elif '<' in m:
                   a, b = m.split('<')
                   events.append("M_%s_%s" % (b, a))
               else:
                   sys.exit("[X] %r is not a valid migration parameter.")
        return events

def main(params):
    try:
        start_time = timer()
        print("[+] Running 'gimble model'")
        args = docopt(__doc__)
        parameterObj = ParameterObj(params, args)
        stateGraph = lib.model.graph_generator(parameterObj)
        if parameterObj.model_output:
            stateGraph.write_model(parameterObj)
        if parameterObj.graph_output:
            stateGraph.plot_dot_graph(parameterObj)
    except KeyboardInterrupt:
        stderr.write("\n[X] Interrupted by user after %s seconds!\n" % (timer() - start_time))
        exit(-1)

###############################################################################

if __name__ == '__main__':
    main()