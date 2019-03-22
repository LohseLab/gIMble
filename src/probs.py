#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""usage: gIMble probs -o STR -g FILE [-M FLOAT] [-E FLOAT] [-h|--help]

    Options:
        -h --help                       show this
        -g, --graph FILE                State graph to analyse
        -o, --out_prefix STR            Prefix for output files
        -M, --migration_rate FLOAT      Migration rate [default: 0.5]
        -E, --exodus_rate FLOAT         Exodus rate [default: 0.75]
        
"""

from docopt import docopt
from collections import Counter, defaultdict
import networkx as nx
import matplotlib as mat
mat.use("agg")
from timeit import default_timer as timer
from sys import stderr, exit
from ast import literal_eval
from pandas import DataFrame

'''
[Problems]
- Ipython notebook does not work on my machine, make general implementation of task
- sage needs python2 !?!
    - hacky compiling with python3 : https://wiki.sagemath.org/Python3-compatible%20code
    - official version soon: 
        - https://trac.sagemath.org/ticket/26212
        - https://trac.sagemath.org/ticket/15530

[To do]
- change path parser to work on source/sink √
- Datastructure to know which edges have a given sub-state => allows knowing paths on which mutations can be placed
'''

################################### Functions #################################

def create_csv(out_f, header, rows, sep):
    df = DataFrame(rows, columns=header)
    print(df)
    df.to_csv(out_f, index=False, sep=sep)

def pairwise(iterable):
    it = iter(iterable)
    a = next(it, None)
    for b in it:
        yield (a, b)
        a = b

class StatePathObj(object):
    def __init__(self, idx, nodes):
        self.idx = idx 
        self.path = [] # node of each step

def get_state_paths(state_graph):
    '''
    Networkx documentation : https://networkx.github.io/documentation/stable/index.html
    '''
    start_time = timer()
    # out_edges_counter = Counter([source for source, sink, count in state_graph.edges.data('count') for i in range(count)]) # previous approach

    # Get total outgoing edges by event_type for each node 
    event_types = {
        '[C]' : 'C', '{C}' : 'C', 'C' : 'C',
        'M' : 'M', 'R' : 'R', 'E' : 'E'
    }
    total_by_event_by_node_idx = defaultdict(Counter)
    for source, sink, data in state_graph.edges.data():
        event = event_types[data['event']]
        total_by_event_by_node_idx[source][event] += data['count']
    
    # Convert to DoD so that can be JSON'ed
    total_by_event_by_node_idx = {node_idx: dict(total_by_event) for node_idx, total_by_event in total_by_event_by_node_idx.items()}
    print(total_by_event_by_node_idx)
    # find source / sink
    node_idx_by_meta = {meta: node_idx for node_idx, meta in nx.get_node_attributes(state_graph, 'meta').items() if meta}

    header = ['idx', 'path', 'probability', 'events']
    paths = []
    for idx, path in enumerate(nx.all_simple_paths(state_graph, source=node_idx_by_meta['source'], target=node_idx_by_meta['sink'])):
        #print(">> PATH %s %s" % (idx, path))
        count_sum, total_sum = 1, 1
        events = []
        for source, sink in pairwise(path):
            event = event_types[state_graph[source][sink]['event']]
            count_sum *= state_graph[source][sink]['count']
            #print("%s -> %s" % (source, sink), state_graph[source][sink])
            events.append(state_graph[source][sink]['event'])
            total_sum *= total_by_event_by_node_idx[source][event]
            #print("[E] %s : p = %s" % (state_graph[source][sink]['event'], state_graph[source][sink]['count'] / out_edges_counter[source]))
        #print("[P] P = %s" % (count_sum / total_sum))
        paths.append([idx, "->".join([str(node) for node in path]), count_sum / total_sum, "->".join(events)])
    print("[+] Found %s paths (sum(P)=%s) in %s seconds." % (len(paths), sum([path[2] for path in paths]), timer() - start_time))
    return (header, paths)

###############################################################################

################################### Classes ###################################

class ParameterObj(object):
    def __init__(self, args):
        self.graph_file = args['--graph']
        self.out_prefix = args['--out_prefix']
        
    def write_paths(self, header, paths):
        start_time = timer()
        out_f = "%s.paths.txt" % (self.get_basename())
        create_csv(out_f, header, paths, sep="\t")
        print("[+] Written paths into file %s in %s seconds." % (out_f, timer() - start_time))

    def write_graph(self, state_graph):
        start_time = timer()
        out_f = "%s.gml" % (self.get_basename())
        # states have to be stringified ...
        for idx, node in state_graph.nodes(data=True):
            # print(idx, node)
            node['state'] = str(node.get('state', None))
        nx.write_gml(state_graph, out_f)
        print("[+] Written graph into file %s in %s seconds." % (out_f, timer() - start_time))
        return out_f

    def read_graph(self):
        start_time = timer()
        state_graph = nx.read_gml(self.graph_file)
        for idx, node in state_graph.nodes(data=True):
            node['state'] = literal_eval(node['state'])
            if 'pos' in node:
                node['pos'] = tuple(node['pos'])
        print("[+] Read graph from file %s in %s seconds." % (self.graph_file, timer() - start_time))
        return state_graph

###############################################################################

################################### Main ######################################

def main():
    try:
        main_time = timer()
        args = docopt(__doc__)
        print(args)
        parameterObj = ParameterObj(args)
        state_graph = parameterObj.read_graph()
        get_state_paths(state_graph)
    except KeyboardInterrupt:
        stderr.write("\n[X] Interrupted by user after %s seconds!\n" % (timer() - main_time))
        exit(-1)
        
###############################################################################

if __name__ == '__main__':
    main()