#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""usage: gIMble model -o STR -s STR [-j STR -m STR -p INT -c --nograph --nomodel] [-h|--help]

    Options:
        -h --help                       show this
        -p, --ploidy INT                Ploidy of samples [default: 2]
        -s, --pop_ids STR               User defined samples, e.g. for populations 'A' and 'B'
        -j, --join_string STR           Newick string of population-join history (Ne changes)
        -m, --migration_string STR      Migration string
        -c, --complete_labelling        Complete labelling of StateGraph (as opposed to partial)
        -o, --out_prefix STR            Prefix for output files
        --nomodel                       No model output
        --nograph                       No graph output 
"""

from tqdm import tqdm
from docopt import docopt
from timeit import default_timer as timer
from sys import stderr, exit
import itertools
import collections
import networkx as nx
import ast
import re
import sys
import lib.model
import lib.functions

'''
conda install -c conda-forge networkx matplotlib docopt tqdm pandas numpy psutil pygraphviz

./gIMble model -s A,B -p 4 -m 'A>B' -j '(A,B)' -o 'samples_A_B.migration_A>B.join_AB' --nograph
[+] StateGraph
[+] Nodes = 394
[+] Edges = 1693
[+] Origins = 2 [A&B=[a+a+a+a+b+b+b+b], B=[a+a+a+a+b+b+b+b]]
[+] Calculating model ...
[>] Created: 'samples_A_B.migration_A>B.join_AB.model.tsv'
[+] Paths
[+] 733026 => all
[+] 418375 => A&B=[a+a+a+a+b+b+b+b]
[+] 314651 => B=[a+a+a+a+b+b+b+b]
... model file is 330Mb

./gIMble model -s A,B,C -p 2 -m 'A>B,A&B>C' -j '((A,B),C)' -o 'samples_A_B.migration_A>B.join_AB'
[+] StateGraph
[+] Nodes = 326
[+] Edges = 1248
[+] Origins = 3 [C=[a+a+b+b+c+c], A&B&C=[a+a+b+b+c+c], B=[a+a+b+b];C=[c+c]]
[+] Calculating model ...
[%] Processed paths : : 168k [00:14, 11.3k/s]
[>] Created: 'samples_A_B.migration_A>B.join_AB.model.tsv'
[+] Paths
[+] 168250 => all
[+] 96178 => A&B&C=[a+a+b+b+c+c]
[+] 122 => B=[a+a+b+b];C=[c+c]
[+] 71950 => C=[a+a+b+b+c+c]
[+] Plotting graph ...

... model file is 83Mb

'''

MUTYPE_BY_LINEAGE = {
    ('aa') : 'fixed',
    ('bb') : 'fixed',
    ('a') : 'hetA',
    ('abb') : 'hetA',
    ('aab') : 'hetB',
    ('b') : 'hetB',
    ('ab') : 'hetAB'
    }

EDGE_COLOUR_BY_EVENT = {
    'C': 'blue'
}

############################## JOIN input



###################################

################################### Migration input

'''
parse_join_events('C,D,(A,B)')          => ['A&B', 'A&B&C&D']
parse_join_events('(C,D),(A,B)')        => ['C&D', 'A&B', 'A&B&C&D']
parse_join_events('A,((C,D),B)')        => ['C&D', 'B&C&D', 'A&B&C&D']
parse_join_events('A,B')                => ['A&B']
parse_join_events('(A,B),C')            => ['A&B', 'A&B&C']

parse_migration_events('B>A')           => ['B>A']
'B>A,A>B')       => ['B>A', 'A>B']
'A>B,B>A,C>A&B') => ['A>B', 'B>A', 'C>A&B']
parse_migration_events('B>A')           => ['B>A']
'B>A,A>B')       => ['B>A', 'A>B']
'A>B,B>A,C>A&B') => ['A>B', 'B>A', 'C>A&B']
parse_migration_events('B>A')           => ['B>A']
'B>A,A>B')       => ['B>A', 'A>B']
'A>B,B>A,C>A&B') => ['A>B', 'B>A', 'C>A&B']
'''



def pairwise(iterable):
    it = iter(iterable)
    a = next(it, None)
    for b in it:
        yield (a, b)
        a = b

###################################

class PopObj(object):
    def __init__(self, pop=[]):
        self.list = sorted(pop)
        self.string = "[%s]" % ",".join("'%s'" % allele for allele in sorted(pop))
    
    def __str__(self):
        return self.string

    def __iter__(self):
        for allele in self.list:
            yield allele

    def __repr__(self):
        return self.string

    def __len__(self):
        return len(self.list)

    def __hash__(self):
        return hash(self.string)

    def __getitem__(self, i):
        return self.list[i]

    def __add__(self, other):
        return PopObj(sorted(self.list + other.list))

    def __radd__(self, other):
        if other == 0:
            return self
        else:
            return self.__add__(other)

class StateObj(object):
    def __init__(self, value):
        self.popObj_by_pop_id = self._get_dict(value)
        self.label = self._get_label(value)
        self.count_by_ancestor_labels_by_event = collections.defaultdict(lambda: collections.Counter())
        self.lineage_counter = self._get_lineage_counter()
    
    def _get_lineage_counter(self):
        lineage_counter = []
        for popObj in self.popObj_by_pop_id.values():
            for lineage in popObj.list:
                lineage_counter.append(lineage.replace("+", ""))            
        return collections.Counter(lineage_counter)

    def _get_label(self, value):
        if isinstance(value, dict):
            return ";".join(["%s=%s" % (pop_id, popObj) for pop_id, popObj in sorted(self.popObj_by_pop_id.items())])
        elif isinstance(value, str):
            return value
        else:
            sys.exit("[X] value must be 'str' or 'dict', is: %r" % value)
    
    def _get_dict(self, value):
        if isinstance(value, dict):
            if all([isinstance(popObj, PopObj) for pop_id, popObj in value.items()]) == True:
                return value
            else:
                sys.exit("[X] 'dict' must contain PopObjs, contains: %r" % str([type(popObj) for pop_id, popObj in value.items()]))
        elif isinstance(value, str):
            string = value.replace("[", "").replace("]", "").replace("'", "")
            popObj_by_pop_id = {}
            if ";" in string:
                for pop in string.split(";"):
                    pop_id, list_string = pop.split("=")
                    popObj_by_pop_id[pop_id] = PopObj(list_string.split(","))
            else:
                pop_id, list_string = string.split("=")
                popObj_by_pop_id[pop_id] = PopObj(list_string.split(","))
            return popObj_by_pop_id
        else:
            sys.exit("[X] value must be 'str' or 'dict', is: %r" % value)

    def __str__(self):
        return self.label

    def __repr__(self):
        return self.label

    def __hash__(self):
        return hash(self.label)

    def __eq__(self, other):
        return self.label == other.label

    def __len__(self):
        return len(self.popObj_by_pop_id)

    def get_ancestors(self, parameterObj):
        for event in parameterObj.set_of_events:
            if event.startswith("C"):
                self._get_coalescence_ancestors(event)
            elif event.startswith("M"):
                self._get_migration_ancestors(event)
            else:
                self._get_join_ancestors(event)

    def _get_join_ancestors(self, event):
        #print("# join...")
        for pop_1_id, pop_2_id in itertools.combinations(sorted(self.popObj_by_pop_id.keys()), 2):
            _event = "&".join(sorted([pop_1_id, pop_2_id]))
            if _event == event: 
                ancestor_popObj_by_pop_id = {pop_id: popObj for pop_id, popObj in self.popObj_by_pop_id.items() if pop_id not in set([pop_1_id, pop_2_id])}
                ancestor_popObj_by_pop_id[event] = self.popObj_by_pop_id[pop_1_id] + self.popObj_by_pop_id[pop_2_id] # summation of PopObjs
                ancestor_stateObj = StateObj(ancestor_popObj_by_pop_id)
                self.count_by_ancestor_labels_by_event[event][ancestor_stateObj.label] += 1
    
    def _get_coalescence_ancestors(self, event):
        #print("# coalesce...")
        pop_id = event.replace("C=", "")
        popObj = self.popObj_by_pop_id.get(pop_id, None)
        if not popObj is None:
            other_popObj_by_pop_id = {other_pop_id: other_popObj for other_pop_id, other_popObj in self.popObj_by_pop_id.items() if not other_pop_id == pop_id}
            for i, j in itertools.combinations(range(len(popObj)), 2):
                ancestor_popObj_by_pop_id = {other_pop_id: other_popObj for other_pop_id, other_popObj in other_popObj_by_pop_id.items()}
                coalesced_alleles = ["+".join(sorted(popObj[i].split("+") + popObj[j].split("+")))]
                ancestor_popObj_by_pop_id[pop_id] = PopObj(coalesced_alleles) + PopObj([allele for k, allele in enumerate(popObj) if k not in set([i, j])])
                ancestor_stateObj = StateObj(ancestor_popObj_by_pop_id)
                self.count_by_ancestor_labels_by_event[event][ancestor_stateObj.label] += 1
    
    def _get_migration_ancestors(self, event):
        #print("# migrate...")
        pop_1_id, pop_2_id = event.replace("M=", "").split(">")
        if all([pop_id in self.popObj_by_pop_id for pop_id in [pop_1_id, pop_2_id]]) == True:
            for i, allele in enumerate(self.popObj_by_pop_id[pop_1_id]):
                ancestor_popObj_by_pop_id = {pop_id: popObj for pop_id, popObj in self.popObj_by_pop_id.items() if pop_id not in set([pop_1_id, pop_2_id])}
                source_popObj = PopObj([allele for k, allele in enumerate(self.popObj_by_pop_id[pop_1_id]) if not k == i])
                if source_popObj:
                    ancestor_popObj_by_pop_id[pop_1_id] = source_popObj
                ancestor_popObj_by_pop_id[pop_2_id] = self.popObj_by_pop_id[pop_2_id] + PopObj([allele])
                ancestor_stateObj = StateObj(ancestor_popObj_by_pop_id)
                self.count_by_ancestor_labels_by_event[event][ancestor_stateObj.label] += 1

class StateGraph(object):
    def __init__(self):
        # default is DiGraph
        self.graph = nx.MultiDiGraph(label='', meta='', mutype='', key='') 
        self.node_idx_by_label = {}
        self.label_by_node_idx = {}
        self.origin_idx_by_label = {}
        self.sample_label = None
        self.events_counter = collections.Counter()
        self.events_counter_by_node_idx = collections.defaultdict(collections.Counter)
        self.edges_by_event = collections.defaultdict(list)

    def get_origins(self):
        self.origin_idx_by_label = {self.graph.nodes[node_idx]['label']: node_idx for node_idx, out_degree in self.graph.out_degree() if out_degree == 0}

    def write_model(self, parameterObj):
        print("[+] Calculating model ...")
        #all_simple_paths_by_origin_idx = collections.defaultdict(list)
        events = sorted(self.events_counter.keys())
        header = ['origin', 'idx', 'label'] + events
        rows = []
        path_counter_by_origin = collections.Counter()
        max_length_label = max([len(label) for label in self.origin_idx_by_label.keys()])
        for origin_label, origin_idx in sorted(self.origin_idx_by_label.items()):
            idx_paths = [(path_idx, path) for path_idx, path in enumerate(nx.all_simple_paths(self.graph, source=0, target=origin_idx))]
            for path_idx, path in tqdm(idx_paths, total=len(idx_paths), desc="[+]\tPaths to %s" % origin_label.rjust(max_length_label), unit='', ncols=100, unit_scale=True):
                path_counter_by_origin[origin_label] += 1
                path_counter_by_origin['all'] += 1
                for node_idx, _ in pairwise(path):
                    row = [origin_idx]
                    row.append(node_idx)
                    row.append(self.label_by_node_idx[node_idx])
                    for event in events:
                        row.append(self.events_counter_by_node_idx[node_idx][event])
                    rows.append(row)
                last_row = [origin_idx, origin_idx, self.label_by_node_idx[node_idx]] + [0] * len(events)
                rows.append(last_row)
        lib.functions.create_csv("%s.model.tsv" % parameterObj.out_prefix, None, header, rows, sep="\t")
        print("[+] Paths \n[+]\t%s => all\n[+]\t%s" % (
            path_counter_by_origin['all'], 
            "\n[+]\t".join(["%s => %s" % (count, origin) for origin, count in sorted(path_counter_by_origin.items()) if not origin == 'all'])))

    def plot_dot_graph(self, parameterObj):
        print("[+] Plotting graph ...")
        G = nx.drawing.nx_agraph.to_agraph(self.graph) 
        G.graph_attr.update(rankdir='BT') # 'LR'
        G.node_attr['shape'] = 'record' # box
        G.node_attr['style'] = 'filled'
        G.node_attr['fillcolor'] = 'white'
        G.node_attr['color'] = 'white'
        #for label, origin_idx in self.origin_idx_by_label.items():
        #    G.get_node(origin_idx).attr['fillcolor']='white'
        #    G.get_node(origin_idx).attr['color']='black'
        G.get_node(0).attr['rank'] = 'source'
        #G.get_node(0).attr['fillcolor']='white'
        #G.get_node(0).attr['color']='black'
        #G.get_node(0).attr['shape']='record'
        colours = ['black', 'grey75', 'grey95', 'grey50', 'grey25']
        colour_by_pops = {pop: colours[idx] for idx, pop in enumerate(sorted(set([label.split("=")[0] for label in self.node_idx_by_label.keys()])))}
        print(colour_by_pops)
        for node_idx, label in self.label_by_node_idx.items():
            newlabel = ['<<table border="0" cellspacing="0">']
            for pop in label.split(";"):
                newlabel.append('<tr><td port="port%s" border="2" color="%s">%s</td></tr>' % (node_idx, colour_by_pops[pop.split("=")[0]], pop.replace("&", "+")))
            newlabel.append('</table>>')
            #G.get_node(node_idx).attr['label']="|".join(["<f%s> %s" % (idx, pop_string) for idx, pop_string in enumerate(label.split(";"))])
            G.get_node(node_idx).attr['label']="".join(newlabel)
        for event, edgelist in self.edges_by_event.items():
            if event.startswith("C="):
                colour = 'dodgerblue'
            elif event.startswith("M="):
                colour = 'orange'
            else:
                colour = 'slateblue'
            for source, sink, count in edgelist:
                G.get_edge(source, sink).attr['color'] = colour
                G.get_edge(source, sink).attr['label'] = count
                G.get_edge(source, sink).attr['fontcolor'] = colour
                G.get_edge(source, sink).attr['penwidth'] = count
        G.layout()
        #G.draw("%s.pdf" % parameterObj.out_prefix, prog='dot', args="-Gratio='compress' -Gsize='7.5,10' -Gconcentrate=true")
        out_f = "%s.pdf" % parameterObj.out_prefix
        print("[>]\tCreated: %r" % out_f)
        G.draw(out_f, prog='dot', args="-Gratio='compress' -Gsize='7.5,10'")

    def is_directed_acyclic_graph(self):
        return nx.is_directed_acyclic_graph(self.graph)

    def add_state(self, stateObj_label):
        if not stateObj_label in self.node_idx_by_label:
            node_idx = len(self.node_idx_by_label)
            self.node_idx_by_label[stateObj_label] = node_idx
            self.label_by_node_idx[node_idx] = stateObj_label.replace("'", "")
            self.graph.add_node(node_idx, label=stateObj_label.replace("'", ""), mutype='', key='')
            if node_idx == 0:
                 self.sample_label = stateObj_label

    def add_ancestry(self, stateObj_label, ancestorObj_label, event, count):
        source_idx = self.node_idx_by_label[stateObj_label]
        sink_idx = self.node_idx_by_label[ancestorObj_label]
        self.events_counter_by_node_idx[source_idx][event] += count
        self.events_counter[event] += count 
        self.graph.add_edge(source_idx, sink_idx, count=count, event=event)
        self.edges_by_event[event].append((source_idx, sink_idx, count))

def graph_generator(parameterObj):
    state_graph = StateGraph()
    queue = [parameterObj.sample_stateObj]
    while 1:
        try:
            focal_stateObj = queue.pop()
            state_graph.add_state(focal_stateObj.label)
            focal_stateObj.get_ancestors(parameterObj)
            for event in focal_stateObj.count_by_ancestor_labels_by_event.keys():
                for ancestor_label, count in focal_stateObj.count_by_ancestor_labels_by_event[event].items():
                    if not ancestor_label in state_graph.node_idx_by_label:
                        state_graph.add_state(ancestor_label)
                        queue.append(StateObj(ancestor_label))
                    state_graph.add_ancestry(focal_stateObj.label, ancestor_label, event, count)
        except IndexError:
            break
    state_graph.get_origins()
    print("[+] StateGraph\n[+]\tNodes = %s\n[+]\tEdges = %s\n[+]\tOrigins = %s [%s]" % (
        len(state_graph.graph.nodes), 
        len(state_graph.graph.edges),
        len(state_graph.origin_idx_by_label),
        ", ".join(state_graph.origin_idx_by_label.keys())))
    return state_graph

class ParameterObj(object):
    def __init__(self, args):
        self.ploidy = int(args['--ploidy'])
        self.pop_ids = sorted(args['--pop_ids'].split(','))
        self.completely_labelled = args['--complete_labelling']
        self.out_prefix = args['--out_prefix']
        self.join_events = self._parse_join_events(args['--join_string']) if args['--join_string'] else []
        self.migration_events = self._parse_migration_events(args['--migration_string']) if args['--migration_string'] else []
        self.sample_stateObj = self._get_sample_stateObj()
        self.set_of_events = self._get_set_of_events() 
        self.model_output = False if args['--nomodel'] else True 
        self.graph_output = False if args['--nograph'] else True 

    def _get_set_of_events(self):
        return set(["C=%s" % pop_id for pop_id in self.pop_ids + self.join_events] + self.join_events + self.migration_events)

    def _get_sample_stateObj(self):
        if self.completely_labelled:
            return StateObj({pop.upper(): PopObj(["%s%s" % (pop.lower(), i+1) for i in range(self.ploidy)]) for pop in self.pop_ids}) 
        return StateObj({pop.upper(): PopObj(["%s" % pop.lower() for i in range(self.ploidy)]) for pop in self.pop_ids}) 

    def _parse_join_events(self, join_string):
        char_counter = collections.Counter(join_string)
        if not all([char_counter[pop_id] <= 1 for pop_id in self.pop_ids]) == True:
            return "[X] Duplicated population IDs"
        return lib.model.joins(ast.literal_eval(re.sub(r'([a-zA-Z]+)',r'"\1"', join_string)))

    def _parse_migration_events(self, migration_string):
        events = []
        for m in migration_string.split(','):
           if '>' in m:
               a, b = m.split('>')
               events.append("M=%s>%s" % (a, b))
           elif '<' in m:
               a, b = m.split('<')
               events.append("M=%s>%s" % (b, a))
           else:
               sys.exit("[X] %r is not a valid migration parameter.")
        return events

def main(run_params):
    try:
        main_time = timer()
        args = docopt(__doc__)
        parameterObj = ParameterObj(args)
        stateGraph = graph_generator(parameterObj)
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