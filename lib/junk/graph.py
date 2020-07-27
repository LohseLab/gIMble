#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""usage: gIMble graph -o STR -s STR [-p INT] [-M] [-e -E] [-P -G] [-h|--help]

    Options:
        -h --help                       show this
        -p, --ploidy INT                Ploidy of samples [default: 2]
        -s, --samples STR               User defined samples, e.g.:
                                            - [['b'], ['a'] : 2 pops, 1 diploid sample in each
                                            - [['b'], {'a'} : 2 pops, 1 diploid sample in each, migration ({a} -> [b])
                                            - [{'b'}, ['a'] : 2 pops, 1 diploid sample in each, migration ({b} -> [a])
        -M, --migration                 Allow migration between populations 
        -e, --exodus                    Allow exodus between populations from {}->[]
        -E, --reverse_exodus            Allow exodus between populations from []->{}
        -o, --out_prefix STR            Prefix for output files
        -P, --plot_mutypes              plot mutypes in graph
        -G, --no_graph_labels           Do not plot graphlabels
"""

from docopt import docopt
from itertools import combinations, repeat
from collections import Counter, defaultdict, OrderedDict
from networkx.drawing.nx_agraph import graphviz_layout
import networkx as nx
import matplotlib as mat
mat.use("agg")
from matplotlib.lines import Line2D
import matplotlib.pyplot as plt
from timeit import default_timer as timer
from sys import stderr, exit
from ast import literal_eval
from pandas import DataFrame

'''
[Testing]
./gIMble graph -s "[{'a'}, ['b']]" -p 2 -M -o output/master && \
./gIMble graph -s "[{'a'}, ['b']]" -p 2 -e -o output/master && \
./gIMble graph -s "[{'a'}, ['b']]" -p 2 -E -o output/master && \
./gIMble graph -s "[{'a'}, ['b']]" -p 2 -M -e -o output/master && \
./gIMble graph -s "[{'a'}, ['b']]" -p 2 -M -E -o output/master


python2 src/probs.py probs -p output/master.2_pop.2_ploidy.M.paths.txt -t 4 -A 1.0 -D 1.0 -M 2.34 -m 1.2 -T 1.4 -P 30 ; \
python2 src/probs.py probs -p output/master.2_pop.2_ploidy.E.paths.txt -t 4 -A 1.0 -D 1.0 -M 2.34 -m 1.2 -T 1.4 -P 30 ; \
python2 src/probs.py probs -p output/master.2_pop.2_ploidy.R.paths.txt -t 4 -A 1.0 -D 1.0 -M 2.34 -m 1.2 -T 1.4 -P 30 ; \
python2 src/probs.py probs -p output/master.2_pop.2_ploidy.M+E.paths.txt -t 4 -A 1.0 -D 1.0 -M 2.34 -m 1.2 -T 1.4 -P 30 ; \
python2 src/probs.py probs -p output/master.2_pop.2_ploidy.M+R.paths.txt -t 4 -A 1.0 -D 1.0 -M 2.34 -m 1.2 -T 1.4 -P 30 ; \
python2 src/probs.py probs -p output/master.2_pop.2_ploidy.M.paths.txt -t 4 -A 1.0 -D 1.3 -M 2.34 -m 1.2 -T 1.4 -P 30 ; \
python2 src/probs.py probs -p output/master.2_pop.2_ploidy.E.paths.txt -t 4 -A 1.0 -D 1.3 -M 2.34 -m 1.2 -T 1.4 -P 30 ; \
python2 src/probs.py probs -p output/master.2_pop.2_ploidy.R.paths.txt -t 4 -A 1.0 -D 1.3 -M 2.34 -m 1.2 -T 1.4 -P 30 ; \
python2 src/probs.py probs -p output/master.2_pop.2_ploidy.M+E.paths.txt -t 4 -A 1.0 -D 1.3 -M 2.34 -m 1.2 -T 1.4 -P 30 ; \
python2 src/probs.py probs -p output/master.2_pop.2_ploidy.M+R.paths.txt -t 4 -A 1.0 -D 1.3 -M 2.34 -m 1.2 -T 1.4 -P 30 

[To do]                      
plot:
    - improve setting of figure dimensions based on graph complexity
    - legend with edge types
    - account for overlapping edges (using width/styles/colours)
    - adjust alpha/colour of edge labels

    

'''
################################### Constants #################################

mutype_by_lineage = {
    ('aa') : 'fixed',
    ('bb') : 'fixed',
    ('a') : 'hetA',
    ('abb') : 'hetA',
    ('aab') : 'hetB',
    ('b') : 'hetB',
    ('ab') : 'hetAB'
    }

ANCESTOR_STRING = {True: "A:", False: "D:"}
MIGRANT_BRACKETS = {True: ["{", "}"], False: ["[", "]"]}

label_text_colour = {
    'fixed' : 'purple',
    'hetA' : 'limegreen',
    'hetB' : 'lightblue',
    'hetAB' : 'yellowgreen',
    'None' : 'darkgrey'
    }

###############################################################################


################################### Functions #################################

def create_csv(out_f, header, rows, sep):
    df = DataFrame(rows, columns=header)
    print(df)
    df.to_csv(out_f, index=False, sep=sep)

def get_colours_text(label):
    items = []
    for c in label.get_text():
        if c.isalpha():
            if items[-1].isalpha():
                items[-1] = "%s%s" % (items[-1], c)
                continue
        items.append(c)
    colours = [label_text_colour[mutype_by_lineage.get(s, 'None')] for s in items]
    return (items, colours) 

def pairwise(iterable):
    it = iter(iterable)
    a = next(it, None)
    for b in it:
        yield (a, b)
        a = b

def plot_state_graph(state_graph, out_f, parameterObj):
    edge_colour_by_event = {
        'C_derived' : 'orange',
        'C_ancestor' : 'gold',
        'C' : 'gold',
        'M' : 'deeppink',
        'E' : 'lightskyblue',
        'R' : 'lightskyblue'
    }
    edge_width_by_event = {
        'C_derived' : 1,
        'C_ancestor' : 1,
        'C' : 1,
        'M' : 1,
        'E' : 1,
        'R' : 1
    }
    edge_style_by_event = {
        'C_derived' : 'solid',
        'C_ancestor' : 'solid',
        'C' : 'solid',
        'M' : 'dashed',
        'E' : 'dotted',
        'R' : 'dotted'
    }
    legend_by_event = {
        'C_derived' : 'C_derived',
        'C_ancestor' : 'C_ancestor',
        'C' : 'Coalescence',
        'M' : 'Migration',
        'E' : 'Exodus',
        'R' : 'RExodus'
    }
    for idx, (x, y) in graphviz_layout(state_graph, prog='dot',  args="-Grotation=90").items():
        state_graph.node[idx]['pos'] = (x, y)
    node_pos_by_idx = nx.get_node_attributes(state_graph, 'pos')
    edges_by_event = defaultdict(list)
    for sink, source in state_graph.edges():
        event = state_graph[sink][source]['event']
        edges_by_event[event].append((sink, source))
    rows = len(set([y for x, y in node_pos_by_idx.values()]))
    columns = len(set([x for x, y in node_pos_by_idx.values()]))
    label_font_size = 12
    state_label_y_offset = 16
    fig = plt.figure(figsize=(columns + 2, rows + 2), dpi=10, frameon=False)
    node_size = 150
    ax = fig.add_subplot(111)
    
    event_handles, event_labels = [], []
    for event, edgelist in edges_by_event.items():
        lines = nx.draw_networkx_edges(state_graph, pos=node_pos_by_idx, edge_color=edge_colour_by_event[event], width=edge_width_by_event[event], edgelist=edgelist, alpha=1.0)
        edge_labels = {(sink, source): state_graph[sink][source]['count'] for sink, source in edgelist}
        nx.draw_networkx_edge_labels(state_graph, pos=node_pos_by_idx, edgelist=edgelist, font_color='black', edge_labels=edge_labels, label_pos=0.42, font_size=label_font_size, rotate=False, bbox=dict(facecolor=edge_colour_by_event[event], edgecolor='black'))
        for line in lines:
            line.set_linestyle(edge_style_by_event[event])
            line.set_arrowstyle('-|>', head_length=0.7)
        event_handles.append(Line2D([0, 1], [0, 1], lw=edge_width_by_event[event], color=edge_colour_by_event[event], linestyle=edge_style_by_event[event]))
        event_labels.append(legend_by_event[event])
    nx.draw_networkx_nodes(state_graph, node_pos_by_idx, node_color='gainsboro', node_size=node_size)
    label_pos_by_idx = {}
    label_by_idx = {}
    state_label_by_idx = nx.get_node_attributes(state_graph, 'label')    
    label_by_idx = state_label_by_idx
    label_pos_by_idx = {idx: (pos[0], pos[1] - state_label_y_offset) for idx, pos in node_pos_by_idx.items()}
    nx.draw_networkx_labels(state_graph, label_pos_by_idx, labels=label_by_idx, horizontalalignment='center', font_size=label_font_size, bbox=dict(edgecolor='black', facecolor='lightgrey', boxstyle='round', alpha=0.25))
    ax.invert_yaxis()
    ax.legend(event_handles, event_labels)
    plt.axis('off')
    plt.tight_layout()
    fig.savefig(out_f, dpi=300, orientation='landscape', format="png")
    plt.close(fig)

def build_state_graph(parameterObj):
    start_time = timer()
    sampleStateObj = parameterObj.sampleStateObj
    lcaStateObj = parameterObj.lcaStateObj
    print("[+] Constructing graph backwards in time: Sample %r -----> LCA %r" % (sampleStateObj.label, lcaStateObj.label))
    #state_graph = nx.DiGraph(
    #                label='', \
    #                meta='', \
    #                mutype=''\
    #                )
    state_graph = nx.MultiDiGraph(
                    label='', \
                    meta='', \
                    mutype='' \
                    )
    edges = []                                 
    node_idx_by_label = {sampleStateObj.label: sampleStateObj.idx}
    stateObj_by_idx = {} 
    queue = [sampleStateObj]
    meta = {
        sampleStateObj.label : 'source',
        lcaStateObj.label : 'sink',
    }
    while 1:
        if len(queue) == 0:
            break
        current_stateObj = queue.pop(0)
        if not current_stateObj.label == lcaStateObj.label:

            ancestorStateObj_by_label, node_idx_by_label = current_stateObj.get_ancestorStateObjs(\
                                                                                    node_idx_by_label, \
                                                                                    migration=parameterObj.migration, \
                                                                                    exodus=parameterObj.exodus, \
                                                                                    reverse_exodus=parameterObj.reverse_exodus)
            for label, ancestorStateObj in ancestorStateObj_by_label.items():
                if not label == lcaStateObj.label:
                    queue.append(ancestorStateObj)
                else:
                    lcaStateObj.idx = ancestorStateObj.idx
        state_graph.add_node(current_stateObj.idx, label=current_stateObj.label, mutype=current_stateObj.mutypes, meta=meta.get(current_stateObj.label))
        stateObj_by_idx[current_stateObj.idx] = current_stateObj
        for node_idx in current_stateObj.count_by_event_by_node_idx:
            for event, count in current_stateObj.count_by_event_by_node_idx[node_idx].items():
                edges.append((current_stateObj.idx, node_idx, {'count': count, 'event': event}))
    state_graph.add_edges_from(edges)
    print("[+] Generated graph with %s states, %s events in %s seconds." % (state_graph.number_of_nodes(), state_graph.number_of_edges(), timer() - start_time))
    return (state_graph, stateObj_by_idx)

def get_state_paths(state_graph, parameterObj):
    start_time = timer()
    # Get total outgoing edges by event_type for each node 
    total_by_event_by_node_idx = defaultdict(Counter)
    for source, sink, data in state_graph.edges.data():
        total_by_event_by_node_idx[source][data['event']] += data['count']
    mutype_by_node_idx = {node_idx: mutype for node_idx, mutype in nx.get_node_attributes(state_graph, 'mutype').items()}
    header = ['path_idx', 'node', 'event', 'event_count', 'C_ancestor', 'C_derived', 'Ms', 'Es', 'hetA', 'fixed', 'hetB', 'hetAB']
    paths = []
    path_count = 0
    # this only works for one (!) M+E event per path (this has to be generalised  if we do >2 pops...)
    path_counter = Counter()
    for path_idx, path in enumerate(nx.all_simple_paths(state_graph, source=parameterObj.sampleStateObj.idx, target=parameterObj.lcaStateObj.idx)):
        path_id = tuple(path)
        step_idx = 0
        path_count += 1
        for source, sink in pairwise(path):
            data = state_graph[source][sink].get(path_counter[path_id], state_graph[source][sink][0])
            line = []
            line.append(path_idx)
            line.append(source)
            if data['event'] == 'R':
                line.append('E')
            else:
                line.append(data['event'])
            line.append(data['count'])
            line.append(total_by_event_by_node_idx[source]['C_ancestor'])
            line.append(total_by_event_by_node_idx[source]['C_derived'])
            line.append(total_by_event_by_node_idx[source]['M'])
            if parameterObj.reverse_exodus:
                line.append(total_by_event_by_node_idx[source]['R'])
            else:
                line.append(total_by_event_by_node_idx[source]['E'])
            line.append(mutype_by_node_idx[source]['hetA'])
            line.append(mutype_by_node_idx[source]['fixed'])
            line.append(mutype_by_node_idx[source]['hetB'])
            line.append(mutype_by_node_idx[source]['hetAB'])
            paths.append(line)
            step_idx += 1
        path_counter[path_id] += 1
        paths.append([path_idx, 'LCA', 'LCA', 0, 0, 0, 0, 0, 0, 0, 0, 0])
    print("[+] Found %s paths in %s seconds." % (sum(path_counter.values()), timer() - start_time))
    return (header, paths)
        


    #for path_idx, path in enumerate(nx.all_simple_paths(state_graph, source=node_idx_by_meta['source'], target=node_idx_by_meta['sink'])):
    #unique_single_paths = set([tuple(path) for path in nx.all_simple_paths(state_graph, source=parameterObj.sampleStateObj.idx, target=parameterObj.lcaStateObj.idx)])

    #combined_single_paths = []
    #for path in unique_single_paths:
    #    multipaths = [True if len(state_graph[source][sink]) > 1 else False for source, sink in pairwise(path)]
    #    print(multipaths)
    #    if any(multipaths):
    #        multi_path_idxs = [idx for idx, multipath in enumerate(multipaths) if idx == True]

    #        print(multi_path_idxs)
    #    print(source, sink, state_graph[source][sink])
            

    exit()
    for path_idx, path in enumerate(nx.all_simple_paths(state_graph, source=parameterObj.sampleStateObj.idx, target=parameterObj.lcaStateObj.idx)):

        path_count += 1
        step_idx = 0
        print(path)
        for source, sink in pairwise(path):
            print(source,sink)
            line = []
            line.append(path_idx)
            line.append(source)
            print(state_graph[source])
            print(data)
            if data['event'] == 'R':
                line.append('E')
            else:
                line.append(data['event'])
            line.append(data['count'])
            line.append(total_by_event_by_node_idx[source]['C_ancestor'])
            line.append(total_by_event_by_node_idx[source]['C_derived'])
            line.append(total_by_event_by_node_idx[source]['M'])
            if parameterObj.reverse_exodus:
                line.append(total_by_event_by_node_idx[source]['R'])
            else:
                line.append(total_by_event_by_node_idx[source]['E'])
            line.append(mutype_by_node_idx[source]['hetA'])
            line.append(mutype_by_node_idx[source]['fixed'])
            line.append(mutype_by_node_idx[source]['hetB'])
            line.append(mutype_by_node_idx[source]['hetAB'])
            paths.append(line)
            step_idx += 1
        paths.append([path_idx, 'LCA', 'LCA', 0, 0, 0, 0, 0, 0, 0, 0, 0])
    print("[+] Found %s paths in %s seconds." % (path_count, timer() - start_time))
    return (header, paths)

###############################################################################

################################### Classes ##################################


class StateObj(object):
    def __init__(self, popObjs, idx=None):
        self.idx = idx
        self.popObjs = popObjs   
        self.label = " ".join([popObj.label for popObj in popObjs])
        self.population_idx_by_type = self.get_population_idx_by_type()
        self.lineages = self._get_lineages()
        self.mutypes = Counter([mutype_by_lineage.get(lineage) for lineage in self.lineages])
        self.count_by_event_by_node_idx = defaultdict(lambda: Counter()) # outgoing edge_weight by idx (target) and event

    def __eq__(self, other):
        if isinstance(other, StateObj):
            return self.label == other.label
    
    def __str__(self):
        return self.label

    def get_population_idx_by_type(self):
        population_idx_by_type = {}
        for idx, popObj in enumerate(self.popObjs):
            if popObj.migrant == True:
                population_idx_by_type['migrant'] = idx
            else:
                population_idx_by_type['resident'] = idx
            if popObj.ancestor == True:
                population_idx_by_type['ancestor'] = idx
            else:
                population_idx_by_type['derived'] = idx
        return population_idx_by_type
    
    def _get_lineages(self):
        lineages = []
        for popObj in self.popObjs:
            for lineage in popObj.lineages:
                lineages.append(lineage)
        return lineages

    def get_popObj(self, agent):
        # only implemented for 1 population per agent...
        return self.popObjs[self.population_idx_by_type[agent]]

    def get_coalescenceStateObjs(self, node_idx_by_label):
        ancestorStateObj_by_label = {}
        for i, popObj in enumerate(self.popObjs):
            if popObj.lineage_count > 1:
                coalescedPopObjs = popObj.coalesce()
                if popObj.ancestor:
                    event = 'C_ancestor'
                    otherPopObj = self.get_popObj('derived')
                    coalesced_pop_idx = self.population_idx_by_type['ancestor']
                    other_pop_idx = self.population_idx_by_type['derived']
                else:
                    event = 'C_derived'
                    otherPopObj = self.get_popObj('ancestor')
                    coalesced_pop_idx = self.population_idx_by_type['derived']
                    other_pop_idx = self.population_idx_by_type['ancestor']
                for coalescedPopObj in coalescedPopObjs:
                    ancestorPopObjs = [None, None]
                    ancestorPopObjs[coalesced_pop_idx] = coalescedPopObj
                    ancestorPopObjs[other_pop_idx] = otherPopObj
                    ancestor_state_label = " ".join([ancestorPopObj.label for ancestorPopObj in ancestorPopObjs])
                    if not ancestor_state_label in node_idx_by_label:
                        node_idx = max(node_idx_by_label.values()) + 1
                        ancestorStateObj_by_label[ancestor_state_label] = StateObj(ancestorPopObjs, idx=node_idx)
                        node_idx_by_label[ancestor_state_label] = node_idx
                    self.count_by_event_by_node_idx[node_idx_by_label[ancestor_state_label]][event] += 1
        return (ancestorStateObj_by_label, node_idx_by_label)

    def get_migrationStateObjs(self, node_idx_by_label):
        ancestorStateObj_by_label = {}
        migrantPopObj = self.get_popObj('migrant')
        if migrantPopObj.lineage_count > 0:
            event = 'M'
            migrant_idx = self.population_idx_by_type['migrant']
            residentPopObj = self.get_popObj('resident')
            resident_idx = self.population_idx_by_type['resident']
            for i, lineage in enumerate(migrantPopObj.lineages):
                ancestorPopObjs = [None, None]
                ancestorPopObjs[migrant_idx] = PopObj(sorted([migrant for j, migrant in enumerate(migrantPopObj.lineages) if not i == j]), migrant=True, ancestor=migrantPopObj.ancestor)
                ancestorPopObjs[resident_idx] = PopObj(sorted(residentPopObj.lineages + [lineage]), migrant=False, ancestor=residentPopObj.ancestor)
                ancestor_state_label = " ".join([ancestorPopObj.label for ancestorPopObj in ancestorPopObjs])
                if not ancestor_state_label in node_idx_by_label:
                    node_idx = max(node_idx_by_label.values()) + 1
                    ancestorStateObj = StateObj(ancestorPopObjs, idx=node_idx)
                    ancestorStateObj_by_label[ancestorStateObj.label] = ancestorStateObj
                    node_idx_by_label[ancestorStateObj.label] = node_idx
                self.count_by_event_by_node_idx[node_idx_by_label[ancestor_state_label]][event] += 1
        return (ancestorStateObj_by_label, node_idx_by_label)

    def get_exodusStateObjs(self, node_idx_by_label, reverse=False):
        ancestorStateObj_by_label = {}
        event, exodus_popObj, exodus_idx, empty_popObj, empty_idx = None, None, None, None, None
        ancestorPopObj = self.get_popObj('ancestor')
        derivedPopObj = self.get_popObj('derived')
        if derivedPopObj.lineage_count > 0:
            if reverse == False:
                event = 'E'
            else:
                event = 'R'
            exodus_idx = self.population_idx_by_type['ancestor']
            exodus_popObj = PopObj(sorted(self.lineages), migrant=ancestorPopObj.migrant, ancestor=ancestorPopObj.ancestor)
            empty_idx = self.population_idx_by_type['derived']
            empty_popObj = PopObj([], migrant=derivedPopObj.migrant, ancestor=derivedPopObj.ancestor)
        if not event is None:
            ancestorStateObj_by_label = {}
            ancestorPopObjs = [None, None]
            ancestorPopObjs[exodus_idx] = exodus_popObj
            ancestorPopObjs[empty_idx] = empty_popObj
            ancestor_state_label = " ".join([ancestorPopObj.label for ancestorPopObj in ancestorPopObjs])
            if not ancestor_state_label in node_idx_by_label:
                node_idx = max(node_idx_by_label.values()) + 1
                ancestorStateObj = StateObj(ancestorPopObjs, idx=node_idx)
                ancestorStateObj_by_label[ancestorStateObj.label] = ancestorStateObj
                node_idx_by_label[ancestorStateObj.label] = node_idx
            self.count_by_event_by_node_idx[node_idx_by_label[ancestor_state_label]][event] += 1
        return (ancestorStateObj_by_label, node_idx_by_label)

    def get_ancestorStateObjs(self, node_idx_by_label, migration=False, exodus=False, reverse_exodus=False):
        ancestorStateObj_by_label, node_idx_by_label = self.get_coalescenceStateObjs(node_idx_by_label)                
        if migration:
            if not len(self.get_popObj('derived').lineages) == 0:
                migrationStateObjs_by_label, node_idx_by_label = self.get_migrationStateObjs(node_idx_by_label)
                ancestorStateObj_by_label = {**ancestorStateObj_by_label, **migrationStateObjs_by_label}
        if exodus:
            exodusStateObjs_by_label, node_idx_by_label = self.get_exodusStateObjs(node_idx_by_label, reverse=False)
            ancestorStateObj_by_label = {**ancestorStateObj_by_label, **exodusStateObjs_by_label}
        if reverse_exodus:
            rexodusStateObjs_by_label, node_idx_by_label = self.get_exodusStateObjs(node_idx_by_label, reverse=True)
            ancestorStateObj_by_label = {**ancestorStateObj_by_label, **rexodusStateObjs_by_label}
        return (ancestorStateObj_by_label, node_idx_by_label)

class PopObj(object):
    def __init__(self, population, ploidy=None, ancestor=bool, migrant=bool):
        # population must be list
        self.lineages = population if ploidy is None \
                else [lineage for lineages in population for lineage in repeat(lineages, ploidy)]
        self.lineage_count = len(self.lineages)
        self.ancestor = ancestor
        self.migrant = migrant
        if self.lineages:
            self.label = "".join([\
                ANCESTOR_STRING[self.ancestor], \
                MIGRANT_BRACKETS[migrant][0], \
                ",".join(self.lineages), \
                MIGRANT_BRACKETS[migrant][1]])
        else:
            self.label = "".join([\
                ANCESTOR_STRING[self.ancestor], \
                MIGRANT_BRACKETS[migrant][0], 
                MIGRANT_BRACKETS[migrant][1]])
    
    def __hash__(self):
        return hash(self.label)

    def __eq__(self, other):
        if isinstance(other, StateObj):
            return self.label == other.label

    def __str__(self):
        return self.label

    def coalesce(self):
        popObjs = []
        for i, j in combinations(range(len(self.lineages)), 2):
            lineages = [s for idx, s in enumerate(self.lineages) if not idx == i and not idx == j]
            lineages.append("".join(sorted(self.lineages[i] + self.lineages[j])))
            popObjs.append(PopObj(sorted(lineages), ancestor=self.ancestor, migrant=self.migrant))
        return popObjs

class ParameterObj(object):
    def __init__(self, args):
        self.ploidy = int(args['--ploidy'])
        self.migration = args['--migration']
        self.plot_mutypes = args['--plot_mutypes']
        self.no_graph_labels = args['--no_graph_labels']
        self.exodus = args['--exodus']
        self.reverse_exodus = args['--reverse_exodus']
        self.out_prefix = args['--out_prefix']
        self.sample_input = args['--samples'] 
        self.migration_possible = any([self.migration, self.exodus, self.reverse_exodus])
        self.population_count = 0
        self.model_label = self.get_model_label()

    def get_model_label(self):
        models = OrderedDict({
            'M': self.migration, 
            'E': self.exodus, 
            'R': self.reverse_exodus
        })
        return "+".join([model_id for model_id, status in models.items() if status])

    def setup(self):
        print("=======================================")
        print("[+] Sample input: ", self.sample_input)
        print("=======================================")
        try:
            sample_string = literal_eval(self.sample_input)
        except SyntaxError:
            exit("[X] Invalid sample (--sample) string: %s" % self.sample_input)
        # Fail if not sets/lists
        if not isinstance(sample_string, list) or not all([(type(population) in set([list, set])) for population in sample_string]):
            exit("[X] Samples (--samples) must list of lists/sets '[['b'], {'a'}]' !")
        # Fail if migration_possible and only one population 
        if self.migration_possible:
            if not any([isinstance(population, set) for population in sample_string]):
                exit("[X] Must specify migrant population using '{}' when using '--migration' and/or '--exodus'!")
        # Get migrants (can only be set or list)
        migrants = [(True if isinstance(population, set) else False) for population in sample_string]
        # Setting ancestor in samplePopObjs
        migrants_counter = Counter(migrants)
        samplePopObjs = []
        ancestor_idx = None
        if self.reverse_exodus: 
            # Fail if more than one ancestor population
            if not migrants_counter[True] == 1:
                exit("[X] Coalescing into ancestor state based on this model is only possible for 1 migrant population!")
            # Set migrant population to ancestor
            ancestor_idx = migrants.index(True)
        else:
            # Fail if no migration and multiple populations
            if not self.migration_possible:
                if not len(samplePopObjs) == 1:
                    exit("[X] Coalescing into ancestor state from %s populations is not possible without '--migration' and/or '--exodus'!" % len(samplePopObjs))
            # Fail if more than one ancestor population
            if not migrants_counter[False] == 1:
                exit("[X] Coalescing into ancestor state based on this model is only possible for 1 resident population!")
            # Set resident population to ancestor
            ancestor_idx = migrants.index(False)
        for idx, (population, migrant) in enumerate(zip(sample_string, migrants)):
            if idx == ancestor_idx:
                ancestor = True
            else:
                ancestor = False
            samplePopObjs.append(PopObj(list(population), migrant=migrant, ancestor=ancestor, ploidy=self.ploidy))
        self.population_count = len(samplePopObjs)
        # Create lcaPopObjs
        lcaPopObjs = []
        for samplePopObj in samplePopObjs:
            lca_population = []
            if samplePopObj.ancestor:
                # get LCA from popObj.lineage strings
                lca_population = ["".join(sorted(["".join(samplePopObj.lineages) for samplePopObj in samplePopObjs]))]
            lcaPopObj = PopObj(lca_population, ancestor=samplePopObj.ancestor, migrant=samplePopObj.migrant)
            lcaPopObjs.append(lcaPopObj)
        if not self.population_count == 2:
             exit("[X] Beyond this point only implemented for 2 populations!")
        # Create sampleStateObj
        self.sampleStateObj = StateObj(samplePopObjs, idx=0)
        print("[=] Sample state =", self.sampleStateObj.label)
        # Create lca_stateObj 
        self.lcaStateObj = StateObj(lcaPopObjs, idx=None)
        print("[=] LCA state    =", self.lcaStateObj)
        # Fail if population count not 2 ... 
        print("---------------------------------------")
        
        print("[=] Model        = %s" % self.model_label)
        print("---------------------------------------")
        if self.migration:
            print("[=] Migration : %s -> %s" % (self.sampleStateObj.get_popObj('migrant'), self.sampleStateObj.get_popObj('resident')))
        if self.exodus:
            print("[=] Exodus    : %s -> %s" % (self.sampleStateObj.get_popObj('migrant'), self.sampleStateObj.get_popObj('resident')))
        if self.reverse_exodus:
            print("[=] Exodus    : %s <- %s" % (self.sampleStateObj.get_popObj('migrant'), self.sampleStateObj.get_popObj('resident')))
        print("---------------------------------------")
        
    def get_basename(self):
        return "%s.%s_pop.%s_ploidy.%s" % (self.out_prefix, self.population_count, self.ploidy, self.model_label)

    def write_paths(self, header, paths):
        start_time = timer()
        out_f = "%s.paths.txt" % (self.get_basename())
        create_csv(out_f, header, paths, sep="\t")
        print("[+] Written paths into file %s in %s seconds." % (out_f, timer() - start_time))

    def write_nodes(self, header, nodes):
        start_time = timer()
        out_f = "%s.nodes.txt" % (self.get_basename())
        create_csv(out_f, header, nodes, sep="\t")
        print("[+] Written paths into file %s in %s seconds." % (out_f, timer() - start_time))

    def write_graph(self, state_graph):
        start_time = timer()
        out_f = "%s.gml" % (self.get_basename())
        # states have to be stringified ...
        for idx, node in state_graph.nodes(data=True):
            # print(idx, node)
            node['label'] = str(node.get('label', None))
        nx.write_gml(state_graph, out_f)
        print("[+] Written graph into file %s in %s seconds." % (out_f, timer() - start_time))
        return out_f

    def read_graph(self, infile):
        start_time = timer()
        state_graph = nx.read_gml(infile)
        for idx, node in state_graph.nodes(data=True):
            node['label'] = literal_eval(node['label'])
            if 'pos' in node:
                node['pos'] = tuple(node['pos'])
            # print(idx, node)
        print("[+] Read graph from file %s in %s seconds." % (infile, timer() - start_time))
        return state_graph

    def plot_graph(self, state_graph, parameterObj):
        start_time = timer()
        out_f = "%s.graph.png" % self.get_basename()
        plot_state_graph(state_graph=state_graph, out_f=out_f, parameterObj=parameterObj)
        print("[+] Plotted graph into file %s in %s seconds." % (out_f, timer() - start_time))
        return out_f

###############################################################################

################################### Main ######################################

def main():
    try:
        main_time = timer()
        args = docopt(__doc__)
        parameterObj = ParameterObj(args)
        parameterObj.setup()
        state_graph, stateObj_by_idx = build_state_graph(parameterObj)
        header, paths = get_state_paths(state_graph, parameterObj)
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