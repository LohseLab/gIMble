#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""usage: gIMble graph -o STR -s STR [-p INT] [-M] [-E -R] [-P -G] [-h|--help]

    Options:
        -h --help                       show this
        -p, --ploidy INT                Ploidy of samples [default: 2]
        -s, --samples STR               User defined samples, e.g.:
                                            - [['b'], ['a'] : 2 pops, 1 diploid sample in each
                                            - [['b'], {'a'} : 2 pops, 1 diploid sample in each, migration ({a} -> [b])
                                            - [{'b'}, ['a'] : 2 pops, 1 diploid sample in each, migration ({b} -> [a])
        -M, --migration                 Allow migration between populations 
        -E, --exodus                    Allow exodus between populations from {}->[]
        -R, --reverse_exodus            Allow exodus between populations from []->{}
        -o, --out_prefix STR            Prefix for output files
        -P, --plot_mutypes              plot mutypes in graph
        -G, --no_graph_labels           Do not plot graphlabels
"""

from docopt import docopt
from itertools import combinations
from collections import Counter, defaultdict, OrderedDict
from networkx.drawing.nx_agraph import graphviz_layout
import networkx as nx
import matplotlib as mat
mat.use("agg")
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.offsetbox import AnchoredOffsetbox, HPacker, TextArea
import re
from matplotlib import transforms
import matplotlib.pyplot as plt
from timeit import default_timer as timer
from sys import stderr, exit
from ast import literal_eval
from pandas import DataFrame

'''
[Testing]
./gIMble graph -s "[{'a'}, ['b']]" -p 2 -M -o output/master && \
./gIMble graph -s "[{'a'}, ['b']]" -p 2 -E -o output/master && \
./gIMble graph -s "[{'a'}, ['b']]" -p 2 -R -o output/master && \
./gIMble graph -s "[{'a'}, ['b']]" -p 2 -M -E -o output/master && \
./gIMble graph -s "[{'a'}, ['b']]" -p 2 -M -R -o output/master

python2 src/probs.py probs -p output/master.2_samples.2_pop.2_ploidy.E.paths.txt -t 1 -A 1.0 -D 1.3 -M 2.34 -T 1.4 && \
python2 src/probs.py probs -p output/master.2_samples.2_pop.2_ploidy.M+E.paths.txt -t 1 -A 1.0 -D 1.3 -M 2.34 -T 1.4 && \
python2 src/probs.py probs -p output/master.2_samples.2_pop.2_ploidy.M+R.paths.txt -t 1 -A 1.0 -D 1.3 -M 2.34 -T 1.4 && \
python2 src/probs.py probs -p output/master.2_samples.2_pop.2_ploidy.M.paths.txt -t 1 -A 1.0 -D 1.3 -M 2.34 -T 1.4 && \
python2 src/probs.py probs -p output/master.2_samples.2_pop.2_ploidy.R.paths.txt -t 1 -A 1.0 -D 1.3 -M 2.34 -T 1.4 && \
python2 src/probs.py probs -p output/master.2_samples.2_pop.2_ploidy.E.paths.txt -t 1 -A 1.0 -D 1.0 -M 2.34 -T 1.4 && \
python2 src/probs.py probs -p output/master.2_samples.2_pop.2_ploidy.M+E.paths.txt -t 1 -A 1.0 -D 1.0 -M 2.34 -T 1.4 && \
python2 src/probs.py probs -p output/master.2_samples.2_pop.2_ploidy.M+R.paths.txt -t 1 -A 1.0 -D 1.0 -M 2.34 -T 1.4 && \
python2 src/probs.py probs -p output/master.2_samples.2_pop.2_ploidy.M.paths.txt -t 1 -A 1.0 -D 1.0 -M 2.34 -T 1.4 && \
python2 src/probs.py probs -p output/master.2_samples.2_pop.2_ploidy.R.paths.txt -t 1 -A 1.0 -D 1.0 -M 2.34 -T 1.4

python2 src/probs.py probs -p output/master.2_samples.2_pop.2_ploidy.M+R.paths.txt -t 1 -A 1.0 -D 1.3 -M 2.34 -T 1.4 && \
python2 src/probs.py probs -p output/master.2_samples.2_pop.2_ploidy.M+R.paths.txt -t 1 -A 1.0 -D 1.0 -M 2.34 -T 1.4 

[Problems] 
- stable sage needs python2 !?!
    - hacky compiling with python3 : https://wiki.sagemath.org/Python3-compatible%20code
        or https://wiki.sagemath.org/Conda
    - official version soon: 
        - https://trac.sagemath.org/ticket/26212
        - https://trac.sagemath.org/ticket/15530
- path computation can be improved using networkit/graph-tools
- but bottleneck is NOT path computation but creation of nodes in graph
    - can be parallelised, but some sort of LOCK needs to be implemented
    - ancestors are independent, so each node can be "ancestorised" independently and put things to queue
                           
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

def plot_state_graph(state_graph, out_f, mutypes_flag=False, no_graph_labels=False):
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
        'C_derived' : '{Coalescence}',
        'C_ancestor' : '[Coalescence]',
        'C' : 'Coalescence',
        'M' : 'Migration',
        'E' : 'Exodus',
        'R' : 'Exodus'
    }
    for idx, (x, y) in graphviz_layout(state_graph, prog='dot',  args="-Grotation=90").items():
        state_graph.node[idx]['pos'] = (x, y)
    node_pos_by_idx = nx.get_node_attributes(state_graph, 'pos')
    
    edges_by_event = defaultdict(list)
    for sink, source in state_graph.edges():
        event = state_graph[sink][source]['event']
        edges_by_event[event].append((sink, source))
    number_of_nodes = len(state_graph)
    rows = len(set([y for x, y in node_pos_by_idx.values()]))
    columns = len(set([x for x, y in node_pos_by_idx.values()]))
    label_font_size = legend_size = 12
    fig, state_label_y_offset = None, None
    
    # graph_labels = False
    if number_of_nodes >= 30:
        fig = plt.figure(figsize=(min(columns + 1, 8), min(rows + 1, 7)), dpi=10, frameon=False)
        state_label_y_offset = 24
        no_graph_labels = True
    else:
        fig = plt.figure(figsize=(columns + 2, rows + 2), dpi=10, frameon=False)
        state_label_y_offset = 14 + (rows - 4)
    node_size = 150
    ax = fig.add_subplot(111)
    
    event_handles, event_labels = [], []
    for event, edgelist in edges_by_event.items():
        lines = nx.draw_networkx_edges(state_graph, pos=node_pos_by_idx, edge_color=edge_colour_by_event[event], width=edge_width_by_event[event], edgelist=edgelist, alpha=1.0)
        if not no_graph_labels:
            edge_labels = {(sink, source): state_graph[sink][source]['count'] for sink, source in edgelist}
            nx.draw_networkx_edge_labels(state_graph, pos=node_pos_by_idx, edgelist=edgelist, font_color='black', edge_labels=edge_labels, label_pos=0.42, font_size=label_font_size, rotate=False, bbox=dict(facecolor=edge_colour_by_event[event], edgecolor=edge_colour_by_event[event]))
                #line.set_connectionstyle("bar", angle=90, fraction=0.2)  
        for line in lines:
            line.set_linestyle(edge_style_by_event[event])
            line.set_arrowstyle('-|>', head_length=0.7)
        event_handles.append(Line2D([0, 1], [0, 1], lw=edge_width_by_event[event], color=edge_colour_by_event[event], linestyle=edge_style_by_event[event]))
        event_labels.append(legend_by_event[event])
    nx.draw_networkx_nodes(state_graph, node_pos_by_idx, node_color='gainsboro', node_size=node_size)
    
    # node labels
    label_pos_by_idx = {}
    label_by_idx = {}
    state_label_by_idx = nx.get_node_attributes(state_graph, 'state')
    if no_graph_labels:
        node_idx_by_meta = {meta: node_idx for node_idx, meta in nx.get_node_attributes(state_graph, 'meta').items() if meta}
        source_idx = node_idx_by_meta['source']
        sink_idx = node_idx_by_meta['sink']
        label_by_idx[source_idx] = state_label_by_idx[source_idx]
        label_pos_by_idx[source_idx] = (state_graph.nodes[source_idx]['pos'][0], state_graph.nodes[source_idx]['pos'][1] - state_label_y_offset)
        label_by_idx[sink_idx] = state_label_by_idx[sink_idx]
        label_pos_by_idx[sink_idx] = (state_graph.nodes[sink_idx]['pos'][0], state_graph.nodes[sink_idx]['pos'][1] - state_label_y_offset)
    else:
        label_by_idx = state_label_by_idx
        label_pos_by_idx = {idx: (pos[0], pos[1] - state_label_y_offset) for idx, pos in node_pos_by_idx.items()}
        nx.draw_networkx_labels(state_graph, label_pos_by_idx, labels=label_by_idx, horizontalalignment='center', font_size=label_font_size) #, bbox=dict(edgecolor='white', facecolor=state_label_fcs[mutype], boxstyle='round', alpha=0.25))
        # state_label_handles = nx.draw_networkx_labels(state_graph, label_pos_by_idx, labels=label_by_idx, horizontalalignment='center', font_size=label_font_size) #, bbox=dict(edgecolor='white', facecolor=state_label_fcs[mutype], boxstyle='round', alpha=0.25))
        #for label in state_label_handles.values():
        #    pos = label.get_position()
        #    items, colours = get_colours_text(label)
        #    children = []
        #    for item, colour in zip(items, colours):
        #        children.append(TextArea(items, textprops=dict(color=colour, size=label_font_size)))
        #    txt = HPacker(children=children, align="baseline", pad=0, sep=0)
        #    txt.set_offset(ax.transData.transform_point(pos))
        #    ax.add_artist(txt)    
    ax.invert_yaxis()
    plt.axis('off')
    plt.tight_layout()
    fig.savefig(out_f, dpi=300, orientation='landscape', format="png")
    plt.close(fig)

def coalesce(stateObj, migrant=False, resident=False):
    ancestor_states, events, counts = [], [], []
    if migrant:
        event = 'C_derived'
        state = stateObj.get_migrant_list()
        other = stateObj.get_resident_list()
    elif resident:
        event = 'C_ancestor'
        state = stateObj.get_resident_list()
        other = stateObj.get_migrant_list()
    else:
        # single population
        event = 'C_ancestor'
        state = stateObj.state[0]
        other = []
    _ancestor_states = []
    for i, j in combinations(range(len(state)), 2):
        coalescence_state = [s for idx, s in enumerate(state) if not idx == i and not idx == j]
        coalescence_state.append("".join(sorted(state[i] + state[j])))
        ancestor_state = [None, None]
        if migrant:
            ancestor_state[stateObj.mig_idx] = coalescence_state
            ancestor_state[stateObj.res_idx] = other
        elif resident:
            ancestor_state[stateObj.mig_idx] = other
            ancestor_state[stateObj.res_idx] = coalescence_state
        else:
            ancestor_state = [coalescence_state]
        _ancestor_states.append(state_from_lol(ancestor_state))
    for ancestor_state, count in Counter(_ancestor_states).items():
        ancestor_states.append(ancestor_state)
        events.append(event)
        counts.append(count)
    return ancestor_states, events, counts

def get_coalescence_ancestors(stateObj, migration=False):
    if migration == True:
        migrant_ancestor_states, migrant_events, migrant_counts = coalesce(stateObj, migrant=True)
        resident_ancestor_states, resident_events, resident_counts = coalesce(stateObj, resident=True)
        ancestor_states = migrant_ancestor_states + resident_ancestor_states
        events = migrant_events + resident_events
        counts = migrant_counts + resident_counts
        return ancestor_states, events, counts    
    else:
        ancestor_states, events, counts = coalesce(stateObj)
        return ancestor_states, events, counts

def get_ancestors(stateObj, parameterObj):
    results = []
    if parameterObj.migration_possible:
        if not stateObj.resident_state == ():
            if parameterObj.reverse_exodus:
                results.append(get_reverse_exodus_ancestors(stateObj))
            if parameterObj.migration:
                results.append(get_migration_ancestors(stateObj))
        if not stateObj.migrant_state == ():
            if parameterObj.exodus:
                results.append(get_exodus_ancestors(stateObj))
        results.append(get_coalescence_ancestors(stateObj, migration=True))
    else:
        results.append(get_coalescence_ancestors(stateObj, migration=False))
    ancestor_states, events, counts = [], [], []
    for result in results:
        _ancestor_states, _events, _counts = result
        ancestor_states.extend(_ancestor_states)
        events.extend(_events)
        counts.extend(_counts)
    return ancestor_states, events, counts

def get_migration_ancestors(stateObj):
    ancestor_states, events, counts, event = [], [], [], 'M'
    _ancestor_states = []
    migrants = stateObj.get_migrant_list()
    for i in range(len(stateObj.migrant_state)):
        residents = stateObj.get_resident_list()
        residents.append(migrants[i])
        ancestor_state = [None, None]
        ancestor_state[stateObj.mig_idx] = sorted([migrant for j, migrant in enumerate(migrants) if not i == j])
        ancestor_state[stateObj.res_idx] = sorted(residents)
        _ancestor_states.append(state_from_lol(ancestor_state))
    for ancestor_state, count in Counter(_ancestor_states).items():
        ancestor_states.append(ancestor_state)
        events.append(event)
        counts.append(count)
    return ancestor_states, events, counts

def get_reverse_exodus_ancestors(stateObj):
    ancestor_state = [None, None]
    ancestor_state[stateObj.mig_idx] = sorted(stateObj.get_migrant_list() + stateObj.get_resident_list())
    ancestor_state[stateObj.res_idx] = []
    ancestor_state = state_from_lol(ancestor_state)
    return [ancestor_state], ['E'], [1] 

def get_exodus_ancestors(stateObj):
    ancestor_state = [None, None]
    ancestor_state[stateObj.mig_idx] = []
    ancestor_state[stateObj.res_idx] = sorted(stateObj.get_migrant_list() + stateObj.get_resident_list())
    ancestor_state = state_from_lol(ancestor_state)
    return [ancestor_state], ['E'], [1] 

def build_state_graph(sampleStateObj, lcaStateObj, parameterObj):
    start_time = timer()
    print("[+] Constructing graph backwards in time: Sample %r -----> LCA %r" % (sampleStateObj.state, lcaStateObj.state))
    state_graph = nx.DiGraph(state='', meta='', mutype='')
    edges = []                                 
    idx_by_state = {} # {StateObj.state: StateObj.idx}
    stateObj_by_idx = {} # {StateObj.idx: StateObj}
    # queue => takes StateObjs and returns (state, events, counts)
    queue = [sampleStateObj]
    meta = {
        sampleStateObj.state : 'source',
        lcaStateObj.state : 'sink',
    }
    while 1:
        if len(queue) == 0:
            break
        current_stateObj = queue.pop(0)
        # check whether seen before
        current_state_idx = idx_by_state.get(current_stateObj.idx, None)
        if current_state_idx is None:
            # build stateObj
            stateObj_by_idx[current_state_idx] = current_stateObj
            mutypes = Counter([mutype_by_lineage.get(lineage, 'None') for lineage in current_stateObj.lineages])
            state_graph.add_node(current_stateObj.idx, state=str(current_stateObj), mutype=mutypes, meta=meta.get(current_stateObj.state, ''))
            #print("[+]", current_stateObj)
        ancestor_states, events, counts = get_ancestors(current_stateObj, parameterObj) # (state, events, counts, child_states)
        for ancestor_state, event, count in zip(ancestor_states, events, counts):
            #print("\t[A]", ancestor_state, event, count)
            # node
            if not ancestor_state in idx_by_state:
                ancestor_state_idx = len(state_graph)
                ancestor_stateObj = StateObj(idx=ancestor_state_idx, state=ancestor_state, res_idx=parameterObj.resident_idx, mig_idx=parameterObj.migrant_idx)
                #print("\t[N]", ancestor_state, event, count, " => ", ancestor_stateObj)
                mutypes = Counter([mutype_by_lineage.get(lineage, 'None') for lineage in current_stateObj.lineages])
                state_graph.add_node(ancestor_state_idx, state=str(ancestor_stateObj), mutype=mutypes, meta=meta.get(ancestor_state, ''))
                #print("\t\t[+]", ancestor_state_idx, ancestor_state)
                idx_by_state[ancestor_state] = ancestor_state_idx
                stateObj_by_idx[ancestor_state_idx] = ancestor_stateObj
                if ancestor_state == lcaStateObj.state:
                    lcaStateObj.idx = ancestor_state_idx
                else:
                    #print("\t\t[<]", ancestor_stateObj)
                    queue.append(ancestor_stateObj)
            else:
                ancestor_state_idx = idx_by_state[ancestor_state]
                ancestor_stateObj = StateObj(idx=ancestor_state_idx, state=ancestor_state, res_idx=parameterObj.resident_idx, mig_idx=parameterObj.migrant_idx)
                #print("\t[E]", ancestor_state, event, count, " => ", ancestor_stateObj)
            # edges
            #print("\t[Edge]", (current_stateObj.idx, ancestor_stateObj.idx, {'count': count, 'event': event}))
            edges.append((current_stateObj.idx, ancestor_stateObj.idx, {'count': count, 'event': event}))
    state_graph.add_edges_from(edges)
    print("[+] Generated graph with %s states, %s events in %s seconds." % (state_graph.number_of_nodes(), state_graph.number_of_edges(), timer() - start_time))
    return state_graph, lcaStateObj

def get_state_paths(state_graph, sampleStateObj):
    '''
    Networkx documentation : https://networkx.github.io/documentation/stable/index.html
    '''
    start_time = timer()
    # out_edges_counter = Counter([source for source, sink, count in state_graph.edges.data('count') for i in range(count)]) # previous approach

    # Get total outgoing edges by event_type for each node 
    total_by_event_by_node_idx = defaultdict(Counter)
    for source, sink, data in state_graph.edges.data():
        total_by_event_by_node_idx[source][data['event']] += data['count']
    node_idx_by_meta = {meta: node_idx for node_idx, meta in nx.get_node_attributes(state_graph, 'meta').items() if meta}
    mutype_by_node_idx = {node_idx: mutype for node_idx, mutype in nx.get_node_attributes(state_graph, 'mutype').items()}
    header = ['path_idx', 'node', 'event', 'event_count', 'C_ancestor', 'C_derived', 'Ms', 'Es', 'hetA', 'fixed', 'hetB', 'hetAB']
    paths = []
    path_count = 0
    for path_idx, path in enumerate(nx.all_simple_paths(state_graph, source=node_idx_by_meta['source'], target=node_idx_by_meta['sink'])):
        path_count += 1
        step_idx = 0
        for source, sink in pairwise(path):
            line = []
            line.append(path_idx)
            line.append(source)
            line.append(state_graph[source][sink]['event'])
            line.append(state_graph[source][sink]['count'])
            line.append(total_by_event_by_node_idx[source].get('C_ancestor', 0))
            line.append(total_by_event_by_node_idx[source].get('C_derived', 0))
            line.append(total_by_event_by_node_idx[source]['M'])
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

def state_from_lol(pops, ploidy=1):
    state = []
    for pop in pops:
        if len(pop) == 0:
            state.append((),)
        else:
            sub_state = []
            for s in pop:
                for p in range(ploidy):
                    sub_state.append(s)
            state.append(tuple(sorted(sub_state)))
    return tuple(state)

###############################################################################

################################### Classes ###################################

class StateObj(object):
    def __init__(self, idx=None, state='', res_idx=None, mig_idx=None):
        self.idx = idx
        self.state = state   # is a tuple of 'P' tuples (P = Number of populations)
        self.resident_state = None if res_idx is None else state[res_idx]
        self.res_idx = res_idx
        self.migrant_state = None if mig_idx is None else state[mig_idx]
        self.mig_idx = mig_idx
        self.lineages = (str(self.state).replace("'", " ").replace("(", " ").replace(")", " ").replace(",", " ").split())

    def __eq__(self, other):
        if isinstance(other, StateObj):
            return self.state == other.state

    def __str__(self):
        if self.mig_idx is None and self.res_idx is None:
            return str([list(pop) for pop in self.state]).replace('\'', '')
        else:
            mig_string = "{%s}" % ",".join(list(self.migrant_state))
            res_string = "[%s]" % ",".join(list(self.resident_state))
            return "[%s, %s]" % (mig_string, res_string)

    def __hash__(self):
        return hash(self.state)

    def get_list(self):
        return [list(state) for state in self.state]

    def get_migrant_list(self):
        return list(self.migrant_state)    

    def get_resident_list(self):
        return list(self.resident_state)

class ParameterObj(object):
    def __init__(self, args):
        self.ploidy = int(args['--ploidy'])
        self.migration = args['--migration']
        self.plot_mutypes = args['--plot_mutypes']
        self.no_graph_labels = args['--no_graph_labels']
        self.exodus = args['--exodus']
        self.reverse_exodus = args['--reverse_exodus']
        self.out_prefix = args['--out_prefix']
        self.samples_input = literal_eval(args['--samples']) # can't be sorted because can be list + set
        self.migration_possible = any([self.migration, self.exodus, self.reverse_exodus])
        self.count_pops = len(self.samples_input)
        self.count_samples = sum([len(pop) for pop in self.samples_input])
        self.migrant_idx = None
        self.resident_idx = None

    def setup_model(self):
        sample = self.samples_input
        print("[=] Samples =", self.samples_input)
        if self.migration_possible: # either --migration or --exodus or --reverse_exodus
            self.polarise_pops()
        self.sample_state = state_from_lol(sample, ploidy=self.ploidy)
        print("[=] Sample state =", self.sample_state)
        sampleStateObj = StateObj(idx=0, state=self.sample_state, res_idx=self.resident_idx, mig_idx=self.migrant_idx)
        lca_state = self.get_lca_state()
        print("[=] LCA state =", lca_state)
        lcaStateObj = StateObj(idx=None, state=lca_state, res_idx=self.resident_idx, mig_idx=self.migrant_idx)
        return sampleStateObj, lcaStateObj
    
    def get_lca_state(self):
        if not self.migration_possible:
            # assumes that [[ ... ]]
            return state_from_lol([["".join(sorted(self.sample_state[0]))]])
        else:
            migrants = sorted(self.sample_state[self.migrant_idx])
            residents = sorted(self.sample_state[self.resident_idx])
            mrcas = [None, None]
            if self.reverse_exodus:
                mrcas[self.migrant_idx] = ["".join(sorted(migrants + residents))]
                mrcas[self.resident_idx] = []
            else:
                mrcas[self.migrant_idx] = []
                mrcas[self.resident_idx] = ["".join(sorted(migrants + residents))]
            return state_from_lol(mrcas)

    def polarise_pops(self):
        if not self.count_pops == 2:
            exit("[X] '--migration' and/or '--exodus' only implemented for 2 populations. You provided %s population(s)!" % self.count_pops)
        if not self.count_samples == 2:
            exit("[X] '--migration' and/or '--exodus' only implemented for 2 samples. You provided %s sample(s)!" % self.count_samples)
        if not any([isinstance(sample, set) for sample in self.samples_input]):
            exit("[X] Must specify migrant population using '{}' when using '--migration' and/or '--exodus'!")
        if isinstance(self.samples_input[0], set):
            self.migrant_idx, self.resident_idx = 0, 1
        else:
            self.migrant_idx, self.resident_idx = 1, 0
        if self.migration:
            print("[=] Migration : ", self.samples_input[self.migrant_idx], "->", self.samples_input[self.resident_idx])
        if self.exodus:
            print("[=] Exodus    : ", self.samples_input[self.migrant_idx], "->", self.samples_input[self.resident_idx])
        elif self.reverse_exodus:
            print("[=] Exodus    : ", self.samples_input[self.migrant_idx], "<-", self.samples_input[self.resident_idx])
        else:
            pass
        
    def get_basename(self):
        models = OrderedDict({
            'M': self.migration, 
            'E': self.exodus, 
            'R': self.reverse_exodus
        })
        model_string = "+".join([model_id for model_id, status in models.items() if status])
        return "%s.%s_samples.%s_pop.%s_ploidy.%s" % (self.out_prefix, self.count_samples, self.count_pops, self.ploidy, model_string)

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
            node['state'] = str(node.get('state', None))
        nx.write_gml(state_graph, out_f)
        print("[+] Written graph into file %s in %s seconds." % (out_f, timer() - start_time))
        return out_f

    def read_graph(self, infile):
        start_time = timer()
        state_graph = nx.read_gml(infile)
        for idx, node in state_graph.nodes(data=True):
            node['state'] = literal_eval(node['state'])
            if 'pos' in node:
                node['pos'] = tuple(node['pos'])
            # print(idx, node)
        print("[+] Read graph from file %s in %s seconds." % (infile, timer() - start_time))
        return state_graph

    def plot_graph(self, state_graph):
        start_time = timer()
        out_f = "%s.mutypes_%s.no_graph_labels_%s.graph.png" % (self.get_basename(), self.plot_mutypes, self.no_graph_labels)
        plot_state_graph(state_graph, out_f, mutypes_flag=self.plot_mutypes, no_graph_labels=self.no_graph_labels)
        print("[+] Plotted graph into file %s in %s seconds." % (out_f, timer() - start_time))
        return out_f

###############################################################################

################################### Main ######################################

def main():
    try:
        main_time = timer()
        args = docopt(__doc__)
        parameterObj = ParameterObj(args)
        sampleStateObj, lcaStateObj = parameterObj.setup_model()
        state_graph, lcaStateObj = build_state_graph(sampleStateObj, lcaStateObj, parameterObj)
        header, paths = get_state_paths(state_graph, sampleStateObj)
        parameterObj.write_paths(header, paths)
        parameterObj.plot_graph(state_graph)
        # if state_graph.number_of_nodes() < 70:
        #     parameterObj.plot_graph(state_graph)
        # else:
        #     print("[!] Too many nodes (%s). No plot is generated." % state_graph.number_of_nodes())
        parameterObj.write_graph(state_graph)
        #state_graph = parameterObj.read_graph(graph_file)
        
    except KeyboardInterrupt:
        stderr.write("\n[X] Interrupted by user after %s seconds!\n" % (timer() - main_time))
        exit(-1)

###############################################################################

if __name__ == '__main__':
    main()