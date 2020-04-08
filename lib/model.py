import itertools
import collections
import networkx as nx
import matplotlib as mat
mat.use("agg")
from tqdm import tqdm
import ast
import re
import sys
import lib.functions

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

'''

'a' = 'a1' 
'a' = 'a2'
'b' = 'b1'
'b' = 'b2'
'aa' = 'a1a2'
'bb' = 'b1b2'
'ab' = 'a1b1'
'ab' = 'a1b2'
'ab' = 'a2b1'
'ab' = 'a2b2'
'aab' = 'a1a2b1'
'aab' = 'a1a2b2'
'abb' = 'a1b1b2'
'abb' = 'a2b1b2'
'aabb' = 'a1a2b1b2'


Needs debugging output as in:

Event_label: event_type + populations

State: {'A': [], 'B': ['a1','a1','b1+b1']}
    [J_A&B] 
    [M_A>B]
    [C_A]
    [C_B]
    [C_A&B]
'''
MUTYPES = ['hetA', 'fixed', 'hetB', 'hetAB']

###############################################################################

def flat_tuple(obj):
    for x in obj:
        if isinstance(x, str):
            yield x
        elif isinstance(x, tuple):
            yield from flat_tuple(x)
        else:
            raise TypeError

def joins(string, events=None):
    if events is None:
        events = []
    for substring in string:
        if isinstance(substring, tuple):
            if len(substring) == 2:
                joins(substring, events)
            else:
                for pop_1, pop_2 in itertools.combinations(sorted(substring), 2):
                    joins((pop_1, pop_2), events)
    if string:
        events.append("%s" % "&".join(sorted(list(flat_tuple(string)))))
    return events

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
        if len(self.list) == 1 and self.list[0] == '':
            return 0
        return len(self.list)

    def __hash__(self):
        return hash(self.string)

    def __getitem__(self, i):
        return self.list[i]

    def __add__(self, other):
        if self.has_lineages() and other.has_lineages():
            return PopObj(sorted(self.list + other.list))
        elif self.has_lineages():
            return PopObj(self.list)
        elif other.has_lineages():
            return PopObj(other.list)
        else:
            sys.exit("[X] Should not happen...")

    def __radd__(self, other):
        if other == 0:
            return self
        else:
            return self.__add__(other)

    def has_lineages(self):
        if len(self.list) == 0 or (len(self.list) == 1 and self.list[0] == ''):
            return False
        return True 

class StateObj(object):
    def __init__(self, value):
        print("[INIT] %s" % value)
        self.popObj_by_pop_id = self._get_dict(value)
        self.label = self._get_label(value)
        self.graph_label = self._get_graph_label()
        self.lineages = self._get_lineages()
        self.count_by_ancestor_labels_by_event = collections.defaultdict(lambda: collections.Counter())
        self.lineage_counter = self._get_lineage_counter()

    def _get_lineages(self):
        lineages = []
        for pop_id, popObj in sorted(self.popObj_by_pop_id.items()):
            if popObj.has_lineages():
                lineages += popObj.list
        return ",".join(sorted(lineages))

    def _get_lineage_counter(self):
        lineage_counter = []
        for popObj in self.popObj_by_pop_id.values():
            for lineage in popObj.list:
                if not lineage == "": 
                    lineage_counter.append(lineage.replace("+", ""))            
        return collections.Counter(lineage_counter)

    def _get_graph_label(self):
        return ";".join(["%s=%s" % (pop_id, popObj) for pop_id, popObj in sorted(self.popObj_by_pop_id.items()) if popObj.has_lineages()])
    
    def _get_label(self, value):
        if isinstance(value, dict):
            return ";".join(["%s=%s" % (pop_id, popObj) if popObj.has_lineages() else "%s=%s" % (pop_id, []) for pop_id, popObj in sorted(self.popObj_by_pop_id.items())])
        elif isinstance(value, str):
            return value
        else:
            sys.exit("[X] value must be 'str' or 'dict', is: %r" % value)

    def _get_dict(self, value):
        print(value, type(value))
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
        return str(self.__dict__)

    def __repr__(self):
        return str(self.__dict__)

    def __hash__(self):
        return hash(self.label)

    def __eq__(self, other):
        return self.label == other.label

    def __len__(self):
        return len(self.popObj_by_pop_id)

    def get_ancestors(self, parameterObj, lca_label):
        if not self.lineages == lca_label:
            for event in sorted(parameterObj.set_of_events):
                if event.startswith("C"):
                    self._get_coalescence_ancestors(event)
                elif event.startswith("M"):
                    self._get_migration_ancestors(event)
                else:
                    self._get_join_ancestors(event)

    def _get_join_ancestors(self, event):
        print("# join...")
        pop_1_id, pop_2_id = event.replace("J=", "").split("&") 
        if pop_1_id in self.popObj_by_pop_id and pop_2_id in self.popObj_by_pop_id:
            ancestor_popObj_by_pop_id = {pop_id: popObj for pop_id, popObj in self.popObj_by_pop_id.items() if pop_id not in set([pop_1_id, pop_2_id])}
            ancestor_popObj_by_pop_id[event] = self.popObj_by_pop_id[pop_1_id] + self.popObj_by_pop_id[pop_2_id] # summation of PopObjs
            #print("self.popObj_by_pop_id[pop_1_id]", self.popObj_by_pop_id[pop_1_id], "self.popObj_by_pop_id[pop_2_id]", self.popObj_by_pop_id[pop_2_id])
            ancestor_stateObj = StateObj(ancestor_popObj_by_pop_id)
            #print("[J]", ancestor_stateObj)
            self.count_by_ancestor_labels_by_event[event][ancestor_stateObj.label] += 1
    
    def _get_coalescence_ancestors(self, event):
        print("# coalesce...")
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
        print("# migrate...")
        pop_1_id, pop_2_id = event.replace("M=", "").split(">")
        if all([pop_id in self.popObj_by_pop_id for pop_id in [pop_1_id, pop_2_id]]) == True:
            print(event, pop_1_id, pop_2_id, self.popObj_by_pop_id[pop_1_id])
            if len(self.popObj_by_pop_id[pop_1_id]):
                for i, allele in enumerate(self.popObj_by_pop_id[pop_1_id]):
                    ancestor_popObj_by_pop_id = {pop_id: popObj for pop_id, popObj in self.popObj_by_pop_id.items() if pop_id not in set([pop_1_id, pop_2_id])}
                    source_popObj = PopObj([allele for k, allele in enumerate(self.popObj_by_pop_id[pop_1_id]) if not k == i])
                    ancestor_popObj_by_pop_id[pop_1_id] = source_popObj
                    ancestor_popObj_by_pop_id[pop_2_id] = self.popObj_by_pop_id[pop_2_id] + PopObj([allele])
                    ancestor_stateObj = StateObj(ancestor_popObj_by_pop_id)
                    self.count_by_ancestor_labels_by_event[event][ancestor_stateObj.label] += 1

class StateGraph(object):
    def __init__(self):
        # default is DiGraph
        self.graph = nx.MultiDiGraph(label='', meta='', mutype='', key='') 
        self.stateObj_by_node_idx = {}
        self.node_idx_by_label = {}
        self.label_by_node_idx = {}
        self.origin_idx_by_label = {}
        self.sample_label = None
        self.lca_label = None
        self.events_counter = collections.Counter()
        self.events_counter_by_node_idx = collections.defaultdict(collections.Counter)
        self.edges_by_event = collections.defaultdict(list) # for colours when plotting
        self.set_of_lineages = set()
        self.pop_ids = set()

    def get_origins(self):
        self.origin_idx_by_label = {self.graph.nodes[node_idx]['label']: node_idx for node_idx, out_degree in self.graph.out_degree() if out_degree == 0}

    def write_model(self, parameterObj):
        print("[+] Calculating model ...")
        #print(self.events_counter)
        events = sorted(self.events_counter.keys())
        lineages = sorted([lineage for lineage in self.set_of_lineages if not lineage == self.lca_label])
        header = ['path_idx', 'node_id', 'LCA', 'label', 'event', 'count'] + ['total'] + events + lineages 
        rows = []
        path_counter_by_origin = collections.Counter()
        max_length_label = max([len(label) for label in self.origin_idx_by_label.keys()])
        path_idx = 0
        for origin_label, origin_idx in sorted(self.origin_idx_by_label.items()):
            paths = [path for path in nx.all_simple_paths(self.graph, source=0, target=origin_idx)]
            #print('idx_paths', idx_paths)
            for path in tqdm(paths, total=len(paths), desc="[+]\tPaths to %s" % origin_label.rjust(max_length_label), unit='', ncols=100, unit_scale=True):
                #print("\n############## PATH ", path_idx, path)
                path_counter_by_origin[origin_label] += 1
                path_counter_by_origin['all'] += 1
                for source_idx, sink_idx in pairwise(path):
                    current_event = self.graph[source_idx][sink_idx][0]['event']
                    current_count = self.graph[source_idx][sink_idx][0]['count']
                    row = [path_idx, sink_idx, origin_idx, self.label_by_node_idx[source_idx], current_event, current_count, sum(self.events_counter_by_node_idx[source_idx].values())]
                    for event in events:
                        row.append(self.events_counter_by_node_idx[source_idx][event])
                    stateObj = self.stateObj_by_node_idx[source_idx]
                    
                    for lineage in lineages:
                        row.append(stateObj.lineage_counter[lineage])
                    rows.append(row)
                last_row = [path_idx, sink_idx, origin_idx, self.label_by_node_idx[source_idx], 'LCA', 0, 0] + ([0] * len(events)) + ([0] * len(lineages))
                rows.append(last_row)
                path_idx += 1
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
        G.get_node(0).attr['rank'] = 'source'
        for origin_label, origin_idx in self.origin_idx_by_label.items():
            G.get_node(origin_idx).attr['fontcolor'] = 'black'
        colours = ['black', 'slategray', 'lightgrey', 'grey50', 'grey25', 'yellowgreen']
        colour_by_pops = {pop: colours[idx] for idx, pop in enumerate(sorted(self.pop_ids))}
        for node_idx, label in self.label_by_node_idx.items():
            newlabel = ['<<table border="0" cellspacing="1">']
            bg_color = 1
            if label in self.origin_idx_by_label:
                bg_color = 0 #'bgcolor="black"'
            for pop in label.split(";"):
                newlabel.append('<tr><td port="port%s" border="%s" style="solid" color="%s">%s</td></tr>' % (node_idx, bg_color, colour_by_pops[pop.split("=")[0]], pop.replace("&", "&amp;")))
            newlabel.append('</table>>')
            #G.get_node(node_idx).attr['label']="|".join(["<f%s> %s" % (idx, pop_string) for idx, pop_string in enumerate(label.split(";"))])
            G.get_node(node_idx).attr['label']="".join(newlabel)
        for event, edgelist in self.edges_by_event.items():
            if event.startswith("C="):
                colour = 'dodgerblue'
            elif event.startswith("M="):
                colour = 'orange'
            else:
                colour = 'grey'
            for source, sink, count in edgelist:
                G.get_edge(source, sink).attr['color'] = colour
                G.get_edge(source, sink).attr['label'] = count
                G.get_edge(source, sink).attr['fontcolor'] = colour
                G.get_edge(source, sink).attr['penwidth'] = count
        G.layout()
        out_f = "%s.pdf" % parameterObj.out_prefix
        print("[>]\tCreated: %r" % out_f)
        G.draw(out_f, prog='dot', args="-Gratio='compress' -Gsize='7.5,10'")

    def is_directed_acyclic_graph(self):
        return nx.is_directed_acyclic_graph(self.graph)

    def add_state(self, stateObj):
        if not stateObj.label in self.node_idx_by_label:
            #print("[ADDING]", stateObj.label)
            node_idx = len(self.node_idx_by_label)
            self.node_idx_by_label[stateObj.label] = node_idx
            self.stateObj_by_node_idx[node_idx] = stateObj
            self.label_by_node_idx[node_idx] = stateObj.graph_label.replace("'", "") 
            self.graph.add_node(node_idx, label=stateObj.graph_label.replace("'", ""), mutype='', key='')
            for lineage in stateObj.lineage_counter.keys():
                self.set_of_lineages.add(lineage)
            for pop_id in stateObj.popObj_by_pop_id.keys():
                self.pop_ids.add(pop_id)
            if node_idx == 0:
                self.lca_label = "+".join(sorted(["+".join(popObj.list) for pop_id, popObj in stateObj.popObj_by_pop_id.items() if popObj.has_lineages()]))
                self.sample_label = stateObj.label
        #else:
            #print("[NOT ADDING]", stateObj.label)

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
            state_graph.add_state(focal_stateObj)
            #print("[F]", focal_stateObj, state_graph.lca_label)
            focal_stateObj.get_ancestors(parameterObj, state_graph.lca_label)
            for event in focal_stateObj.count_by_ancestor_labels_by_event.keys():
                for ancestor_label, count in focal_stateObj.count_by_ancestor_labels_by_event[event].items():
                    if not ancestor_label in state_graph.node_idx_by_label:
                        ancestor_stateObj = StateObj(ancestor_label)
                        #print("[A]", event, ancestor_stateObj)
                        state_graph.add_state(ancestor_stateObj)
                        queue.append(ancestor_stateObj)
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
        self.samples_by_pop_id = {pop_id: int(samples) for pop_id, samples in zip(args['--pop_ids'].split(','), args['--samples'].split(','))}
        self.out_prefix = args['--out_prefix']
        self.join_events = self._parse_join_events(args['--join_string']) if args['--join_string'] else []
        self.migration_events = self._parse_migration_events(args['--migration_string']) if args['--migration_string'] else []
        self.sample_stateObj = self._get_sample_stateObj()
        self.set_of_events = self._get_set_of_events() 
        self.model_output = False if args['--nomodel'] else True 
        self.graph_output = False if args['--nograph'] else True
        print(self.__dict__) 

    def _get_set_of_events(self):
        return set(["C=%s" % pop_id for pop_id in self.pop_ids + self.join_events] + self.join_events + self.migration_events)

    def _get_sample_stateObj(self):
        pop_dict = {}
        for pop_id in self.pop_ids:
            pop_list = []
            for i in range(self.samples_by_pop_id[pop_id]):
                for j in range(self.ploidy):
                    pop_list.append("%s%s" % (pop_id.lower(), i+1))
            pop_dict[pop_id.upper()] = PopObj(pop_list)
        return StateObj(pop_dict)

    def _parse_join_events(self, join_string):
        char_counter = collections.Counter(join_string)
        if not all([char_counter[pop_id] <= 1 for pop_id in self.pop_ids]) == True:
            return "[X] Duplicated population IDs"
        print(joins(ast.literal_eval(re.sub(r'([a-zA-Z]+)',r'"\1"', join_string))))
        return joins(ast.literal_eval(re.sub(r'([a-zA-Z]+)',r'"\1"', join_string)))

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
