import sys
import os
import yaml
import pandas
import numpy
import multiprocessing
import psutil
from contextlib import contextmanager
from collections import OrderedDict, defaultdict, Counter, deque
from tqdm import tqdm
from math import floor
from cyvcf2 import VCF
#from datetime import datetime
from natsort import natsorted
#import operator
from itertools import combinations, product, chain
import copy
import pickle
#from scipy.interpolate import BSpline

########################### PLOT
import matplotlib as mat
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.ticker import MaxNLocator
from matplotlib import cm
mat.use("agg")
import matplotlib.pyplot as plt
import seaborn as sns
plt.style.use('seaborn-white')
#mat.rcParams['font'] = 16
#print(mat.rcParams.keys())
mat.rcParams['text.color'] = 'grey'
mat.rcParams['axes.edgecolor'] = 'lightgrey'
mat.rcParams['xtick.color'] = 'grey'
mat.rcParams['ytick.color'] = 'grey'
mat.rcParams['grid.color'] = 'lightgrey'
mat.rcParams['font.family'] = 'sans-serif'
legendfontsize = 18
axestickfontsize = 16
axeslabelfontsize = 16
mat.rcParams['xtick.labelsize'] = axestickfontsize
mat.rcParams['ytick.labelsize'] = axestickfontsize
mat.rcParams['figure.frameon'] = False
mat.rcParams['axes.grid'] = False

'''
ax.grid(True)
plt.gca().spines["top"].set_visible(False)  
plt.gca().spines["right"].set_visible(False)
ax.set_frame_on(False)
ax.lines[0].set_visible(False)

[Notes]

- 
- test dataset (only chr1 & chr3)

> awk '$1=="chr1" || $1=="chr3"' ../../input/Hmel2.chromTransfers.2016_09_20.txt | cut -f4 | sort | uniq > chr1_chr3.contigs.txt

# genome file (71 contigs)
> grep -wFf chr1_chr3.contigs.txt ../real/hmel.autosomes.genomefile > test.autosomes.genomefile

# BED (1642471 intervals on 71 contigs)
> grep -wFf chr1_chr3.contigs.txt ../real/hmel.multiinter.samples_as_string.only_intergenic.bed > test.multiinter.samples_as_string.only_intergenic.bed

[To Do]

- change pi_A/B to pop name
- Fst/dxy contig means if longer than x windows

- write proper summary tables
    - pair IDs
    - piA/piB
    - look at div/pi computations and how to list them
- histogram of block/window Fst's








- Plot Fst of pairs in B and AllCallable
- calculate/save/print Pi for each sample  

- report stats on average hetA, hetAb, hetB, fixed
 
- make a tree of pairwise dxy : phylip, which format?
- astral/dadi/twist/discovista
- change algoA so that it gets invoked if BLOCKLENGTH ≥ MIN_INTERVAL_LENGTH
- Check how it works with 1bp blocks?
- Check 64b, 128b
- check difference(BED,BED_ALLCALLABLE) for how many regions covered by each
- does bias only affect intergenic? do same with genic

- Mean_Fst of algoB should be closer to lohse2016 thank all_callable. The more samples are added, the greater the upward bias under all_callable.
=> scatterplot with Fst values.

'''

#def format_bases(*args, block_length):
#    return "%0.0f" % (args[0] * block_length)

class PlotScatterObj():
    def __init__(self, parameterObj, name, x_label, x_key, y_label, y_key, df):
        self.parameterObj = parameterObj
        self.name = name
        self.y_label = y_label
        self.x_label = x_label
        self.x = df[x_key]
        self.y = df[y_key]

    def plot(self, alpha, points):
        fig = plt.figure(figsize=(12,12), dpi=400, frameon=False)
        ax = fig.add_subplot(111)
        #ax.set_ylim((-0.05, 1))
        #ax.set_xlim((0, categories + 1))
        #ax.set_xscale("log")
        plt.ylabel(self.y_label, fontsize=axeslabelfontsize)
        plt.xlabel(self.x_label, fontsize=axeslabelfontsize)
        ax.plot(self.x, self.y, 'o', alpha=alpha, label=points, color='cornflowerblue', markersize=4)
        plt.legend(fontsize = legendfontsize, frameon=True)
        #ax.set_aspect(1.)
        divider = make_axes_locatable(ax)
        axHistx = divider.append_axes("top", 1.2, pad=0.2, sharex=ax)
        axHisty = divider.append_axes("right", 1.2, pad=0.2, sharey=ax)
        axHistx.xaxis.set_tick_params(labelbottom=False)
        axHisty.yaxis.set_tick_params(labelleft=False)
        axHistx.hist(self.x, bins=100, color='cornflowerblue')
        axHisty.hist(self.y, bins=100, color='cornflowerblue', orientation='horizontal')
        fn = "%s.png" % (self.name) 
        if self.parameterObj:
            fn = "%s.%s.png" % (self.parameterObj.outprefix, self.name) 
        plt.tight_layout()
        fig.savefig(fn, format="png")
        plt.close(fig)
        return fn

class PlotSwarmObj():
    def __init__(self, name, x_label, x_key, y_label, y_key, df):
        self.name = name
        self.y_label = y_label
        self.x_label = x_label
        self.x = df[x_key]
        self.y = df[y_key]
        self.df = df

    def plot(self, alpha, points):
        fig = plt.figure(figsize=(12,12), dpi=400, frameon=False)
        plt.ylabel(self.y_label, fontsize=axeslabelfontsize)
        plt.xlabel(self.x_label, fontsize=axeslabelfontsize)
        sns.swarmplot(x=self.x, y=self.y)
        sns.boxplot(x=self.x, y=self.y,
            showcaps=False,boxprops={'facecolor':'None', 'edgecolor':'None'},
            showfliers=False,whiskerprops={'linewidth':0})
        fn = "%s.png" % (self.name) 
        plt.gca().spines["top"].set_visible(False)  
        plt.gca().spines["right"].set_visible(False)
        fig.savefig(fn, format="png")
        plt.close(fig)
        return fn

class PlotHeatmapObj():
    def __init__(self, parameterObj, name, title, df):
        self.parameterObj = parameterObj
        self.title = title
        self.df = df
        self.name = name

    def get_data(self):
        for population_idx, population_id in self.parameterObj.population_by_population_idx.items():
            self.df = self.df.assign(**{population_id : pandas.Series(self.df['pair_idx'].apply(get_pair_id, args=(population_idx, self.parameterObj,)))})
            #self.df.loc[:,population_id] = self.df['pair_idx'].apply(get_pair_id, args=(population_idx, self.parameterObj,))
        df_pivot = self.df.pivot_table(self.name, *sorted(self.parameterObj.population_by_population_idx.values()), fill_value=0.0)
        return df_pivot

    def plot(self):
        fig = plt.figure(figsize=(10,10), dpi=400, frameon=False)
        #
        ax = fig.add_subplot(111)
        df_pivot = self.get_data()
        A_ids = df_pivot.columns.astype(str).tolist()
        B_ids = df_pivot.index.astype(str).tolist()
        vmax = df_pivot.values.max()
        vmin = df_pivot.values.min()
        vrange = numpy.linspace(vmin, vmax, 10, endpoint=True)
        im = ax.imshow(df_pivot, cmap=cm.Oranges, origin='lower')
        ax.set_xticks(numpy.arange(len(A_ids)))
        ax.set_yticks(numpy.arange(len(B_ids)))
        ax.set_xticks(numpy.arange(-.5, len(A_ids), 1), minor=True)
        ax.set_yticks(numpy.arange(-.5, len(B_ids), 1), minor=True)
        ax.grid(which='minor', color='w', linestyle='-', linewidth=2)
        ax.set_xticklabels(A_ids)
        ax.set_yticklabels(B_ids)
        ax.grid(False)
        plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
             rotation_mode="anchor")
        for i in range(len(A_ids)):
            for j in range(len(B_ids)):
                ax.text(i, j, "%.2f" % (df_pivot[A_ids[i]][B_ids[j]] / vmax), ha="center", va="center", color="black", fontsize=14)
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.1)
        cbar = fig.colorbar(im, cax=cax, ticks=vrange)
        cbar.ax.set_title(self.name, loc='left', pad=10)
        cbar.ax.tick_params(labelsize=12)
        cbar.outline.set_visible(False)
        ax.set_title(self.title, pad=9)
        fn = "%s.%s.png" % (self.parameterObj.outprefix, self.name) 
        plt.tight_layout()
        fig.savefig(fn, format="png")
        plt.close(fig)
        return fn 

class PlotHistObj():
    def __init__(self, parameterObj, name, x_label, y_label, df):
        self.parameterObj = parameterObj
        self.name = name
        self.x_label = x_label
        self.y_label = y_label
        self.df = df
        self.cm = plt.get_cmap('Set2')

    def get_data(self, column_ids, desired_bin_size):
        x_min, x_max = 0, 0
        x_by_column_id = {}
        for idx, column_id in enumerate(column_ids):
            data = self.df[lambda df: df.columns[idx]]
            _min = min(data)
            _max = max(data)
            if not x_min or x_min > _min:
                x_min = _min
            if not x_max or x_max < _max:
                x_max = _max
            x_by_column_id[column_id] = data
        min_boundary = -1.0 * (x_min % desired_bin_size - x_min)
        max_boundary = x_max - x_max % desired_bin_size + desired_bin_size
        n_bins = int((max_boundary - min_boundary) / desired_bin_size) + 1
        bins = numpy.linspace(min_boundary, max_boundary, n_bins)
        xlim = (x_min - 0.1, x_max + 0.1)
        return (x_by_column_id, xlim, bins)

    def plot(self):
        fig = plt.figure(figsize=(16,6), dpi=400, frameon=False)
        ax = fig.add_subplot(111)
        column_ids = self.df.columns.tolist()
        x_by_column_id, xlim, bins = self.get_data(column_ids, 0.01)
        colour = 'cornflowerblue'
        for idx, column_id in enumerate(column_ids):
            if len(column_ids) > 1:
                colour = self.cm(idx / len(column_ids))
            _lw=1
            counts, bins, patches = ax.hist(x=x_by_column_id[column_id], label="Distribution of %s" % (column_id), bins=bins, alpha=0.9, lw=_lw, facecolor=colour, edgecolor='white', align='left')
        ax.set_xlim(xlim)
        plt.legend(fontsize = legendfontsize, frameon=True)
        plt.title("Distribution of %s" % (self.name))
        plt.ylabel(self.y_label, fontsize=axeslabelfontsize)
        plt.xlabel(self.x_label, fontsize=axeslabelfontsize)
        fn = "%s.%s.png" % (self.parameterObj.outprefix, self.name) 
        #plt.gca().spines["top"].set_visible(False)  
        #plt.gca().spines["right"].set_visible(False)
        plt.tight_layout()
        fig.savefig(fn, format="png")
        plt.close(fig)
        return fn

class PlotGenomeObj():
    def __init__(self, parameterObj, name, df, length_by_seq_id, subplots, by_population):
        self.parameterObj = parameterObj
        self.name = name
        self.df = df
        self.length_by_seq_id = length_by_seq_id
        self.subplots = subplots
        self.cm = plt.get_cmap('Set2')
        self.by_population = by_population

    def get_data(self, column_ids, row_id):
        offset_by_contig_id = {}
        offset = 0
        x_boundaries = []
        for seq_id, length in self.length_by_seq_id.items():
            offset_by_contig_id[seq_id] = offset
            x_boundaries.append(offset)
            offset += length
        ys_by_column_id_by_window_id = self.df.set_index('window_id').T.to_dict()
        y_by_contig_id_by_column_id = defaultdict(lambda: defaultdict(list))
        x_by_contig_id = defaultdict(list)
        for window_id in natsorted(ys_by_column_id_by_window_id.keys()):
            contig_id, start, end = window_id.split("_")
            position = int(start) + ((int(end) - int(start)) / 2) + offset_by_contig_id[contig_id]
            x_by_contig_id[contig_id].append(position)
            for column_id in column_ids:
                y_by_contig_id_by_column_id[column_id][contig_id].append(float(ys_by_column_id_by_window_id[window_id][column_id]))
        return (x_by_contig_id, y_by_contig_id_by_column_id, x_boundaries)

    def plot(self):
        column_ids = self.df.columns.tolist()[1:]
        row_id = self.df[self.df.columns[0]].tolist()
        fig, axarr = plt.subplots(1, 1, figsize=(16,6), dpi=200, frameon=False)
        if self.subplots:
            fig, axarr = plt.subplots(len(column_ids), 1, sharex=True, figsize=(24,(len(column_ids) * 6)), dpi=200)
        x_by_contig_id, y_by_contig_id_by_column_id, x_boundaries = self.get_data(column_ids, row_id)
        #print(x_by_contig_id, y_by_contig_id_by_column_id, x_boundaries)
        plt.gca().spines["top"].set_visible(False)  
        plt.gca().spines["right"].set_visible(False)
        max_y = 0.0
        min_y = 1.0
        _handles, _labels = [], []
        for idx, column_id in enumerate(column_ids):
            i = 0
            y_list = []
            for seq_id, length in self.length_by_seq_id.items():
                x, y = x_by_contig_id[seq_id], y_by_contig_id_by_column_id[column_id][seq_id]
                #print(seq_id, length)
                if y:
                    y_list.extend(y)
                    _max_y = max(y)
                    if max_y < _max_y:
                        max_y = _max_y
                    if min_y > min(y):
                        min_y = min(y)
                    if self.by_population:
                        colour = self.parameterObj.colour_by_population[self.parameterObj.population_by_sample_id[column_id]]
                        axarr.plot(x, y, color=colour, alpha=0.5, marker='o', markersize=0.2, linewidth=0)
                        axarr.vlines(x_boundaries, 0, max_y, colors=['lightgrey'], linestyles='dashed', linewidth=1)
                        axarr.set(ylabel = column_id)
                        axarr.spines['right'].set_visible(False)
                        axarr.spines['top'].set_visible(False)
                    elif self.subplots:
                        if self.name == 'missing_multiallelic':
                            colour = 'cornflowerblue' if i % 2 else 'yellowgreen'
                        else:
                            colour = 'orange' if i % 2 else 'mediumorchid'
                        axarr[idx].plot(x, y, color=colour, alpha=0.5, marker='o', markersize=0, linewidth=2)
                        axarr[idx].plot(x, y, color='black', alpha=0.5, marker='o', markersize=0.5, linewidth=0)
                        axarr[idx].vlines(x_boundaries, 0, max_y, colors=['lightgrey'], linestyles='solid')
                        axarr[idx].set(ylabel = column_id)
                        axarr[idx].spines['right'].set_visible(False)
                        axarr[idx].spines['top'].set_visible(False)
                    else:
                        colour = self.cm(idx / len(column_ids))
                        axarr.plot(x, y, color=colour, alpha=0.5, marker='o', markersize=0, linewidth=1)
                        axarr.plot(x, y, color=colour, alpha=1, marker='o', markersize=0.2, linewidth=0)
                        axarr.vlines(x_boundaries, 0, max_y, colors=['lightgrey'], linestyles='dashed', linewidth=1)
                        axarr.set(ylabel = column_id)
                        axarr.spines['right'].set_visible(False)
                        axarr.spines['top'].set_visible(False)
                i += 1
            if self.subplots:
                axarr[idx].axhline(y=numpy.mean(y_list), xmin=0, xmax=1, color='darkgrey', linestyle='--', linewidth=2)
                axarr[idx].set_ylim(min_y, max_y)
            else:
                axarr.set_ylim(min_y - (min_y / 100) , max_y + (max_y / 100))
        if self.by_population:
            for population, colour in self.parameterObj.colour_by_population.items():
                _handles.append(mat.lines.Line2D([0], [0], c=colour, lw=4))
                _labels.append(population)
            plt.ylabel(self.name, fontsize=axeslabelfontsize)
        elif self.subplots:
            plt.ylabel(self.name, fontsize=axeslabelfontsize)
            pass
        else:
            plt.ylabel(self.name, fontsize=axeslabelfontsize)
            for idx, column_id in enumerate(column_ids):
                colour = self.cm(idx / len(column_ids))
                _handles.append(mat.lines.Line2D([0], [0], c=colour, lw=4))
                _labels.append(column_id)
        if not self.subplots:
            plt.legend(fontsize = legendfontsize, \
                handles=_handles, \
                labels=_labels, \
                frameon=True\
                )
        plt.xlabel("Genome coordinate", fontsize=axeslabelfontsize)
        fn = "%s.%s.png" % (self.parameterObj.outprefix, self.name) 
        plt.tight_layout()
        fig.savefig(fn, format="png")
        plt.close(fig)
        return fn


class PlotCovObj():
    def __init__(self, parameterObj, name, y_label, data):
        self.parameterObj = parameterObj
        self.name = name
        self.y_label = y_label
        self.genome_length = None
        self.length_by_pair_count = []
        self.block_length_by_pair_count = []
        self.final_length_by_pair_count = []
        self.length_by_sample_count = []
        self.block_length_by_sample_count = []
        self.final_length_by_sample_count = []
        self.populate(data)
        self.cm = plt.get_cmap('Set2')

    def populate(self, data):
        self.genome_length = sum([y for x, y in data.length_by_sample_count.items()])
        for pair_count in range(self.parameterObj.pairs_count, 0, -1):
            try:
                self.length_by_pair_count.append(self.length_by_pair_count[-1] + data.length_by_pair_count.get(pair_count, 0.0) / self.genome_length)
            except IndexError:
                self.length_by_pair_count.append(data.length_by_pair_count.get(pair_count, 0.0) / self.genome_length)
            try:
                self.block_length_by_pair_count.append(self.block_length_by_pair_count[-1] + data.block_length_by_pair_count.get(pair_count, 0.0) / self.genome_length)
            except IndexError:
                self.block_length_by_pair_count.append(data.block_length_by_pair_count.get(pair_count, 0.0) / self.genome_length)
            try:
                self.final_length_by_pair_count.append(self.final_length_by_pair_count[-1] + data.final_length_by_pair_count.get(pair_count, 0.0) / self.genome_length)
            except IndexError:
                self.final_length_by_pair_count.append(data.block_length_by_pair_count.get(pair_count, 0.0) / self.genome_length)
        for sample_count in range(self.parameterObj.samples_count, 0, -1):
            try:
                self.length_by_sample_count.append(self.length_by_sample_count[-1] + data.length_by_sample_count.get(sample_count, 0.0) / self.genome_length)
            except IndexError:
                self.length_by_sample_count.append(data.length_by_sample_count.get(sample_count, 0.0) / self.genome_length)
            try:
                self.block_length_by_sample_count.append(self.block_length_by_sample_count[-1] + data.block_length_by_sample_count.get(sample_count, 0.0) / self.genome_length)
            except IndexError:
                self.block_length_by_sample_count.append(data.block_length_by_sample_count.get(sample_count, 0.0) / self.genome_length)
            try:
                self.final_length_by_sample_count.append(self.final_length_by_sample_count[-1] + data.final_length_by_sample_count.get(sample_count, 0.0) / self.genome_length) 
            except IndexError:
                self.final_length_by_sample_count.append(data.final_length_by_sample_count.get(sample_count, 0.0) / self.genome_length) 
            
    def plot(self, grouping, xlabel, title):
        categories, y1, y2, y3 = None, None, None, None
        if grouping == "samples":
            categories = self.parameterObj.samples_count
            y1 = self.length_by_sample_count
            y2 = self.block_length_by_sample_count
            y3 = self.final_length_by_sample_count
        else:
            categories = self.parameterObj.pairs_count
            y1 = self.length_by_pair_count
            y2 = self.block_length_by_pair_count
            y3 = self.final_length_by_pair_count
        fig = plt.figure(figsize=(16,6), dpi=400, frameon=False)
        ax = fig.add_subplot(111)
        #x = [int(_x + 1) for _x in range(categories)]
        x = list(range(0, categories, 1))
        colours = [self.cm(idx / 3) for idx in range(3)]
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))
        ax.set_ylim((-0.05, 1.05))
        ax.set_xlim((-1, categories))
        #ax.set_xscale("log")
        plt.title(title)
        plt.ylabel(self.y_label, fontsize=axeslabelfontsize)
        plt.xlabel(xlabel, fontsize=axeslabelfontsize)
        plt.gca().spines["top"].set_visible(False)  
        plt.gca().spines["right"].set_visible(False)
        plt.tight_layout()
        ax.axhline(y=1.0, xmin=0, xmax=1, color='lightgrey', linestyle='--', linewidth=2)
        ax.plot(x, y1, '-ok', label="BED Intervals (>=1 b)", color=colours[0], markersize=6, linewidth=6, markerfacecolor='black')
        plt.legend(fontsize = legendfontsize, loc='lower right', frameon=True)
        fn = "%s.%s.%s.1.png" % (self.parameterObj.outprefix, self.name, grouping) 
        fig.savefig(fn, format="png")
        ax.plot(x, y2, '-ok', label="BED Intervals (>=64 b)", color=colours[1], markersize=6, linewidth=6, markerfacecolor='black')
        plt.legend(fontsize = legendfontsize, loc='lower right', frameon=True)
        fn = "%s.%s.%s.2.png" % (self.parameterObj.outprefix, self.name, grouping) 
        fig.savefig(fn, format="png")
        ax.plot(x, y3, '-ok', label="Blocktools (64 b)", color=colours[2], markersize=6, linewidth=6, markerfacecolor='black')
        plt.legend(fontsize = legendfontsize, loc='lower right', frameon=True)
        fn = "%s.%s.%s.3.png" % (self.parameterObj.outprefix, self.name, grouping) 
        fig.savefig(fn, format="png")
        plt.close(fig)
        return 1

def plot_window_coverage_tsv(parameterObj, sequence_OrdDict):
    df_window_coverage_tsv = pandas.read_csv( \
        parameterObj.window_coverage_tsv_f, \
        sep="\t" \
        )
    df_window_coverage_tsv = df_window_coverage_tsv.dropna()
    df_mean_block_density = df_window_coverage_tsv[df_window_coverage_tsv.columns[3:4].tolist()]
    plotHistObj = PlotHistObj(parameterObj, "mean_window_density", "Value", "Count", df_mean_block_density)
    fn_mean_block_density_hist = plotHistObj.plot()

    sample_cov_df_metrics = df_window_coverage_tsv[df_window_coverage_tsv.columns[0:1].tolist() + df_window_coverage_tsv.columns[-(len(df_window_coverage_tsv.columns)-6):].tolist()]
    plotGenomeObj = PlotGenomeObj(parameterObj, "sample_coverage", sample_cov_df_metrics, sequence_OrdDict, subplots=False, by_population=True)
    fn_sample_cov_genome = plotGenomeObj.plot()
    return (fn_mean_block_density_hist, fn_sample_cov_genome)

def plot_pairs_comparison(parameterObj_1, parameterObj_2):
    df_variant_pairs_1 = pandas.read_csv( \
        parameterObj_1.variant_pairs_tsv_f, \
        sep="\t" \
        )
    df_variant_pairs_1 = df_variant_pairs_1.dropna()
    
    df_variant_pairs_2 = pandas.read_csv( \
        parameterObj_2.variant_pairs_tsv_f, \
        sep="\t" \
        )
    df_variant_pairs_2 = df_variant_pairs_2.dropna()
    
    df_variant_pairs_1_fst = df_variant_pairs_1[['pair_idx', 'f_st']]
    df_variant_pairs_1_fst.insert(0, 'type', os.path.basename(parameterObj_1.outprefix))
    df_variant_pairs_2_fst = df_variant_pairs_2[['pair_idx', 'f_st']]
    df_variant_pairs_2_fst.insert(0, 'type', os.path.basename(parameterObj_2.outprefix))
    df_variant_pairs_fst = pandas.concat([df_variant_pairs_1_fst, df_variant_pairs_2_fst], axis=0)
    plotSwarmObj = PlotSwarmObj('pair_fst_comparison', 'F_st %s' % parameterObj_1.outprefix, 'type', 'F_st %s' % parameterObj_1.outprefix, 'f_st', df_variant_pairs_fst)
    fn = plotSwarmObj.plot(alpha=0.8, points='Pairs')

    df_variant_pairs_1_dxy = df_variant_pairs_1[['pair_idx', 'd_xy']]
    df_variant_pairs_1_dxy.insert(0, 'type', os.path.basename(parameterObj_1.outprefix))
    df_variant_pairs_2_dxy = df_variant_pairs_2[['pair_idx', 'd_xy']]
    df_variant_pairs_2_dxy.insert(0, 'type', os.path.basename(parameterObj_2.outprefix))
    df_variant_pairs_dxy = pandas.concat([df_variant_pairs_1_dxy, df_variant_pairs_2_dxy], axis=0)
    plotSwarmObj = PlotSwarmObj('pair_dxy_comparison', 'D_xy %s' % parameterObj_1.outprefix, 'type', 'D_xy %s' % parameterObj_1.outprefix, 'd_xy', df_variant_pairs_dxy)
    fn = plotSwarmObj.plot(alpha=0.8, points='Pairs')
    return fn 

def plot_variant_pairs_tsv(parameterObj, sequence_OrdDict):
    df_variant_pairs = pandas.read_csv( \
        parameterObj.variant_pairs_tsv_f, \
        sep="\t" \
        )
    df_variant_pairs = df_variant_pairs.dropna()

    df_variant_pairs_piA = df_variant_pairs[df_variant_pairs.columns[0:1].tolist() + df_variant_pairs.columns[9:10].tolist()]
    plotHeatmapObj = PlotHeatmapObj(parameterObj, 'pi_A', 'piA (as proportion of highest value)', df_variant_pairs_piA)
    piA_fn = plotHeatmapObj.plot()

    df_variant_pairs_piB = df_variant_pairs[df_variant_pairs.columns[0:1].tolist() + df_variant_pairs.columns[10:11].tolist()]
    plotHeatmapObj = PlotHeatmapObj(parameterObj, 'pi_B', 'piB (as proportion of highest value)', df_variant_pairs_piB)
    piB_fn = plotHeatmapObj.plot()

    df_variant_pairs_dxy = df_variant_pairs[df_variant_pairs.columns[0:1].tolist() + df_variant_pairs.columns[11:12].tolist()]
    plotHeatmapObj = PlotHeatmapObj(parameterObj, 'd_xy', 'Dxy (as proportion of highest value)', df_variant_pairs_dxy)
    dxy_fn = plotHeatmapObj.plot()

    df_variant_pairs_fst = df_variant_pairs[df_variant_pairs.columns[0:1].tolist() + df_variant_pairs.columns[12:13].tolist()]
    plotHeatmapObj = PlotHeatmapObj(parameterObj, 'f_st', 'Fst (as proportion of highest value)', df_variant_pairs_fst)
    fst_fn = plotHeatmapObj.plot()

    return (piA_fn, piB_fn, dxy_fn, fst_fn)

def plot_coverage(parameterObj, sequence_OrdDict):
    block_span_df = pandas.read_csv( \
        parameterObj.block_pairs_f, \
        sep="\t" \
        ) 
    plotHeatmapObj = PlotHeatmapObj(parameterObj, 'bases', 'Blocked span (as proportion of highest value)', block_span_df)
    fn_block_span_tsv = plotHeatmapObj.plot()
    block_summary_df = pandas.read_csv( \
        parameterObj.block_summary_f, \
        sep="\t" \
        )
    final_length_by_sample_count = dict(block_summary_df.groupby('count_samples')['length'].sum())
    final_length_by_pair_count = dict(block_summary_df.groupby('count_pairs')['length'].sum())
    coverage_dict = parse_yaml(parameterObj.bed_coverage_f)
    coverageObj = CoverageObj()
    coverageObj.set_from_dict(coverage_dict)
    coverageObj.add_final_lengths(final_length_by_sample_count, final_length_by_pair_count)
    plotCovObj = PlotCovObj(parameterObj, 'coverage', 'Proportion of genome', coverageObj)
    plotCovObj.plot('samples', xlabel='Number of missing Samples', title="Proportion of genome visible when allowing for missing Samples")
    plotCovObj.plot('pairs', xlabel='Number of missing Pairs', title="Proportion of genome visible when allowing for missing Pairs")
    return fn_block_span_tsv

        
def plot_window_variant_tsv(parameterObj, sequence_OrdDict):
    df_window_variant_tsv = pandas.read_csv( \
        parameterObj.window_variant_tsv_f, \
        sep="\t" \
        )
    df_window_variant_tsv = df_window_variant_tsv.dropna()

    dxy_fst_df = df_window_variant_tsv[['window_id', 'd_xy', 'f_st']]
    plotGenomeObj = PlotGenomeObj(parameterObj, "dxy_fst", dxy_fst_df, sequence_OrdDict, subplots=True, by_population=False)
    dxy_fst_fn = plotGenomeObj.plot()

    piA_piB_df = df_window_variant_tsv[['window_id', 'pi_A', 'pi_B']]
    plotGenomeObj = PlotGenomeObj(parameterObj, "piA_piB", piA_piB_df, sequence_OrdDict, subplots=True, by_population=False)
    piA_piB_fn = plotGenomeObj.plot()

    profile_df = df_window_fst_vs_missing = df_window_variant_tsv[['window_id', 'missing', 'multiallelic']]
    plotGenomeObj = PlotGenomeObj(parameterObj, "missing_multiallelic", profile_df, sequence_OrdDict, subplots=True, by_population=False)
    tuple_fn = plotGenomeObj.plot()
    
    profile_df = df_window_fst_vs_missing = df_window_variant_tsv[['window_id', 'hetA', 'hetB', 'hetAB', 'fixed']]
    plotGenomeObj = PlotGenomeObj(parameterObj, "tuple", profile_df, sequence_OrdDict, subplots=False, by_population=False)
    tuple_fn = plotGenomeObj.plot()

    df_window_fst_vs_missing = df_window_variant_tsv[['missing', 'f_st']]
    plotScatterObj = PlotScatterObj(parameterObj, 'window_fst_vs_missing', 'Proportion of bases with "missing" genotypes', 'missing', 'F_st across all blocks', 'f_st', df_window_fst_vs_missing)
    plotScatterObj.plot(alpha=0.1, points='Window of %s blocks' % parameterObj.window_size)

    df_window_fst_vs_multiallelic = df_window_variant_tsv[['multiallelic', 'f_st']]
    plotScatterObj = PlotScatterObj(parameterObj, 'window_fst_vs_multiallelic', 'Proportion of bases with "multiallelic" genotypes', 'multiallelic', 'F_st across all blocks', 'f_st', df_window_fst_vs_multiallelic)
    plotScatterObj.plot(alpha=0.1, points='Window of %s blocks' % parameterObj.window_size)

    df_window_piA_vs_missing = df_window_variant_tsv[['missing', 'pi_A']]
    plotScatterObj = PlotScatterObj(parameterObj, 'window_piA_vs_missing', 'Proportion of bases with "missing" genotypes', 'missing', 'Pi_A across all blocks', 'pi_A', df_window_piA_vs_missing)
    plotScatterObj.plot(alpha=0.1, points='Window of %s blocks' % parameterObj.window_size)

    df_window_piA_vs_multiallelic = df_window_variant_tsv[['multiallelic', 'pi_A']]
    plotScatterObj = PlotScatterObj(parameterObj, 'window_piA_vs_multiallelic', 'Proportion of bases with "multiallelic" genotypes', 'multiallelic', 'Pi_A across all blocks', 'pi_A', df_window_piA_vs_multiallelic)
    plotScatterObj.plot(alpha=0.1, points='Window of %s blocks' % parameterObj.window_size)

    df_window_dxy_vs_missing = df_window_variant_tsv[['missing', 'd_xy']]
    plotScatterObj = PlotScatterObj(parameterObj, 'window_dxy_vs_missing', 'Proportion of bases with "missing" genotypes', 'missing', 'D_xy across all blocks', 'd_xy', df_window_dxy_vs_missing)
    plotScatterObj.plot(alpha=0.1, points='Window of %s blocks' % parameterObj.window_size)

    df_window_dxy_vs_multiallelic = df_window_variant_tsv[['multiallelic', 'd_xy']]
    plotScatterObj = PlotScatterObj(parameterObj, 'window_dxy_vs_multiallelic', 'Proportion of bases with "multiallelic" genotypes', 'multiallelic', 'D_xy across all blocks', 'd_xy', df_window_dxy_vs_multiallelic)
    plotScatterObj.plot(alpha=0.1, points='Window of %s blocks' % parameterObj.window_size)

    df_window_fixed_vs_missing = df_window_variant_tsv[['missing', 'fixed']]
    plotScatterObj = PlotScatterObj(parameterObj, 'window_fixed_vs_missing', 'Proportion of bases with "missing" genotypes', 'missing', 'Fixed mutations across all blocks', 'fixed', df_window_fixed_vs_missing)
    plotScatterObj.plot(alpha=0.1, points='Window of %s blocks' % parameterObj.window_size)

    df_window_fixed_vs_multiallelic = df_window_variant_tsv[['multiallelic', 'fixed']]
    plotScatterObj = PlotScatterObj(parameterObj, 'window_fixed_vs_multiallelic', 'Proportion of bases with "multiallelic" genotypes', 'multiallelic', 'Fixed mutations across all blocks', 'fixed', df_window_fixed_vs_multiallelic)
    plotScatterObj.plot(alpha=0.1, points='Window of %s blocks' % parameterObj.window_size)

    df_window_hetAB_vs_missing = df_window_variant_tsv[['missing', 'hetAB']]
    plotScatterObj = PlotScatterObj(parameterObj, 'window_hetAB_vs_missing', 'Proportion of bases with "missing" genotypes', 'missing', 'HetAB mutations across all blocks', 'hetAB', df_window_hetAB_vs_missing)
    plotScatterObj.plot(alpha=0.1, points='Window of %s blocks' % parameterObj.window_size)

    df_window_hetAB_vs_multiallelic = df_window_variant_tsv[['multiallelic', 'hetAB']]
    plotScatterObj = PlotScatterObj(parameterObj, 'window_hetAB_vs_multiallelic', 'Proportion of bases with "multiallelic" genotypes', 'multiallelic', 'HetAB mutations across all blocks', 'hetAB', df_window_hetAB_vs_multiallelic)
    plotScatterObj.plot(alpha=0.1, points='Window of %s blocks' % parameterObj.window_size)

    return (dxy_fst_fn, piA_piB_fn, tuple_fn)


#########################

CONFIG_BY_ZYGOSITY = {
    'MIS' : {
        'HOM' : 'missing',
        'HET' : 'missing',
        'MIS' : 'missing'
    },
    'HOM' : {
        'HOM' : 'fixed', 
        'HET' : 'hetB',
        'MIS' : 'missing'
    },
    'HET' : {
        'HOM' : 'hetA',
        'HET' : 'hetAB',
        'MIS' : 'missing'
    }
}

def pairs_to_samples(pair_idxs, parameterObj):
    return set(chain.from_iterable([parameterObj.sample_idxs_by_pair_idx[pair_idx] for pair_idx in pair_idxs]))

def get_pair_id(pair_id, idx, parameterObj):
    return parameterObj.pair_ids_by_pair_idx[pair_id][idx]

def write_yaml(data, yaml_f):
    with open(yaml_f, 'w') as yaml_fh:
        # yaml.safe_dump(data, yaml_fh) # this does not work with fancy things such as frozenset-keys in dicts ...
        yaml.safe_dump(data, yaml_fh)

def memory_usage_psutil():
    # return the memory usage in MB
    process = psutil.Process(os.getpid())
    mem = process.memory_info()[0] / float(2 ** 20)
    return mem

@contextmanager
def poolcontext(*args, **kwargs):
    pool = multiprocessing.Pool(*args, **kwargs)
    yield pool
    pool.terminate()

def check_file(infile):
    if infile:
        if not os.path.exists(os.path.expanduser(infile)):
            return None
        return os.path.expanduser(infile)
    else:
        return None

def parse_parameters(args):
    _args = args
    if args['--yaml']:
        args = parse_yaml(args['--yaml'])['args']
        for k, v in _args.items():
            if v:
                args[k] = v
    try:
        parameterObj = ParameterObj(args)
    except PopsCountException:
        sys.exit("[X] Only works on two populations ...")
    return parameterObj

def parse_yaml(yaml_f):
    check_file(yaml_f)
    with open(yaml_f, 'r') as stream:
        # config = yaml.safe_load(stream) # this does not work with fancy things such as frozenset-keys in dicts ...
        config = yaml.safe_load(stream) # this does not work with fancy things such as frozenset-keys in dicts ...
    return config

def parse_genome_f(parameterObj):
    df = pandas.read_csv(parameterObj.genome_f, sep="\t", usecols=[0, 1], names=['chrom', 'length'], header=None)
    return OrderedDict({chrom: int(length) for chrom, length in df.values.tolist()})

def parse_coordinate_transformation_f(parameterObj):
    mapping_df = pandas.read_csv(\
        parameterObj.coordinate_f, \
        sep="\t", \
        usecols=[0, 1, 2, 3, 4, 5, 6], \
        names=['chrom', 'chrom_start', 'chrom_end', 'seq_id', 'seq_start', 'seq_end', 'seq_orientation'], \
        skiprows=2, \
        header=None, \
        dtype={ \
            'chrom': 'category', \
            'chrom_start': numpy.int, \
            'chrom_end': numpy.int, \
            'seq_id': 'category', \
            'seq_start': numpy.int, \
            'seq_end': numpy.int, \
            'seq_orientation': 'category' \
            })
    # convert coordinates to 0-based (BED)
    mapping_df['chrom_start'] = mapping_df['chrom_start'] - 1
    mapping_df['seq_start'] = mapping_df['seq_start'] - 1
    coordinateTransformObj = CoordinateTransformObj()
    for chrom, chrom_start, chrom_end, seq_id, seq_start, seq_end, seq_orientation in tqdm(mapping_df.values.tolist(), total=len(mapping_df.index), desc="[%] ", ncols=200):
        coordinateTransformObj.add_tuple(seq_id, int(seq_start), int(seq_end), seq_orientation, chrom, int(chrom_start), int(chrom_end))
    return coordinateTransformObj

def parse_multibed_f(parameterObj, sequence_OrdDict):
    df = pandas.read_csv( \
        parameterObj.multibed_f, \
        sep="\t", \
        usecols=[0, 1, 2, 4], \
        names=['chrom', 'start', 'end', 'samples'], \
        skiprows=1, \
        header=None, \
        dtype={ \
            'chrom': 'category', \
            'start': numpy.int, \
            'end': numpy.int, \
            'samples': 'category'\
        })
    # filter rows based on chrom
    df = df[df['chrom'].isin(sequence_OrdDict)]
    # compute length column
    df['length'] =  df['end'] - df['start']
    # filter intervals shorter than MIN_INTERVAL_LEN
    # df = df[df['length'] >= parameterObj.min_interval_len]
    # samples to frozenset
    df['samples_idxs'] = df['samples'].apply(\
        generate_sample_idxs, \
        sample_idx_by_sample_id=parameterObj.sample_idx_by_sample_id\
        )
    # generate pairs in dataframe based on combo-intersect
    df['pair_idxs'] = df['samples'].apply(\
        generate_pair_idxs, \
        pops_count=parameterObj.pops_count, \
        pair_idx_by_pair_ids=parameterObj.pair_idx_by_pair_ids\
        )
    # Drop intervals that don't affect pairs
    df = df.dropna()
    
    # compute distance to next interval
    df['distance'] = numpy.where((df['chrom'] == df['chrom'].shift(-1)), df['start'].shift(-1) - df['end'], numpy.nan)
    print(df)
    #spanObj = SpanObj('raw')
    idx = 0
    regionBatchObjs = deque()
    regionBatchObj = RegionBatchObj(idx)
    coverageObj = CoverageObj()
    coverageObj.set_genome_length(sequence_OrdDict) # better: record BED intervals for each contig
    for chrom, start, end, samples, length, sample_idxs, pair_idxs, distance in tqdm(df.values.tolist(), total=len(df.index), desc="[%] ", ncols=200):
        #if not numpy.isnan(pair_idxs):
        pair_count = len(pair_idxs)
        coverageObj.add_pair_region(pair_count, length)
        if length >= parameterObj.block_length:
            coverageObj.add_pair_region_block_size(pair_count, int(length / parameterObj.block_length) * parameterObj.block_length)
        if length >= parameterObj.min_interval_len:
            bedObj = BedObj(chrom, start, end, pair_idxs, length) 
            if not regionBatchObj.contig_id:
                regionBatchObj.add_bedObj_to_batch(bedObj)
            else:
                if chrom == regionBatchObj.contig_id and not numpy.isnan(distance) and int(distance) <= parameterObj.max_interval_distance:
                    regionBatchObj.add_bedObj_to_batch(bedObj)
                else:
                    regionBatchObjs.append(regionBatchObj)
                    idx += 1
                    regionBatchObj = RegionBatchObj(idx)
                    regionBatchObj.add_bedObj_to_batch(bedObj)
            #if numpy.isnan(distance) or int(distance) > parameterObj.max_interval_distance:
        
        sample_count = len(sample_idxs)
        coverageObj.add_sample_region(sample_count, length)
        if length >= parameterObj.block_length:
            coverageObj.add_sample_region_block_size(sample_count, int(length / parameterObj.block_length) * parameterObj.block_length)
    if regionBatchObj.contig_id:
        regionBatchObjs.append(regionBatchObj)
    #for x in regionBatchObjs:
    #    print(x)
        #for y in x.bedObjs:
        #    print(y)

    # awk '{if($3-$2>=65){a[$4]+=$3-$2;}}END{for (i in a)print i", "a[i];}' data/input/tiny/hmel.autosomes.tiny.multiinter.samples_as_string.bed
    return regionBatchObjs, coverageObj

class CoverageObj():
    
    def __init__(self):
        self.length_by_pair_count = {}
        self.block_length_by_pair_count = {}
        self.length_by_sample_count = {}
        self.block_length_by_sample_count = {}
        self.genome_length = None
        self.final_length_by_sample_count = {}
        self.final_length_by_pair_count = {}

    def set_from_dict(self, coverage_dict):
        self.length_by_pair_count = coverage_dict["length_by_pair_count"]
        self.block_length_by_pair_count = coverage_dict["block_length_by_pair_count"]
        self.length_by_sample_count = coverage_dict["length_by_sample_count"]
        self.block_length_by_sample_count = coverage_dict["block_length_by_sample_count"]
        self.genome_length = coverage_dict["genome_length"]

    def add_final_lengths(self, final_length_by_sample_count, final_length_by_pair_count):
        self.final_length_by_sample_count = final_length_by_sample_count
        self.final_length_by_pair_count = final_length_by_pair_count

    def set_genome_length(self, sequence_OrdDict):
        self.genome_length = sum(sequence_OrdDict.values())

    def add_sample_region_block_size(self, sample_count, length):
        self.block_length_by_sample_count[sample_count] = self.block_length_by_sample_count.get(sample_count, 0) + length

    def add_pair_region_block_size(self, pair_count, length):
        self.block_length_by_pair_count[pair_count] = self.block_length_by_pair_count.get(pair_count, 0) + length

    def add_sample_region(self, sample_count, length):
        self.length_by_sample_count[sample_count] = self.length_by_sample_count.get(sample_count, 0) + length

    def add_pair_region(self, pair_count, length):
        self.length_by_pair_count[pair_count] = self.length_by_pair_count.get(pair_count, 0) + length

def load_profileObjs(parameterObj, blockDataObj):
    profile_df = pandas.read_csv(\
        parameterObj.variant_blocks_tsv_f, \
        sep="\t", \
        usecols=[0, 1, 2, 3, 4, 5, 6, 7], \
        names=['block_id', 'pair_idx', 'fixed', 'hetA', 'hetAB', 'hetB', 'multiallelic', 'missing'], \
        skiprows=1, \
        header=None, \
        dtype={ \
            'block_id': 'category', \
            'pair_idx': numpy.int, \
            'fixed': numpy.int, \
            'hetA': numpy.int, \
            'hetAB': numpy.int, \
            'hetB': numpy.int, \
            'multiallelic': numpy.int, \
            'missing': numpy.int \
            } \
        )
    for block_id, pair_idx, fixed, hetA, hetAB, hetB, multiallelic, missing in tqdm(profile_df.values.tolist(), total=len(profile_df.index), desc="[%] ", ncols=200):
        profileObj = ProfileObj((fixed, hetA, hetAB, hetB, missing, multiallelic))
        try:
            blockObj = blockDataObj.blockObjs_by_block_id[block_id]
            blockObj.profileObj_by_pair_idx[pair_idx] = profileObj
        except KeyError:
            pass
    return blockDataObj

def load_blockDataObj(parameterObj):
    bed_f = parameterObj.block_bed_f
    if parameterObj.new_bed_f:
        bed_f = parameterObj.new_bed_f
    bed_df = pandas.read_csv(\
        bed_f, \
        sep="\t", \
        usecols=[0, 1, 2, 3], \
        names=['chrom', 'start', 'end', 'block_id'], \
        skiprows=1, \
        header=None, \
        dtype={ \
            'chrom': 'category', \
            'start': numpy.int, \
            'end': numpy.int, \
            'block_id': 'category', \
            } \
        )
    tsv_df = pandas.read_csv(\
        parameterObj.block_summary_f, \
        sep="\t", \
        usecols=[0, 1, 2, 5, 6], \
        names=['block_id', 'length', 'span', 'sample_idxs', 'pair_idxs'], \
        skiprows=1, \
        header=None, \
        dtype={ \
            'block_id': 'category', \
            'length': numpy.int, \
            'span': numpy.int, \
            'sample_idxs' : 'category', \
            'pair_idxs': 'category'
            } \
        )
    bed_tuples_by_block_id = defaultdict(list)
    for chrom, start, end, block_id in bed_df.values.tolist():
        bed_tuples_by_block_id[block_id].append((chrom, start, end))
    blockObjs = []
    for block_id, length, span, sample_idxs, pair_idxs in tqdm(tsv_df.values.tolist(), total=len(tsv_df.index), desc="[%] ", ncols=200):
        if block_id in bed_tuples_by_block_id: # only those in BED file are instanciated !!!
            pair_idxs = [int(pair_idx) for pair_idx in pair_idxs.split(",")]
            blockObj = BlockObj(block_id, parameterObj.block_length)
            #print("#", blockObj)
            for bed_tuple in bed_tuples_by_block_id[block_id]:
                #bedObj = BedObj(bed_tuple[0], bed_tuple[1], bed_tuple[2], pair_idxs, length)
                bedObj = BedObj(bed_tuple[0], bed_tuple[1], bed_tuple[2], pair_idxs, bed_tuple[2] - bed_tuple[1])
                #print("=>", bedObj)
                #_bed = blockObj.add_bedObj(bedObj, parameterObj)
                blockObj.add_bedObj(bedObj, parameterObj)
                #print("<=", _bed)
                #print("==", blockObj)
            #print("#", blockObj.bed_tuples)
            blockObj.profileObj_by_pair_idx = {pair_idx: ProfileObj((0, 0, 0, 0, 0, 0)) for pair_idx in pair_idxs}
            blockObjs.append(blockObj)
    blockDataObj = BlockDataObj(parameterObj) 
    blockDataObj.add_blockObjs(blockObjs)
    return blockDataObj

def generate_pair_idxs(sample_string, **kwargs):
    # generates pairs of interval and returns them as set of idx
    pair_idxs = frozenset(filter(lambda x: x >= 0, [kwargs['pair_idx_by_pair_ids'].get(frozenset(x), -1) for x in combinations(sample_string.split(","), kwargs['pops_count'])])) 
    if pair_idxs:
        return pair_idxs
    return numpy.nan

def generate_sample_idxs(sample_string, **kwargs):
    sample_idxs = frozenset([kwargs['sample_idx_by_sample_id'][sample_id] for sample_id in sample_string.split(",") if sample_id in kwargs['sample_idx_by_sample_id']])
    if sample_idxs:
        return sample_idxs
    else:
        return numpy.nan

def make_blocks(parameterObj, regionBatchObjs):
    if parameterObj.algorithm == "A":
        if parameterObj.min_interval_len >= parameterObj.block_length:
            algorithm = block_algorithm_a 
        else:
            sys.exit('[X] Algorithm A condition not met: MIN_INTERVAL_LEN (%s) >= block_length (%s)' % (parameterObj.min_interval_len, parameterObj.block_length))    
    elif parameterObj.algorithm == "B": 
        if parameterObj.min_interval_len < parameterObj.block_length:
            algorithm = block_algorithm_b
        else:
            sys.exit('[X] Algorithm B condition not met: MIN_INTERVAL_LEN (%s) < block_length (%s)' % (parameterObj.min_interval_len, parameterObj.block_length))    
    elif parameterObj.algorithm == "D":
        algorithm = block_algorithm_d
    else:
        sys.exit('[X] You broke my software ...')
    _temp = []
    print("[+] Analysing %s regionBatchObjs on %s threads ..." % (len(regionBatchObjs), parameterObj.threads))
    if parameterObj.threads < 2:
        with tqdm(total=len(regionBatchObjs), desc="[%] ", ncols=200, unit_scale=True) as pbar:
            for regionBatchObj in regionBatchObjs:
                if regionBatchObj.length() >= parameterObj.block_length:
                    blockObjs = algorithm((regionBatchObj, parameterObj))
                    #print(type(blockObjs))
                    _temp.append(blockObjs)
                    #for blockObj in blockObjs:
                    #    print("[DONE]", blockObj)
                pbar.update()
    else:
        # if multiple threads then arguments have to be passed to algorithm
        params = [(regionBatchObj, parameterObj) for regionBatchObj in regionBatchObjs if regionBatchObj.length() >= parameterObj.block_length]
        with poolcontext(processes=parameterObj.threads) as pool:
            with tqdm(total=len(params), desc="[%] ", ncols=200, unit_scale=True) as pbar:
                for blockObjs in pool.imap_unordered(algorithm, params):
                    #print(type(blockObjs))
                    _temp.append(blockObjs)
                    #for blockObj in blockObjs:
                        #print("[DONE]", blockObj)
                    pbar.update()
    #print("[+] Generated %s blocks ..." % len(_temp))
    flat_list = [_list for sublist in _temp for _list in sublist]
    blockDataObj = BlockDataObj(parameterObj) 
    blockDataObj.add_blockObjs(flat_list)
    return blockDataObj

def block_algorithm_a(params):
    regionBatchObj, parameterObj = params
    block_idx = 0
    blockObjs = deque()
    #print("[R]", regionBatchObj)
    while 1:
        try:
            bedObj = regionBatchObj.bedObjs.popleft()
            fraction = int(floor(bedObj.length / parameterObj.block_length))
            #print("[BED]", fraction, bedObj)
            _start = bedObj.start
            _end = _start + parameterObj.block_length
            for j in range(0, fraction):
                block_id = "%s.i%s.b%s" % (regionBatchObj.contig_id, regionBatchObj.idx, block_idx)
                blockObj = BlockObj(block_id, parameterObj.block_length)
                bedObj = BedObj( \
                    bedObj.chrom, \
                    _start, \
                    _end, \
                    bedObj.pair_idxs, \
                    parameterObj.block_length)
                blockObj.add_bedObj(bedObj, parameterObj)
                #print(blockObj)
                blockObjs.append(blockObj)
                block_idx += 1
                _start += parameterObj.block_length
                _end += parameterObj.block_length
        except IndexError:
            break
    return blockObjs

def block_algorithm_b(params):
    '''
    - able to jump gaps in intervals
    - NOT able to jump intervals
    '''

    regionBatchObj, parameterObj = params
    block_idx = 0
    blockObjs = deque()
    block_id = "%s.i%s.b%s" % (regionBatchObj.contig_id, regionBatchObj.idx, block_idx)
    #print("\n[NEW]", blockObj)
    blockObj = BlockObj(block_id, parameterObj.block_length)
    while 1:
        try:
            bedObj = regionBatchObj.bedObjs.popleft()
            #print("->", bedObj)
            bedObj_score = (len(bedObj.pair_idxs) / parameterObj.pairs_count) * (min(bedObj.length, blockObj.needed) / parameterObj.block_length)
            #bedObj_score = (len(bedObj.pair_idxs) / parameterObj.pairs_count) * parameterObj.block_length
            remainder_bedObj = blockObj.add_bedObj(bedObj, parameterObj)
            if blockObj.score >= bedObj_score:
                if not blockObj.needed:
                    #print("[Done]", blockObj)
                    blockObjs.append(blockObj)
                    #print("\n[NEW]", blockObj)
                    block_idx += 1
                    block_id = "%s.i%s.b%s" % (regionBatchObj.contig_id, regionBatchObj.idx, block_idx)
                    blockObj = BlockObj(block_id, parameterObj.block_length)
                    if remainder_bedObj:
                        #print("[<-]", remainder_bedObj)
                        regionBatchObj.bedObjs.appendleft(remainder_bedObj)
            else: # span violation
                #print("[<-]", bedObj)
                regionBatchObj.bedObjs.appendleft(bedObj)
                #print("\n[NEW]", blockObj)
                blockObj = BlockObj(block_id, parameterObj.block_length)
        except IndexError:
            break
    return blockObjs

def block_algorithm_d(params):
    '''
    - able to jump gaps in intervals
    - NOT able to jump intervals
    - does each pair_idx separately
    '''
    regionBatchObj, parameterObj = params
    blockObjs = []
    for pair_idx in parameterObj.pair_idxs:
        bedObjs = deque([bedObj for bedObj in copy.deepcopy(regionBatchObj.bedObjs) if pair_idx in bedObj.pair_idxs])
        if bedObjs:
            block_idx = 0
            block_id = "%s.i%s.b%s.p%s" % (regionBatchObj.contig_id, regionBatchObj.idx, block_idx, pair_idx)
            blockObj = BlockObj(block_id, parameterObj.block_length)
            blockObj.pair_idxs = set([pair_idx])
            #  print("\n### PAIR", pair_idx, "beds", len(bedObjs))
            while 1:
                try:
                    bedObj = bedObjs.popleft()
                    #  print("[->]", bedObj)
                    remainder_bedObj = blockObj.add_bedObj(bedObj, parameterObj)
                    if blockObj.score: # no span-violation 
                        ##  print("[BLOCK]", blockObj)
                        if not blockObj.needed:
                            #  print("[Done]", blockObj)
                            blockObjs.append(blockObj)
                            block_idx += 1
                            if remainder_bedObj:
                                bedObjs.appendleft(remainder_bedObj)
                                #  print("[<-]", remainder_bedObj)                    
                            block_id = "%s.i%s.b%s.p%s" % (regionBatchObj.contig_id, regionBatchObj.idx, block_idx, pair_idx)
                            blockObj = BlockObj(block_id, parameterObj.block_length)
                            blockObj.pair_idxs = set([pair_idx])
                            #  print("\n[NEW]", blockObj)
                    else:
                        bedObjs.appendleft(bedObj)
                        # if there is "span violation"
                        blockObj = BlockObj(block_id, parameterObj.block_length)
                        blockObj.pair_idxs = set([pair_idx])
                        #  print("\n[NEW]", blockObj)
                except IndexError:
                    break
    return blockObjs

def parse_vcf(parameterObj, blockDataObj):
    genotypes_by_block_id = {}
    if parameterObj.threads < 2:
        with tqdm(total=len(blockDataObj.blockObjs), desc="[%] ", ncols=200, unit_scale=True) as pbar:
            for blockObj in blockDataObj.blockObjs:
                block_id, genotypes_by_sample_idx = fetch_variants((blockObj, parameterObj))
                genotypes_by_block_id[block_id] = genotypes_by_sample_idx
                pbar.update()
    else:
        # if multiple threads then arguments have to be passed to algorithm
        params = [(blockObj, parameterObj) for blockObj in blockDataObj.blockObjs]
        with poolcontext(processes=parameterObj.threads) as pool:
            with tqdm(total=len(params), desc="[%] ", ncols=200, unit_scale=True) as pbar:
                for block_id, genotypes_by_sample_idx in pool.imap_unordered(fetch_variants, params):    
                    genotypes_by_block_id[block_id] = genotypes_by_sample_idx
                    pbar.update()
    for blockObj in blockDataObj.blockObjs:
        blockObj.genotypes_by_sample_idx = genotypes_by_block_id[blockObj.block_id]
    return blockDataObj

def fetch_variants(param):
    blockObj, parameterObj = param
    genotypes_by_sample_idx = defaultdict(list)
    sample_ids = [parameterObj.sample_id_by_sample_idx[sample_idx] for sample_idx in blockObj.sample_idxs]
    #print("# Fetch:", sample_ids)
    vcf_reader = VCF(parameterObj.vcf_f, samples=sample_ids)
    #print(vcf_reader.samples)
    #vcf_reader.set_samples(sample_ids)
    #print(vcf_reader.samples)
    for contig_id, start, end in blockObj.bed_tuples:
        #for record in vcf_reader('%s:%s-%s' % (blockObj.contig_id, blockObj.start, blockObj.end)): # BUG3
        #print(contig_id, start, end)
        for record in vcf_reader('%s:%s-%s' % (contig_id, start, end)):
            #print(record.start, record.end)
            if record.is_snp:
                #print("# GTs:", record.genotypes)
                #for sample_id, genotype in zip(sample_ids, record.genotypes):
                for sample_id, genotype in zip(vcf_reader.samples, record.genotypes):
                    #print(sample_id, genotype)
                    genotypes_by_sample_idx[parameterObj.sample_idx_by_sample_id[sample_id]].append([genotype[0], genotype[1]])
    return blockObj.block_id, genotypes_by_sample_idx

def analyse_variants(parameterObj, blockDataObj):
    profileObj_by_block_id = {}
    if parameterObj.threads < 2:
        with tqdm(total=len(blockDataObj), desc="[%] ", ncols=200, unit_scale=True) as pbar:
            for blockObj in blockDataObj.blockObjs:
                block_id, profileObj_by_pair_idx = infer_configurations((blockObj, parameterObj))
                profileObj_by_block_id[block_id] = profileObj_by_pair_idx
                pbar.update()
    else:
        # if multiple threads then arguments have to be passed to algorithm
        params = [(blockObj, parameterObj) for blockObj in blockDataObj.blockObjs]
        with poolcontext(processes=parameterObj.threads) as pool:
            with tqdm(total=len(blockDataObj), desc="[%] ", ncols=200, unit_scale=True) as pbar:
                for block_id, profileObj_by_pair_idx in pool.imap_unordered(infer_configurations, params):
                    profileObj_by_block_id[block_id] = profileObj_by_pair_idx
                    pbar.update()
    for blockObj in blockDataObj.blockObjs:
        blockObj.profileObj_by_pair_idx = profileObj_by_block_id[blockObj.block_id]
    return blockDataObj


def infer_configurations(params):
    blockObj, parameterObj = params
    profileObj_by_pair_idx = {}
    #print(blockObj.block_id)
    #print(blockObj.bed_tuples)
    if blockObj.genotypes_by_sample_idx:
        for pair_idx in blockObj.pair_idxs:
            #print("# PAIR:", pair_idx)
            profile_dict = {} 
            sample_idx_A, sample_idx_B = parameterObj.sample_idxs_by_pair_idx[pair_idx]
            #print(parameterObj.sample_id_by_sample_idx[sample_idx_A], parameterObj.sample_id_by_sample_idx[sample_idx_B])
            #print("## GTs: ", len(blockObj.genotypes_by_sample_idx[sample_idx_A]), len(blockObj.genotypes_by_sample_idx[sample_idx_B]))
            for gt_A, gt_B in zip(blockObj.genotypes_by_sample_idx[sample_idx_A], blockObj.genotypes_by_sample_idx[sample_idx_B]):
                genotype_set = set(gt_A + gt_B)
                config = None
                if -1 in genotype_set:
                    config = 'missing'
                else:
                    if len(genotype_set) == 2:
                        config = CONFIG_BY_ZYGOSITY \
                                        [get_zygosity(gt_A)] \
                                        [get_zygosity(gt_B)]
                    elif len(genotype_set) > 2:
                        config = 'multiallelic'
                        #config = 'invariant'
                    else:  # len(genotype_set) < 2:
                        pass
                if config:
                    profile_dict[config] = profile_dict.get(config, 0) + 1
                #print("A =", gt_A, "B =", gt_B, "=>", config)
            #profileObj_by_pair_idx[pair_idx] = ProfileObj(list(profile_dict.values())) # as of Python 3.7, insertion order is maintained in dict.values() (https://docs.python.org/3.7/library/stdtypes.html#dict.values)
            profileObj_by_pair_idx[pair_idx] = ProfileObj((\
                                                profile_dict.get('fixed', 0), \
                                                profile_dict.get('hetA', 0), \
                                                profile_dict.get('hetAB', 0), \
                                                profile_dict.get('hetB', 0), \
                                                profile_dict.get('missing', 0), \
                                                profile_dict.get('multiallelic', 0) \
                                               ))
            #print("### Profile", profileObj_by_pair_idx[pair_idx])
    return blockObj.block_id, profileObj_by_pair_idx

def get_zygosity(gt):
    if gt[0] == gt[1]:
        return 'HOM'
    return 'HET'

def transform_coordinates(parameterObj, blockDataObj, coordinateTransformObj):
    blockObjs = []
    params = [(blockObj, coordinateTransformObj) for blockObj in blockDataObj.blockObjs]
    if parameterObj.threads < 2:
        with tqdm(total=len(params), desc="[%] ", ncols=200, unit_scale=True) as pbar:
            for param in params:
                #print(blockObj, blockObj.void) 
                blockObj = transform_coordinates_blockObj(param)
                #print(blockObj, blockObj.void)
                blockObjs.append(blockObj)
                pbar.update()
    else:
        # if multiple threads then arguments have to be passed to algorithm
        with poolcontext(processes=parameterObj.threads) as pool:
            with tqdm(total=len(params), desc="[%] ", ncols=200, unit_scale=True) as pbar:
                for blockObj in pool.imap_unordered(transform_coordinates_blockObj, params):
                    blockObjs.append(blockObj)
                    pbar.update()
    _blockDataObj = BlockDataObj(parameterObj) 
    _blockDataObj.add_blockObjs(blockObjs)
    return _blockDataObj

def transform_coordinates_blockObj(params):
    blockObj, coordinateTransformObj = params
    #print(">", blockObj.contig_id, blockObj.start, blockObj.end)
    _contig_id, _start, _end = coordinateTransformObj.transform_coordinates(blockObj.contig_id, int(blockObj.start), int(blockObj.end))
    #print("<", _contig_id, _start, _end)
    _bed_tuples = []
    if _contig_id:
        blockObj.contig_id = _contig_id
        blockObj.start = _start
        blockObj.end = _end
        for bed_tuple in blockObj.bed_tuples:
            bed_contig_id, bed_start, bed_end = coordinateTransformObj.transform_coordinates(bed_tuple[0], int(bed_tuple[1]), int(bed_tuple[2]))
            _bed_tuples.append((bed_contig_id, bed_start, bed_end))
        blockObj.bed_tuples = _bed_tuples
    else:
        blockObj.void = True
    return blockObj

def make_windows(parameterObj, blockDataObj):
    windowDataObj = WindowDataObj()
    _lol = []
    # split big chroms further since otherwise

    params = []
    # for each chrom
    for start, end in blockDataObj.blockObj_idxs:
        #print("# Length: %s from %s to %s" % (len(blockDataObj.blockObjs[start:end]), start, end))
        chunk = parameterObj.window_size 
        if len(blockDataObj.blockObjs[start:end]) > chunk:
            init, stop = 0, 0
            _buffer = 0
            for i in range(start, end, chunk):
                init = i + _buffer
                stop = i + 2 * chunk + _buffer
                #print(init, stop)
                if stop > end:
                    stop = end
                params.append((parameterObj, blockDataObj.blockObjs[init:stop]))
                _buffer += parameterObj.window_overlap
        else:
            params.append((parameterObj, blockDataObj.blockObjs[start:end]))
    if parameterObj.threads < 2:
        with tqdm(total=len(params), desc="[%] ", ncols=200, unit_scale=True) as pbar:
            for param in params:
                _windowObjs = window_algorithm(param)
                if _windowObjs:
                    _lol.append(_windowObjs)
                pbar.update()
    else:
        with poolcontext(processes=parameterObj.threads) as pool:
             with tqdm(total=len(params), desc="[%] ", ncols=200, unit_scale=True) as pbar:
                for _windowObjs in pool.imap_unordered(window_algorithm, params):
                    if _windowObjs:
                        _lol.append(_windowObjs)
                    pbar.update()
    flat_list = list(chain.from_iterable(_lol))
    windowDataObj.add_windowObjs(flat_list)
    print(memory_usage_psutil())
    return windowDataObj

def window_algorithm(params):
    _windowObjs = []
    parameterObj, blockObjs = params[0], params[1]
    for i in range(0, len(blockObjs), parameterObj.window_overlap):
        if (len(blockObjs) - i) < parameterObj.window_size:
            break
        else:
            windowObj = WindowObj(blockObjs[0].contig_id, blockObjs[i : i + parameterObj.window_size], parameterObj.block_length)
            #print(i, i + parameterObj.window_size, windowObj)
            _windowObjs.append(windowObj)
    return _windowObjs


def analyse_windows(parameterObj, windowDataObj):
    if not windowDataObj.windowObjs:
        sys.exit("[X] No windows found.")
    params = [(parameterObj, windowObj) for windowObj in windowDataObj.windowObjs.values()]
    if parameterObj.threads < 2:
        with tqdm(total=len(params), desc="[%] ", ncols=200, unit_scale=True) as pbar:
            for windowObj in params:
                windowObj = compute_window_metrics(windowObj)
                windowDataObj.windowObjs[windowObj.window_id] = windowObj
                pbar.update()
    else:
        with poolcontext(processes=parameterObj.threads) as pool:
            with tqdm(total=len(params), desc="[%] ", ncols=200, unit_scale=True) as pbar:
                for windowObj in pool.imap_unordered(compute_window_metrics, params):
                    windowDataObj.windowObjs[windowObj.window_id] = windowObj
                    pbar.update()
    return windowDataObj

def compute_window_metrics(params):
    parameterObj, windowObj = params
    windowObj.profileObj = sum(windowObj.profileObjs) / sum(windowObj.profile_weights)
    windowObj.metrics = calculate_metrics(windowObj.profileObj, parameterObj.block_length)
    return windowObj

def calculate_metrics(profileObj, total_sites):
    # if missing, assume invariant
    pi_A = float("%.8f" % ((profileObj.hetA + profileObj.hetAB) / total_sites))
    pi_B = float("%.8f" % ((profileObj.hetB + profileObj.hetAB) / total_sites))
    d_xy = float("%.8f" % ((((profileObj.hetA + profileObj.hetB + profileObj.hetAB) / 2.0) + profileObj.fixed) / total_sites))
    mean_pi = (pi_A + pi_B) / 2.0
    total_pi = (d_xy + mean_pi) / 2.0 # special case of pairwise Fst
    f_st = numpy.nan
    if (total_pi):
        f_st = float("%.8f" % ((total_pi - mean_pi) / total_pi)) # special case of pairwise Fst
    return {
            "pi_A" : pi_A, \
            "pi_B" : pi_B, \
            "d_xy" : d_xy, \
            "f_st" : f_st \
            }

class CoordinateTransformObj(object):
    def __init__(self):
        self.tuples_by_seq_id = defaultdict(list)
    
    def add_tuple(self, seq_id, seq_start, seq_end, seq_orientation, chrom, chrom_start, chrom_end):
        self.tuples_by_seq_id[seq_id].append((seq_id, seq_start, seq_end, seq_orientation, chrom, chrom_start, chrom_end))

    def transform_coordinates(self, _seq_id, _seq_start, _seq_end):
        new_seq_id, new_start, new_end = None, None, None
        if _seq_id in self.tuples_by_seq_id:
            for _tuple in self.tuples_by_seq_id[_seq_id]:
                if _seq_start >= _tuple[1] and _seq_end <= _tuple[2]:
                    new_seq_id = _tuple[4] 
                    if _tuple[3] == "+":
                        new_start = _tuple[5] + (_seq_start - _tuple[1])
                        new_end = _tuple[5] + (_seq_end - _tuple[1])
                        return (new_seq_id, new_start, new_end)
                    else:
                        new_start = _tuple[5] + (_tuple[2] - _seq_end)
                        new_end = _tuple[5] + (_tuple[2] - _seq_start)
                        return (new_seq_id, new_start, new_end)
        return (new_seq_id, new_start, new_end)

class ParameterObj(object):
    '''
    Class containing all parameters necessary for:
        - blocks
        - variants
        - windows
        - plots
        - picking up after each checkpoint (!!!)
    '''
    def __init__(self, args):
        # input files
        self.genome_f = check_file(args.get('--genome', None))
        self.multibed_f = check_file(args.get('--multibed', None))
        self.vcf_f = check_file(args.get('--vcf', None))
        self.coordinate_f = check_file(args.get('--coordinates', None))
        self.new_bed_f = check_file(args.get('--bed', None))
        # analysis parameters
        self.threads = int(args['--threads']) - 1 if '--threads' in args else 1
        self.algorithm = args['--algorithm']
        self.block_length = int(args['--block_length'])
        self.min_interval_len = int(args['--min_interval_len'])
        self.max_interval_distance = int(args['--max_interval_distance'])
        self.block_span = self.max_interval_distance

        #self.modifier = float(args['--modifier'])
        self.window_size = int(args['--window_size']) if '--window_size' in args else 500
        self.window_overlap = int(args['--overlap']) if '--overlap' in args else 100
        
        # Samples/Pops

        self.sample_ids_by_population = args['--populations']
        print("sample_ids_by_population", self.sample_ids_by_population)
        self.population_by_population_idx = OrderedDict({idx: population_id for idx, population_id in enumerate(natsorted(self.sample_ids_by_population.keys()))})
        print("population_by_population_idx", self.population_by_population_idx)
        self.population_by_sample_id = self.get_population_by_sample_id()
        print("population_by_sample_id", self.population_by_sample_id)
        self.colour_by_population = {population : colour for population, colour in zip(args['--populations'], ['gold', 'purple'])}
        print("colour_by_population", self.colour_by_population)
        self.colour_default = 'purple'
        self.pops_count = len(self.sample_ids_by_population)
        if not self.pops_count == 2:
            raise PopsCountException()
        self.sample_ids = [sample_id for sample_ids in sorted(self.sample_ids_by_population.values()) for sample_id in sample_ids]
        self.samples_count = len(self.sample_ids)
        self.sample_id_by_sample_idx = {idx: sample_id for idx, sample_id in enumerate(self.sample_ids)}
        self.sample_idx_by_sample_id = {sample_id: idx for idx, sample_id in enumerate(self.sample_ids)}
        self.pair_ids = [(x) for x in product(*sorted(self.sample_ids_by_population.values()))]
        print("pair_ids", self.pair_ids)
        self.pairs_count = len(self.pair_ids)
        self.pair_idxs = [pair_idx for pair_idx, pair_id in enumerate(self.pair_ids)]
        print("pair_idxs", self.pair_idxs)
        self.pair_idx_by_pair_ids = {frozenset(pair_id): pair_idx for pair_idx, pair_id in zip(self.pair_idxs, self.pair_ids)}
        print("pair_idx_by_pair_ids", self.pair_idx_by_pair_ids)
        self.pair_ids_by_pair_idx = {pair_idx: pair_id for pair_idx, pair_id in zip(self.pair_idxs, self.pair_ids)}
        print("pair_ids_by_pair_idx", self.pair_ids_by_pair_idx)
        self.sample_idxs_by_pair_idx = {pair_idx: (self.sample_idx_by_sample_id[pair_id[0]], self.sample_idx_by_sample_id[pair_id[1]]) for pair_idx, pair_id in zip(self.pair_idxs, self.pair_ids)}
        print("sample_idxs_by_pair_idx", self.sample_idxs_by_pair_idx)
        # Output files
        self.outprefix = args['--outprefix']
        self.block_bed_f = "%s.block.bed" % (self.outprefix)
        self.block_bed_fix_f = "%s.block.fix.bed" % (self.outprefix)
        self.block_void_bed_f = "%s.block.void.bed" % (self.outprefix)
        self.block_summary_f = "%s.block.summary.tsv" % (self.outprefix)
        self.sample_ids_f = "%s.sample_ids.tsv" % (self.outprefix)
        self.pair_ids_f = "%s.pair_ids.tsv" % (self.outprefix)
        self.block_pairs_f = "%s.block.pairs.tsv" % (self.outprefix)
        self.block_samples_f = "%s.block.samples.tsv" % (self.outprefix)
        self.bed_coverage_f = "%s.bed_coverage.txt" % (self.outprefix)
        self.variant_blocks_tsv_f = "%s.variant.blocks.tsv" % (self.outprefix)
        self.variant_pairs_tsv_f = "%s.variant.pairs.tsv" % (self.outprefix)
        self.variant_pairs_sfs_tally_f = "%s.variant.pairs.sfs_tally.pickle" % (self.outprefix)
        self.variant_total_sfs_tally_f = "%s.variant.total.sfs_tally.pickle" % (self.outprefix)
        self.window_coverage_tsv_f = "%s.window.coverage.tsv" % (self.outprefix)
        self.window_bed_f = "%s.window.bed" % (self.outprefix)
        self.window_variant_tsv_f = "%s.window.variant.tsv" % (self.outprefix)
        self.window_sfs_tally_f = "%s.window.sfs_tally.pickle" % (self.outprefix)

    def get_population_by_sample_id(self):
        _population_by_sample_id = {}
        for population_id, sample_ids in sorted(self.sample_ids_by_population.items()):
            for sample_id in sample_ids:
                _population_by_sample_id[sample_id] = population_id
        return _population_by_sample_id

    def write_sample_ids(self):
        sample_ids_line = []
        sample_ids_line.append("\t".join(["sample_idx", "sample_id", "population"]))
        for sample_id in self.sample_ids:
            sample_ids_line.append("\t".join([str(x) for x in [ \
                self.sample_idx_by_sample_id[sample_id], \
                sample_id, \
                self.population_by_sample_id[sample_id] \
                ]]) \
            )
        with open(self.sample_ids_f, 'w') as sample_ids_fh:
            sample_ids_fh.write("\n".join(sample_ids_line))
        return self.sample_ids_f

    def write_pair_ids(self):
        pair_ids_line = []
        pair_ids_line.append("\t".join(["pair_idx", "pair_ids"]))
        for pair_idx in self.pair_idxs:
            pair_ids_line.append("\t".join([str(x) for x in [ \
                pair_idx, \
                ",".join(self.pair_ids_by_pair_idx[pair_idx])]]) \
            )
        with open(self.pair_ids_f, 'w') as pair_ids_fh:
            pair_ids_fh.write("\n".join(pair_ids_line))
        return self.pair_ids_f

class RegionBatchObj(object):
    __slots__ = ["contig_id", "bedObjs", "idx"]

    def __init__(self, idx):
        self.contig_id = None
        self.bedObjs = deque()
        self.idx = idx

    def __nonzero__(self):
        if self.contig_id:
            return True
        return False

    def __str__(self):
        try:
            return "[I] : I=%s\tc=%s\tS=%s\tE=%s\tbedObjs=%s\tlength=%s\tspan=%s" % (self.idx, self.contig_id, self.bedObjs[0].start, self.bedObjs[-1].end, len(self.bedObjs), self.length(), self.span())
        except IndexError:
            return "[I] : I=%s\tc=%s\tS=%s\tE=%s\tbedObjs=%s\tlength=%s\tspan=%s" % (self.idx, self.contig_id, "?", "?", len(self.bedObjs), self.length(), self.span())

    def add_bedObj_to_batch(self, bedObj):
        if self.contig_id == None:
            self.contig_id = bedObj.chrom
        self.bedObjs.append(bedObj)

    def length(self):
        try:
            return sum([bedObj.length for bedObj in self.bedObjs])
        except TypeError:
            return 0

    def span(self):
        try:
            return (self.bedObjs[-1].end - self.bedObjs[0].start)
        except IndexError:
            return 0

class BedObj(object):
    __slots__ = ['chrom', 'start', 'end', 'pair_idxs', 'length']
    
    def __init__(self, chrom, start, end, pair_idxs, length):
        self.chrom = chrom
        self.start = int(start)
        self.end = int(end)
        self.length = int(length)
        self.pair_idxs = set(pair_idxs)

    def __str__(self):
        return "\t".join([self.chrom, str(self.start), str(self.end), str(self.length), str(self.pair_idxs)]) 

class BlockObj(object):

    __slots__ = [\
        "contig_id", \
        "block_id", \
        "pair_idxs", \
        "sample_idxs", \
        "start", \
        "end", \
        "length", \
        "span", \
        "score", \
        "needed", \
        "bed_tuples", \
        "genotypes_by_sample_idx", \
        "profileObj_by_pair_idx", \
        "void" \
        ]

    def __init__(self, block_id, block_length):
        self.contig_id = block_id.split(".")[0]
        self.block_id = block_id
        self.pair_idxs = None
        self.sample_idxs = None
        self.start = None
        self.end = None
        self.length = 0
        self.span = 0 
        self.score = 0.0
        self.needed = block_length
        self.bed_tuples = [] # list of tuples (contig_id, start, end) of consecutive regions

        self.genotypes_by_sample_idx = {}  # dict of lists
        self.profileObj_by_pair_idx = {}  # dict of ProfileObjs 

        self.void = False

    def __str__(self):
        return "[B] ID=%s %s %s %s LEN=%s SPAN=%s SCORE=%s [P]=%s" % (self.block_id, self.contig_id, self.start, self.end, self.length, self.span, self.score, self.pair_idxs)

    def __nonzero__(self):
        if self.length:
            return True
        return False

    def add_bedObj(self, bedObj, parameterObj):
        '''
        Function for adding a bedObj to the blockObj

        [parameters]
            - bedObj to be added
            - parameterObj

        [returns]
            a) None (if bedObj has been consumed)
            b) bedObj
                b1) original bedObj (if span-violation)
                b2) remainder bedObj (if bedObj.length > blockObj.needed)
        
        [comments]
            - span-violation:
                if blockObj.span > parameterObj.block_span:
                    - original bedObj is returned, blockObj.score is set to 0.0
            - blockObj.needed: allows distinction between 
                a) finished block: blockObj.needed = 0
                b) virgin block: blockObj.needed = parameterObj.block_length
                c) started block: 0 < blockObj.needed < parameterObj.block_length
            - blockObj.score: 
                a) if blockObj.needed == parameterObj.block_length (virgin block):
                    blockObj.score = (bedObj.pair_idxs / parameterObj.pairs_count) * (min(bedObj.length, required_length) / parameterObj.block_length)
                b) if blockObj.needed < parameterObj.block_length (non-virgin block):
                    blockObj.score = (len(blockObj.pair_idxs.intersection(bedObj.pair_idxs)) / parameterObj.pairs_count) * (min(bedObj.length, required_length) / parameterObj.block_length)
                c) if span-violation:
                    blockObj.score = 0.0
        '''
        interval_length = min(self.needed, bedObj.length)
        block_end = bedObj.start + interval_length
        try:
            _span = block_end - self.start # TypeError: int() argument must be a string, a bytes-like object or a number, not 'NoneType' 
        except TypeError:
            self.start = bedObj.start
            _span = block_end - self.start
        if _span > parameterObj.block_span:
            self.score = 0.0
            return bedObj
        try:
            self.pair_idxs = self.pair_idxs.intersection(bedObj.pair_idxs) # AttributeError: 'NoneType' object has no attribute 'intersection'
        except AttributeError:
            self.pair_idxs = bedObj.pair_idxs

        self.end = block_end 
        self.span = _span
        self.length += interval_length
        self.needed -= interval_length
        self.score = (len(self.pair_idxs) / parameterObj.pairs_count) * (self.length / parameterObj.block_length)
        self.sample_idxs = pairs_to_samples(self.pair_idxs, parameterObj)   
        try:
            last_tuple = self.bed_tuples.pop()
            if (last_tuple[2] - bedObj.start) == 0: # no gap
                self.bed_tuples.append((bedObj.chrom, last_tuple[1], block_end)) 
            else: # gap
                self.bed_tuples.append(last_tuple)
                self.bed_tuples.append((bedObj.chrom, bedObj.start, block_end)) 
        except IndexError:
            self.bed_tuples.append((bedObj.chrom, bedObj.start, block_end))
        self.contig_id = bedObj.chrom
        if interval_length < bedObj.length:
            return BedObj(bedObj.chrom, (bedObj.start + interval_length), bedObj.end, bedObj.pair_idxs, (bedObj.length - interval_length))    
        else:
            return None
            

class BlockDataObj(object):
    __slots__ = ["blockObjs", "blockObjs_by_block_id", "blockObj_idxs", "idx_list_by_pair_idx"]

    def __init__(self, parameterObj):
        self.blockObjs = []
        self.blockObjs_by_block_id = {}
        self.blockObj_idxs = [] # list of start_idxs of contigs, hacky ...
        self.idx_list_by_pair_idx = OrderedDict({pair_idx : [] for pair_idx in parameterObj.pair_idxs})

    def __len__(self):
        return len(self.blockObjs)

    def write_profiles_by_pairs(self, parameterObj):
        data_variant_blocks = []
        data_variant_blocks.append("%s" % ("\t".join(["block_id", "pair_idx", "fixed", "hetA", "hetAB", "hetB", "missing", "multiallelic"])))
        profileObjs_by_pair_idx = defaultdict(list)
        for blockObj in self.blockObjs:
            for pair_idx, profileObj in blockObj.profileObj_by_pair_idx.items():
                #print("\t".join([str(x) for x in profileObj.tuple()]) )
                profileObjs_by_pair_idx[pair_idx].append(profileObj)
                data_variant_blocks.append("\t".join([ \
                    blockObj.block_id, \
                    str(pair_idx), \
                    "\t".join([str(x) for x in profileObj.tuple()]) \
                ]))
        
        fn_variant_blocks_tsv = parameterObj.variant_blocks_tsv_f
        with open(fn_variant_blocks_tsv, 'w') as fh_variant_blocks_tsv:
            fh_variant_blocks_tsv.write("\n".join(data_variant_blocks) + "\n")
        return fn_variant_blocks_tsv, profileObjs_by_pair_idx

    def add_blockObjs(self, blockObjs):
        '''
        takes list of blockObjs
        ''' 
        blockObjs.sort(key=lambda i: (i.contig_id, i.start)) # so that they are ordered by contig/start ... 
        last_contig_id, start, end = blockObjs[0].contig_id, 0, 0
        for idx, blockObj in enumerate(blockObjs):
            if not blockObj.contig_id == last_contig_id:
                end = idx
                self.blockObj_idxs.append((start, end))
                start = idx
                last_contig_id = blockObj.contig_id
            self.blockObjs.append(blockObj)
            self.blockObjs_by_block_id[blockObj.block_id] = blockObj
            for pair_idx in blockObj.pair_idxs:
                self.idx_list_by_pair_idx[pair_idx].append(len(self.blockObjs) - 1)
        self.blockObj_idxs.append((start, len(self.blockObjs)))
        
    def write_block_pairs(self, parameterObj):
        lines_block_pairs = []
        lines_block_pairs.append("%s" % ("\t".join(["pair_idx", "bases", "sample_ids"])))
        lines_block_samples = []
        lines_block_samples.append("%s" % ("\t".join(["sample_idx", "sample_id", "bases_mean", "bases_min", "bases_max"])))
        blocks_by_sample_idx = defaultdict(list)
        for pair_idx, idx_list in self.idx_list_by_pair_idx.items():
            bases = (len(idx_list) * parameterObj.block_length)
            sample_idxs = parameterObj.sample_idxs_by_pair_idx[pair_idx]
            for sample_idx in sample_idxs:
                blocks_by_sample_idx[sample_idx].append(bases)
            lines_block_pairs.append("\t".join([str(x) for x in [pair_idx, bases, ",".join(parameterObj.pair_ids_by_pair_idx[pair_idx])]]))
        for sample_idx, bases_list in sorted(blocks_by_sample_idx.items()):
            lines_block_samples.append("\t".join([str(x) for x in [sample_idx, parameterObj.sample_id_by_sample_idx[sample_idx], "%.2f" % numpy.mean(bases_list), numpy.min(bases_list), numpy.max(bases_list)]]))
        fn_block_pairs = parameterObj.block_pairs_f
        with open(fn_block_pairs, 'w') as fh_block_pairs:
            fh_block_pairs.write("\n".join(lines_block_pairs))
        fn_block_samples = parameterObj.block_samples_f
        with open(fn_block_samples, 'w') as fh_block_samples:
            fh_block_samples.write("\n".join(lines_block_samples))
        return fn_block_pairs

    def write_block_bed(self, parameterObj):
        header_bed = "%s" % ("\t".join(["CHROM", "START", "END", "block_id"]))
        lines_bed = []
        lines_bed.append(header_bed)
        void_lines_bed = []
        void_lines_bed.append(header_bed)
        for blockObj in sorted(self.blockObjs, key=lambda i: (i.contig_id, i.start)):
            if blockObj.void:
                for bed_tuple in blockObj.bed_tuples:
                    void_lines_bed.append("\t".join([str(x) for x in [bed_tuple[0], bed_tuple[1], bed_tuple[2], blockObj.block_id]]))
            else:
                for bed_tuple in blockObj.bed_tuples:
                    lines_bed.append("\t".join([str(x) for x in [bed_tuple[0], bed_tuple[1], bed_tuple[2], blockObj.block_id]]))
        fn_bed = parameterObj.block_bed_f
        if parameterObj.coordinate_f:
            fn_bed = parameterObj.block_bed_fix_f
        with open(fn_bed, 'w') as fh_bed:
            fh_bed.write("\n".join(lines_bed))
        if len(void_lines_bed) > 1:
            fn_block_bed_void = parameterObj.block_void_bed_f
            with open(fn_block_bed_void, 'w') as fh_block_bed_void:
                fh_block_bed_void.write("\n".join(void_lines_bed))
            return fn_bed, fn_block_bed_void
        return fn_bed, None

    def write_block_summary(self, parameterObj):
        lines_tsv = []
        lines_tsv.append("%s" % ("\t".join(["block_id", "length", "span", "count_samples", "count_pairs", "samples", "pairs"])))
        for blockObj in sorted(self.blockObjs, key=lambda i: (i.contig_id, i.start)):
            lines_tsv.append("\t".join([str(x) for x in [blockObj.block_id, blockObj.length, blockObj.span, len(blockObj.sample_idxs), len(blockObj.pair_idxs), ",".join([str(x) for x in sorted(blockObj.sample_idxs)]), ",".join([str(x) for x in sorted(blockObj.pair_idxs)])]]))
        fn_tsv = parameterObj.block_summary_f
        with open(fn_tsv, 'w') as fh_tsv:
            fh_tsv.write("\n".join(lines_tsv))
        return fn_tsv    

    def write_variant_blocks(self, parameterObj):
        data_variant_blocks = []
        data_variant_blocks.append("%s" % ("\t".join(["block_id", "pair_idx", "fixed", "hetA", "hetAB", "hetB", "missing", "multiallelic"])))
        mutuples_by_pair_idx = defaultdict(list)
        profileObjs_by_pair_idx = defaultdict(list)
        for blockObj in self.blockObjs:
            for pair_idx, profileObj in blockObj.profileObj_by_pair_idx.items():
                #print("\t".join([str(x) for x in profileObj.tuple()]) )
                mutuples_by_pair_idx[pair_idx].append(profileObj.mutuple())
                profileObjs_by_pair_idx[pair_idx].append(profileObj)
                data_variant_blocks.append("\t".join([ \
                    blockObj.block_id, \
                    str(pair_idx), \
                    "\t".join([str(x) for x in profileObj.tuple()]) \
                ]))
        
        fn_variant_blocks_tsv = parameterObj.variant_blocks_tsv_f
        with open(fn_variant_blocks_tsv, 'w') as fh_variant_blocks_tsv:
            fh_variant_blocks_tsv.write("\n".join(data_variant_blocks) + "\n")
        
        sfs_by_pair_idx = {pair_idx : dict(Counter(mutuples_by_pair_idx[pair_idx])) for pair_idx in parameterObj.pair_idxs}
        with open(parameterObj.variant_pairs_sfs_tally_f, "wb") as pickle_out:
            pickle.dump(sfs_by_pair_idx, pickle_out)
        sfs_total = dict(sum([Counter(mutuples_by_pair_idx[pair_idx]) for pair_idx in parameterObj.pair_idxs], Counter()))
        print(sfs_total)
        with open(parameterObj.variant_total_sfs_tally_f, "wb") as pickle_out:
            pickle.dump(sfs_total, pickle_out)
        return fn_variant_blocks_tsv, profileObjs_by_pair_idx

    def write_variant_summary(self, parameterObj, profileObjs_by_pair_idx):
        data_variant_summary = []
        data_variant_summary.append("%s" % ("\t".join(["pair_idx", "blocks", "bases", "fixed", "hetA", "hetAB", "hetB", "missing", "multiallelic", "pi_A", "pi_B", "d_xy", "f_st"])))
        for pair_idx in parameterObj.pair_idxs:
            profileObjs = profileObjs_by_pair_idx[pair_idx]
            bases = len(profileObjs) * parameterObj.block_length
            normed_profileObj = sum(profileObjs) / bases
            metrics = calculate_metrics(normed_profileObj, 1)
            data_variant_summary.append("\t".join([ \
                str(pair_idx), \
                str(len(profileObjs)), \
                str(bases), \
                "\t".join([str(x) for x in normed_profileObj.tuple()]), \
                str(metrics.get('pi_A', "NA")), \
                str(metrics.get('pi_B', "NA")), \
                str(metrics.get('d_xy', "NA")), \
                str(metrics.get('f_st', "NA")) \
                ]))
        fn_profiles_summary_tsv = parameterObj.variant_pairs_tsv_f
        with open(fn_profiles_summary_tsv, 'w') as fh_profiles_summary_tsv:
            fh_profiles_summary_tsv.write("\n".join(data_variant_summary) + "\n")
        return fn_profiles_summary_tsv
        
class WindowDataObj(object):
    __slots__ = ["windowObjs"]

    def __init__(self):
        self.windowObjs = OrderedDict()

    def __len__(self):
        return len(self.windowObjs)

    def add_windowObjs(self, windowObjs):
        for windowObj in windowObjs:
            self.windowObjs[windowObj.window_id] = windowObj

    def write_window_output(self, parameterObj):
        # plot actual profiles of windows
        fn_window_metrics_tsv = parameterObj.window_coverage_tsv_f
        data_window_metrics_tsv = []
        data_window_metrics_tsv.append("%s" % "\t".join(["window_id", "length", "span", "mean_block_density", "mean_sample_count", "mean_pair_count"] + parameterObj.sample_ids))

        fn_window_bed = parameterObj.window_bed_f
        data_window_bed = []
        data_window_bed.append("%s" % "\t".join(["CHROM", "START", "END", "window_id"]))
        
        fn_window_sfs_tally = parameterObj.window_sfs_tally_f
        data_window_sfs_tally = {}

        fn_window_profiles_tsv = parameterObj.window_variant_tsv_f
        data_window_profiles_tsv = []
        data_window_profiles_tsv.append("%s" % ("\t".join([ \
            "window_id", \
            "hetA", \
            "hetB", \
            "hetAB", \
            "fixed", \
            "multiallelic", \
            "missing", \
            "pi_A", \
            "pi_B", \
            "d_xy", \
            "f_st" \
            ])))
        for window_id, windowObj in tqdm(self.windowObjs.items(), total=len(self.windowObjs), desc="[%] ", ncols=200, unit_scale=True):
            sample_covs = []
            data_window_sfs_tally[window_id] = dict(windowObj.sfs_tally)
            for sample_id in parameterObj.sample_ids:
                _sample_idx = parameterObj.sample_idx_by_sample_id[sample_id]
                sample_covs.append("%.2f" % (windowObj.cov_by_sample_idx.get(_sample_idx, 0) / parameterObj.window_size))
            data_window_metrics_tsv.append( \
                "%s" % "\t".join([\
                    window_id, \
                    str(windowObj.length), \
                    str(windowObj.span), \
                    "%.4f" % (windowObj.length / windowObj.span), \
                    "%.2f" % (sum([cov for sample_idx, cov in windowObj.cov_by_sample_idx.items()]) / (parameterObj.window_size * parameterObj.samples_count)), \
                    "%.2f" % (sum([cov for pair_idx, cov in windowObj.cov_by_pair_idx.items()]) / (parameterObj.window_size * parameterObj.pairs_count)), \
                    "\t".join(sample_covs) \
                ]))
            for contig_id, start, end in windowObj.bed_tuples:
                data_window_bed.append("\t".join([contig_id, str(start), str(end), window_id]))
            data_window_profiles_tsv.append("\t".join([ \
                window_id, \
                str(windowObj.profileObj.hetA / parameterObj.block_length), \
                str(windowObj.profileObj.hetB / parameterObj.block_length), \
                str(windowObj.profileObj.hetAB / parameterObj.block_length), \
                str(windowObj.profileObj.fixed / parameterObj.block_length), \
                str(windowObj.profileObj.multiallelic / parameterObj.block_length), \
                str(windowObj.profileObj.missing / parameterObj.block_length), \
                str(windowObj.metrics['pi_A']), \
                str(windowObj.metrics['pi_B']), \
                str(windowObj.metrics['d_xy']), \
                str(windowObj.metrics['f_st']) \
            ]))
        #with open(fn_window_sfs_tally, 'w') as fh_window_sfs_tally:
        #    fh_window_sfs_tally.write("\n".join(data_window_sfs_tally))
        #numpy.save(fn_window_sfs_tally, data_window_sfs_tally)
        print(fn_window_sfs_tally)
        with open(fn_window_sfs_tally, "wb") as fh_window_sfs_tally:
            pickle.dump(data_window_sfs_tally, fh_window_sfs_tally)
        with open(fn_window_metrics_tsv, 'w') as fh_window_metrics_tsv:
            fh_window_metrics_tsv.write("\n".join(data_window_metrics_tsv))
        with open(fn_window_bed, 'w') as fh_window_bed:
            fh_window_bed.write("\n".join(data_window_bed))
        with open(fn_window_profiles_tsv, 'w') as fh_window_profiles_tsv:
            fh_window_profiles_tsv.write("\n".join(data_window_profiles_tsv))
        return fn_window_metrics_tsv, fn_window_bed, fn_window_profiles_tsv, fn_window_sfs_tally

class WindowObj(object):
    __slots__ = ["contig_id", \
                 "start", \
                 "end", \
                 "centre", \
                 "length", \
                 "span", \
                 "window_id", \
                 "profileObj", \
                 "profileObjs", \
                 "profile_weights", \
                 "metrics", \
                 "sfs_tally", \
                 "cov_by_sample_idx", \
                 "cov_by_pair_idx", \
                 "bed_tuples"]

    def __init__(self, contig_id, blockObjs, block_length):
        self.contig_id = contig_id
        self.start = blockObjs[0].start
        self.end = blockObjs[-1].end
        self.length = block_length * len(blockObjs)
        self.span = blockObjs[-1].end - blockObjs[0].start
        self.window_id = "%s_%s_%s" % (contig_id, blockObjs[0].start, blockObjs[-1].end)
        # => populate
        self.centre = None
        self.profileObj = None
        self.profileObjs = []
        self.profile_weights = []
        self.metrics = {}
        self.sfs_tally = None
        self.cov_by_sample_idx = {}
        self.cov_by_pair_idx = {}
        self.bed_tuples = []
        self.populate(blockObjs, block_length)

    def __str__(self):
        return self.window_id

    def populate(self, blockObjs, block_length):
        sfs_list = []
        centre_list = []
        for blockObj in blockObjs:
            self.profile_weights.append(len(blockObj.pair_idxs))
            block_centre = blockObj.start + (block_length / 2)
            centre_list.append(block_centre)
            sum_profileObj = []
            for pair_idx in blockObj.pair_idxs:
                sum_profileObj.append(blockObj.profileObj_by_pair_idx[pair_idx])
                sfs_list.append(blockObj.profileObj_by_pair_idx[pair_idx].mutuple())
                self.cov_by_pair_idx[pair_idx] = self.cov_by_pair_idx.get(pair_idx, 0) + 1
            self.profileObjs.append(sum(sum_profileObj))
            for sample_idx in blockObj.sample_idxs:
                self.cov_by_sample_idx[sample_idx] = self.cov_by_sample_idx.get(sample_idx, 0) + 1
            for bed_tuple in blockObj.bed_tuples:
                self.bed_tuples.append(bed_tuple)
        self.sfs_tally = Counter(sfs_list)
        self.centre = numpy.median(centre_list)

class ProfileObj():
    
    __slots__ = ['fixed', 'hetA', 'hetAB', 'hetB', 'missing', 'multiallelic']

    '''
    haploid : ['fixed', 'missing', 'multiallelic']
    n haploids : ['fixed', 'missing', 'multiallelic']
    '''
    def __init__(self, profile_tuple):
        self.fixed = profile_tuple[0]
        self.hetA = profile_tuple[1]
        self.hetAB = profile_tuple[2]
        self.hetB = profile_tuple[3]
        self.missing = profile_tuple[4]
        self.multiallelic = profile_tuple[5]

    def __radd__(self, other):       
        return self.__add__(other)   

    def __add__(self, other):
        if isinstance(other, ProfileObj):
            return ProfileObj(([x + y for x, y in zip(self.tuple(), other.tuple())]))
        else:
            return ProfileObj((0, 0, 0, 0, 0, 0, 0))

    def mutuple(self):
        # This one gets written to SFStallies ...
        return (self.hetA, self.fixed, self.hetB, self.hetAB) # 0.8.0
        # return (self.fixed, self.hetA, self.hetAB, self.hetB) # 0.7.0

    def tuple(self):
        return (self.fixed, self.hetA, self.hetAB, self.hetB, self.missing, self.multiallelic)

    def __str__(self):
        return str(self.tuple())

    def __mul__(self, integer):
        return ProfileObj([x * integer for x in self.tuple()])

    def __truediv__(self, integer):
        return ProfileObj([x / integer for x in self.tuple()])

class BlockLengthException(Exception):
    pass

class PopsCountException(Exception):
    pass

if __name__ == "__main__":
    pass