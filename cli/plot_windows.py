#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""usage: plot_windows.py                -b FILE -g FILE [-h|--help]
                                            
    Options:
        -b, --bed_f FILE                 BED file
        -g, --genome_f FILE              GENOME file
        -h --help                        show this
"""

'''
ylim_min for fst IS 0
dxy
het_a + het_b
het_a + het_b + dxy

'''

import pandas as pd
from timeit import default_timer as timer
from docopt import docopt
import collections
import matplotlib.pyplot as plt
import matplotlib as mat
import sys
mat.use("agg")

COLOR_HISTOGRAM = 'orange'
COLORS = ['deeppink', 'dodgerblue']
COLOR_AXES = 'grey'

LEGEND_FONTSIZE = 10
AXES_TICKS_FONTSIZE = 8
AXES_LABELS_FONTSIZE = 8

# MATPLOTLIB PARAMS
mat.rcParams['text.color'] = COLOR_AXES
mat.rcParams['axes.edgecolor'] = COLOR_AXES
mat.rcParams['xtick.color'] = COLOR_AXES
mat.rcParams['ytick.color'] = COLOR_AXES
mat.rcParams['grid.color'] = COLOR_AXES
mat.rcParams['font.family'] = 'sans-serif'
mat.rcParams['xtick.labelsize'] = AXES_LABELS_FONTSIZE
mat.rcParams['ytick.labelsize'] = AXES_LABELS_FONTSIZE
mat.rcParams['figure.frameon'] = False
mat.rcParams['axes.grid'] = False

def plot_popgen_scan(bed_df, offset_by_sequence_id, bg_rectangles, out_prefix):
    bed_df['x'] = bed_df['pos_mean'] + bed_df['sequence'].map(offset_by_sequence_id)
    fig = plt.figure(figsize=(16,4), dpi=200, frameon=True)
    ax = fig.add_subplot(111)
    for axvspan_min, axvspan_max in bg_rectangles:
        plt.axvspan(axvspan_min, axvspan_max, color='whitesmoke', zorder=-1)
    ax.plot(bed_df['x'], bed_df['f_st'], color='black', alpha=1, linestyle='-', linewidth=0.1) 
    ax.scatter(bed_df['x'], bed_df['f_st'], color='black', s=0.1, alpha=1) 
    plt.tight_layout()
    out_f = '%s.f_st.scan.png' % out_prefix
    fig.savefig(out_f)
    print("[+] Created %s" % out_f)
    plt.close(fig)

def plot_barchart_window_span(bed_df, out_prefix):
    fig = plt.figure(figsize=(8,4), dpi=200, frameon=True)
    ax = fig.add_subplot(111)  
    ax.hist(bed_df['span'], bins=100, color='lightgrey', alpha=0.8, edgecolor='white') 
    ax.set_xlabel('window span (b)')
    ax.set_ylabel('count')
    plt.axvline(bed_df['span'].median(), color='black', linestyle='dashed', linewidth=1)
    min_ylim, max_ylim = plt.ylim()
    plt.text(bed_df['span'].median()*1.1, max_ylim*0.9, 'Median: {:.2f}'.format(bed_df['span'].median()), color='black', fontsize=8)
    plt.tight_layout()
    out_f = '%s.barchart_window_span.png' % out_prefix
    fig.savefig(out_f)
    print("[+] Created %s" % out_f)
    plt.close(fig)

def plot_pi_genome_scan(window_df, out_f, sequenceObjs):
    offset_by_sequence_id = {}
    offset = 0
    x_boundaries = []
    for sequenceObj in sequenceObjs:
        offset_by_sequence_id[sequenceObj.id] = offset
        x_boundaries.append(offset)
        offset += sequenceObj.length

    x_boundaries.append(offset)
    #print([(sequenceObj.id, sequenceObj.length) for sequenceObj in sequenceObjs])
    #print(x_boundaries)
    fig = plt.figure(figsize=(18,4), dpi=200, frameon=True)
    #connecting dots
    ax = fig.add_subplot(111)  
    window_df['rel_pos'] = window_df['centre'] + window_df['sequence_id'].map(offset_by_sequence_id)
    window_df.sort_values(['rel_pos'], inplace=True)
    #print(window_df)
    pi_A_key = list(window_df.columns)[6]
    pi_B_key = list(window_df.columns)[7]
    ax.plot(window_df['rel_pos'], window_df[pi_A_key], color='orange', alpha=0.8, linestyle='-', linewidth=1, label=pi_A_key.replace('pi_', ''))
    ax.plot(window_df['rel_pos'], window_df[pi_B_key], color='dodgerblue', alpha=0.8, linestyle='-', linewidth=1, label=pi_B_key.replace('pi_', ''))
    y_lim = (min(window_df[pi_A_key].min(), window_df[pi_B_key].min()), max(window_df[pi_A_key].max(), window_df[pi_B_key].max()))
    ax.vlines(x_boundaries, y_lim[0], y_lim[1], colors=['lightgrey'], linestyles='dashed', linewidth=1)
    ax.set_ylim(y_lim)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.legend(numpoints=1)
    plt.ylabel('Pi')
    plt.xlabel("Genome coordinate")
    print("[>] Created: %r" % str(out_f))
    plt.tight_layout()
    fig.savefig(out_f, format="png")
    plt.close(fig)

def plot_scatter(bed_df, out_prefix):
    fig = plt.figure(figsize=(18,4), dpi=200, frameon=True)
    ax = fig.add_subplot(111)  
    scatter_parameter_cols = bed_df.columns[4:-1]
    lncl_col = bed_df.columns[-1]
    fig, ax = plt.subplots(nrows=len(scatter_parameter_cols), ncols=1, sharey=False, sharex=True, figsize=(18,2*len(scatter_parameter_cols)))
    for idx, parameter in enumerate(scatter_parameter_cols):
        ax[idx].scatter(bed_df[lncl_col], bed_df[parameter], color='lightgrey', alpha=0.2, label=parameter) 
        ax[idx].set_ylabel(parameter)
    ax[-1].set_xlabel(lncl_col)
    plt.tight_layout()
    out_f = '%s.scatter.png' % out_prefix
    fig.savefig(out_f)
    print("[+] Created %s" % out_f)
    plt.close(fig)

def parse_bed(bed_f):
    print("[+] Parsing BED file ...")
    bed_df = pd.read_csv(
        bed_f, 
        sep="\t",  
        header=1)
    bed_df.columns = [col.lstrip("#").strip() for col in bed_df.columns]
    bed_df['span'] = bed_df['end'] - bed_df['start']
    out_prefix = ".".join(bed_f.split(".")[:-1])
    return (bed_df, out_prefix)

def get_coordinates(genome_f):
    genome_df = pd.read_csv(genome_f, sep="\t", usecols=[0,1], names=['sequence_id', 'sequence_length'], dtype={'sequence_id': 'category', 'sequence_length': 'int64'}, header=None)
    offset_by_sequence_id = {}
    offset = 0
    bg_rectangles = []
    for idx, (sequence_id, sequence_length) in enumerate(list(genome_df.itertuples(index=False, name=None))):
        offset_by_sequence_id[sequence_id] = offset
        if idx % 2 != 0:
            bg_rectangles.append((offset, offset + sequence_length))
        offset += sequence_length
    return (offset_by_sequence_id, bg_rectangles)

def plot_distributions(bed_df, out_prefix):
    parameter_cols = bed_df.columns[4:]
    fig, ax = plt.subplots(nrows=len(parameter_cols), ncols=1, sharey=False, sharex=False, figsize=(5,2*len(parameter_cols)))
    for idx, parameter in enumerate(parameter_cols):
        bins = len(bed_df[parameter].unique())
        ax[idx].hist(bed_df[parameter], bins=bins, color='lightgrey', alpha=0.8, linestyle='-', linewidth=1, label=parameter) 
        ax[idx].set_xlabel(parameter)
        ax[idx].set_ylabel('Count')
    plt.tight_layout()
    out_f = '%s.hist.png' % out_prefix
    fig.savefig(out_f)
    print("[+] Created %s" % out_f)
    plt.close(fig)

if __name__ == '__main__':
    __version__ = 0.1
    try:
        start_time = timer()
        args = docopt(__doc__)
        offset_by_sequence_id, bg_rectangles = get_coordinates(args['--genome_f'])
        bed_df, out_prefix = parse_bed(args['--bed_f'])
        plot_barchart_window_span(bed_df, out_prefix)
        plot_popgen_scan(bed_df, offset_by_sequence_id, bg_rectangles, out_prefix) 
        
    except KeyboardInterrupt:
        exit(-1)