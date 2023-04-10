#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""usage: plot_gridsearch               -b FILE [-f INT] [-l INT] [-h|--help]
                                            
    Options:
        -b, --bed_f FILE                 BED file
        -f, --column_first INT           First column (0-based) to include in barchart [default: 4]
        -l, --column_last INT            Last column (0-based) to include in barchart 
        -h --help                        show this
"""

import pandas as pd
from timeit import default_timer as timer
from docopt import docopt
import matplotlib.pyplot as plt
import matplotlib as mat
import sys
import collections
mat.use("agg")

COLOR_HISTOGRAM = 'orange'
COLORS = ['deeppink', 'dodgerblue']
COLOR_AXES = 'grey'

LEGEND_FONTSIZE = 10
AXES_TICKS_FONTSIZE = 8
AXES_LABELS_FONTSIZE = 8
MAX_BINS = 500

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

def plot_windows(bed_df, out_prefix):
    sequence_boundaries = bed_df[bed_df['sequence'] != bed_df['sequence'].shift(-1)].index
    parameter_cols = bed_df.columns[4:]
    fig, ax = plt.subplots(nrows=len(parameter_cols), ncols=1, sharey=False, sharex=True, figsize=(18,2*len(parameter_cols)))
    for idx, parameter in enumerate(parameter_cols):
        ax[idx].plot(bed_df[parameter], color='lightgrey', alpha=0.8, linestyle='-', linewidth=1, label=parameter) 
        ax[idx].set_ylabel(parameter)
        y_lim = (bed_df[parameter].min(), bed_df[parameter].max())
        ax[idx].vlines(sequence_boundaries, y_lim[0], y_lim[1], colors=['black'], linestyles='dashed', linewidth=1)
    ax[-1].set_xlabel('window_idx')
    plt.tight_layout()
    out_f = '%s.scan.png' % out_prefix
    fig.savefig(out_f)
    print("[+] Created %s" % out_f)
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
    print("[+] Parsing BED f'ile ...")
    bed_df = pd.read_csv(
        bed_f, 
        sep="\t",  
        header=0)
    bed_df.columns = [col.lstrip("#").strip() for col in bed_df.columns]
    out_prefix = ".".join(bed_f.split(".")[:-1])
    return (bed_df, out_prefix)

def check_parameters(args):
    args_clean = {}
    if not args['--column_last'] is None:
        try:
            args['--column_last'] = int(args['--column_last'])
        except ValueError:
            sys.exit('[X] Invalid value for --column_last : %r' % args['--column_last'])
    try:
        args['--column_first'] = int(args['--column_first'])
    except ValueError:
        sys.exit('[X] Invalid value for --column_first : %r' % args['--column_first'])
    return args

def plot_distributions(bed_df, out_prefix, first_col, last_col):
    parameter_cols = bed_df.columns[first_col:last_col]
    bins_by_parameter = {}
    for idx, parameter in enumerate(parameter_cols):
        bins = len(bed_df[parameter].unique())
        if bins == 1:
            print("[+] Ignoring column %r ... constant value" % parameter)
        elif bins < MAX_BINS:
            print("[+] Plotting column %r" % parameter)
            bins_by_parameter[parameter] = bins
        else:
            print("[+] Plotting column %r ... adjusting bins to 100" % (parameter))
            bins_by_parameter[parameter] = 100
    fig, ax = plt.subplots(nrows=len(bins_by_parameter), ncols=1, sharey=False, sharex=False, figsize=(5,3*len(bins_by_parameter)))
    for idx, (parameter, bins) in enumerate(bins_by_parameter.items()):
        ax[idx].hist(bed_df[parameter], bins=bins, color='lightgrey', alpha=0.8, linestyle='-', linewidth=1, label=parameter) 
        ax[idx].set_xlabel(parameter)
        ax[idx].set_ylabel('Count')
    plt.tight_layout()
    out_f = '%s.hist.png' % out_prefix
    fig.savefig(out_f)
    print("[+] Created %s" % out_f)
    plt.close(fig)

if __name__ == '__main__':
    __version__ = 0.3
    try:
        start_time = timer()
        args = docopt(__doc__)
        args = check_parameters(args)
        bed_df, out_prefix = parse_bed(args['--bed_f'])
        plot_distributions(bed_df, out_prefix, args['--column_first'], args['--column_last'])
        #plot_windows(bed_df, out_prefix)
        #plot_scatter(bed_df, out_prefix)
    except KeyboardInterrupt:
        exit(-1)