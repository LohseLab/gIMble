import pathlib
import os
from psutil import Process
import multiprocessing
from contextlib import contextmanager

import sys
import collections
import pandas as pd
from tabulate import tabulate

import matplotlib as mat
import matplotlib.pyplot as plt
mat.use("agg")

# CONSTANTS
MUTYPE_ORDER = ['hetA', 'fixed', 'hetB', 'hetAB']

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


def plot_genome_scan(window_df, out_f, sequenceObjs):
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
    y_lim = (0.0, 1.0)
    window_df['rel_pos'] = window_df['centre'] + window_df['sequence_id'].map(offset_by_sequence_id)
    window_df.sort_values(['rel_pos'], inplace=True)
    #print(window_df)
    ax.plot(window_df['rel_pos'], window_df['fst'], color='lightgrey', alpha=0.8, linestyle='-', linewidth=1)
    scatter = ax.scatter(window_df['rel_pos'], window_df['fst'], c=window_df['dxy'], alpha=1.0, cmap='PiYG_r', edgecolors='white', marker='o', s=40, linewidth=0.2)
    cbar = fig.colorbar(scatter, ax=ax)
    cbar.ax.set_title('D_xy')
    ax.vlines(x_boundaries, 0.0, 1.0, colors=['lightgrey'], linestyles='dashed', linewidth=1)
    ax.set_ylim(y_lim)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    plt.ylabel('F_st')
    plt.xlabel("Genome coordinate")
    print("[>] Created: %r" % str(out_f))
    ax.autoscale_view(tight=None, scalex=True, scaley=True)
    fig.savefig(out_f, format="png")
    plt.close(fig)

def plot_parameter_scan(composite_likelihood_df, parameterObj):
    fig, ax = plt.subplots(nrows=2, ncols=1, sharey=False, sharex=True, figsize=(18,4))
    #y_lim = (0.0, 1.0)

    best_composite_likelihood_df = composite_likelihood_df[composite_likelihood_df['ismax'] == True]
    theta_ancestral_vals = []
    theta_derived_vals = []
    migration_vals = []
    position_vals = []
    for window_id, cL, theta_ancestral, theta_derived, migration, ismax in best_composite_likelihood_df.values:
        theta_ancestral_vals.append(float(theta_ancestral))
        theta_derived_vals.append(float(theta_derived))
        migration_vals.append(float(migration))
        position_vals.append(parameterObj.window_pos_by_window_id[window_id])
    #print(theta_ancestral_vals)
    ax[0].plot(position_vals, theta_derived_vals, color='dodgerblue', alpha=0.8, linestyle='-', linewidth=1, label='derived')
    ax[0].plot(position_vals, theta_ancestral_vals, color='orange', alpha=0.8, linestyle='-', linewidth=1, label='ancestral')
    ax[0].set_ylabel('theta')
    ax[1].plot(position_vals, migration_vals, color='purple', alpha=0.8, linestyle='-', linewidth=1)
    ax[1].set_ylabel('Migration')
    #scatter = ax.scatter(window_df['rel_pos'], window_df['fst'], c=window_df['dxy'], alpha=1.0, cmap='PiYG_r', edgecolors='white', marker='o', s=40, linewidth=0.2)
    #cbar = fig.colorbar(scatter, ax=ax)
    #cbar.ax.set_title('D_xy')
    ax[0].vlines(parameterObj.x_boundaries, min(min(theta_ancestral_vals), min(theta_derived_vals)), max(max(theta_derived_vals), max(theta_ancestral_vals)), colors=['lightgrey'], linestyles='dashed', linewidth=1)
    ax[1].vlines(parameterObj.x_boundaries, min(migration_vals), max(migration_vals), colors=['lightgrey'], linestyles='dashed', linewidth=1)
    #ax.set_ylim(y_lim)
    ax[0].spines['right'].set_visible(False)
    ax[0].spines['top'].set_visible(False)
    ax[1].spines['right'].set_visible(False)
    ax[1].spines['top'].set_visible(False)
    ax[0].legend(numpoints=1)
    #ax[0].get_yaxis().set_major_formatter(
    #    mat.ticker.FormatStrFormatter("%2f"))
    ax[1].get_yaxis().set_major_formatter(
        mat.ticker.FormatStrFormatter("%.2e"))
    plt.xlabel("Genome coordinate")
    out_f = "%s.parameter_scan.png" % parameterObj.dataset
    #ax.autoscale_view(tight=None, scalex=True, scaley=True)
    plt.tight_layout()
    print("[>] Created: %r" % str(out_f))
    fig.savefig(out_f, format="png")
    plt.close(fig)

def plot_pi_scatter(window_df, out_f):
    fig, ax = plt.subplots(nrows=1, ncols=2, sharey=True, sharex=True, figsize=(18,6))
    #print(list(window_df.columns))
    pi_A_key = list(window_df.columns)[6]
    pi_B_key = list(window_df.columns)[7]
    ax[0].scatter(window_df[pi_A_key], window_df['dxy'], c=window_df['fst'], cmap='PiYG_r', marker='o', s=20, alpha=0.9, label=pi_A_key)
    scatter = ax[1].scatter(window_df[pi_B_key], window_df['dxy'], c=window_df['fst'], cmap='PiYG_r', marker='o', s=20, alpha=0.9, label=pi_B_key)
    cbar = fig.colorbar(scatter, ax=ax)
    cbar.ax.set_title('F_st')
    ax[0].set_ylabel('D_xy')
    ax[0].set_xlabel("%s (Pi)" % pi_A_key.replace('pi_', ''))
    ax[1].set_xlabel("%s (Pi)" % pi_B_key.replace('pi_', ''))
    fig.savefig(out_f, format="png")
    print("[>] Created: %r" % str(out_f))
    plt.close(fig)

def plot_sample_barchart(out_f, x_label, barchart_y_vals, barchart_x_vals, barchart_labels, barchart_colours, barchart_populations):
    fig = plt.figure(figsize=(12,10), dpi=200, frameon=True)
    ax = fig.add_subplot(111)
    ax.bar(barchart_x_vals, barchart_y_vals, color=barchart_colours)
    ax.get_yaxis().set_major_formatter(
        mat.ticker.FuncFormatter(lambda x, p: format(int(x), ',')))
    plt.xticks(barchart_x_vals, barchart_labels, rotation=45)
    plt.ylabel('Bases in blocks')
    plt.xlabel('Sample ID')
    ax.autoscale_view(tight=False, scalex=True, scaley=True)
    population_ids = [x for i, x in enumerate(barchart_populations) if i == barchart_populations.index(x)]
    colours = [x for i, x in enumerate(barchart_colours) if i == barchart_colours.index(x)]
    legend_markers = [plt.Line2D([0,0],[0,0], color=colour, marker='o', linestyle='') for colour in colours]
    ax.legend(legend_markers, population_ids, numpoints=1)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    #plt.tight_layout()
    fig.savefig(out_f, format="png")
    print("[>] Created: %r" % str(out_f))
    plt.close(fig)

def plot_mutuple_barchart(out_f, global_mutuple_counter):
    # y axis is percentage of blocks, x axis is mutype, bar label is count
    
    total_counts = sum(global_mutuple_counter.values())
    y_vals = []
    x_vals = []
    #bar_labels = []
    x_labels = []
    max_idx = 30 
    for idx, (mutuple, count) in enumerate(global_mutuple_counter.most_common()):
        if idx < max_idx:
            x_vals.append(idx)
            #y_vals.append(count/total_counts)
            y_vals.append(count)
            #bar_labels.append(count)
            x_labels.append("[%s]" % ", ".join([str(count) for count in mutuple]))
        elif idx == max_idx:
            x_vals.append(idx)
            #y_vals.append(count/total_counts)
            y_vals.append(count)
            #bar_labels.append(count)
            x_labels.append('other')
        else:
            #y_vals[-1] += count/total_counts
            y_vals[-1] += count
            #bar_labels[-1] += count

    fig = plt.figure(figsize=(16, 4), dpi=200, frameon=True)
    ax = fig.add_subplot(111)
    ax.bar(x_vals, y_vals, color='orange')
    for idx, p in enumerate(ax.patches):
        width, height = p.get_width(), p.get_height()
        x, y = p.get_xy() 
        ax.annotate('{:.0%}'.format(height/total_counts), (x + width / 2, y + height + 0.01), ha='center', va='bottom', fontsize=AXES_LABELS_FONTSIZE)
        #ax.annotate(format_count(bar_labels[idx]), (x + width / 2, y + height + 0.01), ha='center', va='bottom', rotation=45)
    ax.get_yaxis().set_major_formatter(
        mat.ticker.FuncFormatter(lambda x, p: format(int(x), ',')))
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    plt.xticks(x_vals, x_labels, rotation=45, fontsize=8)
    plt.xlabel('Mutuples [hetA, fixed, hetB, hetAB]')
    plt.ylabel('Count')
    #ax.autoscale_view(tight=False, scalex=True, scaley=True)
    plt.tight_layout()
    fig.savefig(out_f, format="png")
    print("[>] Created: %r" % str(out_f))
    plt.close(fig)

def plot_distance_scatter(out_f, x_label, _x_values):
    y_label = 'Counts'
    x_counter = collections.Counter([int(x) for x in _x_values])
    x_vals, y_vals = [], []

    for _x_val, _y_val in sorted(x_counter.items()):
        x_vals.append(_x_val)
        y_vals.append(_y_val)
    fig = plt.figure(figsize=(12,4), dpi=200, frameon=True)
    ax = fig.add_subplot(111)
    ax.scatter(
        x_vals, 
        y_vals, 
        label=x_label,
        alpha=0.2, 
        s=5,
        facecolor=COLOR_HISTOGRAM, 
        marker='o'
        )
    ax.set_yscale("log")
    ax.set_xscale("log")
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.set_xlim(0.9, max(x_vals))
    plt.ylabel(y_label)
    plt.xlabel(x_label)
    #ax.autoscale_view(tight=True, scalex=True, scaley=True)
    plt.tight_layout()
    fig.savefig(out_f, format="png")
    print("[>] Created: %r" % str(out_f))
    plt.close(fig)

def make_path(path):
    path.mkdir(mode=0o777, parents=True)

def check_file(infile):
    path = pathlib.Path(infile).resolve()
    if path.exists():
        return path
    sys.exit("[X] File not found: %r" % str(infile))

def check_path(prefix):
    if prefix:
        path = pathlib.Path(prefix).resolve()
        if prefix.endswith("/"): # assuming dir prefix
            if not path.exists():
                make_path(path)
            return path
        if not path.parent.exists():
            make_path(path)
        return path.parent
    return None

def check_prefix(prefix):
    if prefix:
        path = pathlib.Path(prefix).resolve()
        if prefix.endswith("/"): # assuming dir prefix
            if not path.exists():
                make_path(path)
            return 'gimble'
        if not path.parent.exists():
            make_path(path)
        return path.name
    return 'gimble'

def memory_usage_psutil():
    # return the memory usage in MB
    return Process(os.getpid()).memory_info()[0] / float(2 ** 20)

def create_csv(out_f, path, header, rows, sep="\t"):
    df = pd.DataFrame(rows, columns=header)
    if not path is None:
        outfile = path.joinpath(out_f)
    else:
        outfile = out_f
    df.to_csv(outfile, index=False, sep=sep)
    print("[>]\tCreated: %r" % str(outfile))

def create_h5_from_df(df, out_f, path, key):
    if not path is None:
        outfile = path.joinpath(out_f)
    else:
        outfile = out_f
    df.to_hdf(outfile, mode='w', format='fixed', key=key)
    print("[>]\tCreated: %r" % str(outfile))
    # pair_values = df.loc(axis=0)[:, [pair_idx], :].values

def tabulate_df(df, columns, title):
    print("\n### %s" % title)
    print_df = pd.DataFrame(df.to_records(), columns=columns)
    print(tabulate(print_df, headers = print_df.columns))

def create_hd5(out_f='', path=None, header=[], rows=[], key='key'):
    df = pd.DataFrame(rows, columns=header)
    if not path is None:
        outfile = path.joinpath(out_f)
    else:
        outfile = out_f
    df.to_hdf(outfile, mode='w', format='fixed', key=key)
    print("[>]\tCreated: %r" % str(outfile))

def create_hdf5_store(out_f, path):
    if not path is None:
        outfile = path.joinpath(out_f)
    else:
        outfile = out_f
    print("[>] Created hdf5 store: %r" % str(outfile))
    return pd.HDFStore(outfile, "w", compression='lzo')

def read_h5(h5_file, key=None):
    return pd.read_hdf(h5_file)

def create(out_f, path, header, lines):
    if out_f:
        if lines:
            out_f = path.joinpath(out_f)
            with open(out_f, 'w') as out_fh:
                if header:
                    out_fh.write("%s\n" % header)
                out_fh.write("%s\n" % "\n".join(lines))
            print("[>]\tCreated: %r" % out_f)
    else:
        print("[>]\tCreated: %r" % path)

@contextmanager
def poolcontext(*args, **kwargs):
    pool = multiprocessing.Pool(*args, **kwargs)
    yield pool
    pool.terminate()

def format_bases(bases):
    return "%s b" % format(bases, ',d')

def format_percentage(fraction, precision=2):
    return "{:.{}%}".format(fraction, precision)

def format_fraction(fraction, precision=2):
    return "{:.{}}".format(fraction, precision)

def format_count(count):
    return "%s" % str(format(count, ',d'))