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

# PLOTs

# Histogram of variants (sum of mutype counts) per block
# distribution of window spans
# fix offset of genome scan


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
    fig = plt.figure(figsize=(16,4), dpi=200, frameon=True)
    #connecting dots
    ax = fig.add_subplot(111)  
    y_lim = (0.0, 1.0)
    scatter = ax.scatter(window_df['centre'], window_df['fst'], c=window_df['dxy'], alpha=1.0, cmap='bone', marker='o', s=2, linewidth=0)
    cbar = fig.colorbar(scatter, ax=ax)
    cbar.ax.set_title('D_xy')
    ax.vlines(x_boundaries, 0.0, 1.0, colors=['lightgrey'], linestyles='dashed', linewidth=1)
    ax.set_ylim(y_lim)
    plt.ylabel('F_st')
    plt.xlabel("Genome coordinate")
    print("[>] Created: %r" % str(out_f))
    ax.autoscale_view(tight=None, scalex=True, scaley=True)
    fig.savefig(out_f, format="png")
    plt.close(fig)

def plot_pi_scatter(window_df, out_f):
    fig, ax = plt.subplots(nrows=1, ncols=2, sharey=True, sharex=True, figsize=(12,6))
    #ax.set_ylim((-0.05, 1))
    #ax.set_xlim((0, categories + 1))
    #ax.set_xscale("log")
    print(list(window_df.columns))
    pi_A_key = list(window_df.columns)[5]
    pi_B_key = list(window_df.columns)[6]
    ax[0].scatter(window_df['dxy'], window_df[pi_A_key], c=window_df['fst'], cmap='Oranges', marker='o', s=100, alpha=0.5, label=pi_A_key)
    scatter = ax[1].scatter(window_df['dxy'], window_df[pi_B_key], c=window_df['fst'], cmap='Oranges', marker='o', s=100, alpha=0.5, label=pi_B_key)
    cbar = fig.colorbar(scatter, ax=ax)
    cbar.ax.set_title('F_st')
    ax[0].set_ylabel('D_xy')
    ax[0].set_xlabel(pi_A_key)
    ax[1].set_xlabel(pi_B_key)
    #plt.legend(frameon=True)
    #ax.set_aspect(1.)
    #divider = make_axes_locatable(ax)
    #axHistx = divider.append_axes("top", 1.2, pad=0.2, sharex=ax)
    #axHisty = divider.append_axes("right", 1.2, pad=0.2, sharey=ax)
    #axHistx.xaxis.set_tick_params(labelbottom=False)
    #axHisty.yaxis.set_tick_params(labelleft=False)
    #axHistx.hist(window_df['dxy'], bins=100, color='deeppink')
    #axHistx.hist(window_df['dxy'], bins=100, color='dodgerblue')
    #axHisty.hist(window_df[pi_A_key], bins=100, color='deeppink', orientation='horizontal')
    #axHisty.hist(window_df[pi_B_key], bins=100, color='dodgerblue', orientation='horizontal')
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
    plt.ylabel(y_label)
    plt.xlabel(x_label)
    ax.autoscale_view(tight=None, scalex=True, scaley=True)
    plt.xticks(rotation=90)
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