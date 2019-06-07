#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""usage: hdf5v -f <FILE> [-d <STR> -p <STR> -c -t -5 -h]
    
    -h, --help
    -f, --file=<FILE>
    -d, --table=<STR>                           Name/key of table 
    -p, --prefix=<STR>                          Output prefix
    -c, --csv                                   Convert to CSV
    -t, --tsv                                   Convert to TSV
    -5, --hdf5                                  Export table as separate HDF5
"""

import pandas as pd
import sys
from docopt import docopt
from tabulate import tabulate

def tabulate_df(df, key):
    print("[+] Columns of table %r\n\t[%s]" % (key, ", ".join(list(df.columns))))
    print("[+] Metrics of numeric columns: ")
    print(tabulate(df.describe(), headers = df.columns, tablefmt="orgtbl", floatfmt=".3f"))
    print()
    print("[+] Sum of numeric columns: ")
    sum_df = pd.DataFrame(df.sum(axis=0, numeric_only=True)).T
    print(tabulate(sum_df, headers = sum_df.columns, tablefmt="orgtbl", floatfmt=".3f"))
    print()

def open_store(store_f):
    store = pd.HDFStore(store_f)
    return store, [key.lstrip("/") for key in store.keys()]

def write(df, out_f, key, mode='5'):
    if mode == 'h5':
        out_f = "%s.h5" % out_f
        df.to_hdf(out_f, mode='w', format='fixed', key=key, compression='blosc:zstd')
    elif mode == 'csv':
        out_f = "%s.csv" % out_f
        df.to_csv(out_f, sep=',', index=False)
    elif mode == 'tsv': 
        out_f = "%s.tsv" % out_f
        df.to_csv(out_f, sep='\t', index=False)
    else:
        store.close()
        sys.exit("[X] Mode %r is not supported." % (mode))    
    print("[>]\tCreated: %r" % str(out_f))

if __name__ == '__main__':
    args = docopt(__doc__)
    store_f = args['--file']
    store, keys = open_store(store_f)
    print("[+] Store %r contains the following tables:\n\t%s" % (store_f, [key.lstrip("/") for key in store.keys()]))
    modes = [label for fmt_flag, label in zip([args['--csv'], args['--tsv'], args['--hdf5']], ['csv', 'tsv', 'h5']) if fmt_flag]
    prefix = args['--prefix']
    table = args['--table']
    if table is not None:
        if not table in set(keys):
            store.close()
            sys.exit("[X] Table %r was not found. Please specify one of the following: %s " % (table, keys))    
        else:
            keys = [table]
    for key in keys:
        df = pd.read_hdf(store, key=key)
        out_f = "%s" % (("%s.%s" % (prefix, key) if prefix else key))
        for mode in modes:  
            write(df, out_f, key, mode=mode)
    if len(keys) == 1:
        tabulate_df(df, key)
    store.close()
        