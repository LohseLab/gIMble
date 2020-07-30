#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""usage: plot_windows -5 <FILE> -f <FILE> [-h]
    
    -h, --help
    -5, --hd5_f=<FILE>
    -f, --fasta=<FILE>
"""

import pandas as pd
import sys
from docopt import docopt
import os

def parse_fasta(fasta_f):
    if not os.path.isfile(fasta_f):
        sys.exit("[-] File %r could not be found" % (fasta_f))
    seq_by_seq_id = {}
    with open(fasta_f) as fasta_fh:
        seq_id, seq = '', []
        for l in fasta_fh:
            if l[0] == '>':
                if seq_id:
                    seq_by_seq_id[seq_id] = "".join(seq)
                    seq = []
                seq_id = l[1:].rstrip("\n")
            else:
                seq.append(l.rstrip("\n"))
        seq_by_seq_id[seq_id] = "".join(seq)
    if len(seq_by_seq_id) == 0:
        sys.exit("[-] No sequences in %r" % (fasta_f))
    return seq_by_seq_id

if __name__ == '__main__':
    args = docopt(__doc__)
    store_f = args['--hd5_f']
    seq_by_seq_id = parse_fasta(args['--fasta'])
    df = pd.read_hdf(pd.HDFStore(store_f), key='window_metrics')
    for seq_id, seq in tqdm(seq_by_seq_id, total=len(seq_by_seq_id), desc="[%] ", ncols=100, unit_scale=True):
        
    print(df)
        