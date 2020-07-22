#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""usage: gimble partitioncds               -f FILE -b FILE -v FILE [-o STR] [-h|--help]
                                            
    Options:
        -h --help                                   show this
        -f, --fasta_file FILE                       FASTA file
        -v, --vcf_file FILE                         VCF file (filtered)
        -b, --bed_file FILE                         BED file
        -o, --outprefix STR                         Outprefix [default: gimble]

"""

'''
awk '$8=="CDS" || $8=="cds"' ~/Dropbox/heliconius_full/current/hmel2_5.chromosomes.annotation.bed > ~/Dropbox/heliconius_full/current/hmel2_5.chromosomes.annotation.CDS.bed
paste <(cut -f1,2 -d"_" ~/Dropbox/heliconius_full/current/hmel2_5.chromosomes.annotation.CDS.bed) <(cut -f5,6 ~/Dropbox/heliconius_full/current/hmel2_5.chromosomes.annotation.CDS.bed) > ~/Dropbox/heliconius_full/current/hmel2_5.chromosomes.annotation.CDS.by_transcript_id.be
'''

from timeit import default_timer as timer
from docopt import docopt
#import lib.gimblelog
from lib.gimble import RunObj
import lib.gimble
import warnings
import numpy as np
import sys
import allel
from tqdm import tqdm
import collections
import pandas as pd
import itertools

DEGENERACIES = [0, 2, 3, 4]

AMINOACID_BY_CODON = {
'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
'TAT': 'Y', 'TAC': 'Y', 'TAA': 'X', 'TAG': 'X',
'TGT': 'C', 'TGC': 'C', 'TGA': 'X', 'TGG': 'W',
'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
}
CODON_DEGENERACY_BY_CODON = {
'TTT': '002', 'TTC': '002', 'TTA': '202', 'TTG': '202', 
'TCT': '004', 'TCC': '004', 'TCA': '004', 'TCG': '004', 
'TAT': '002', 'TAC': '002', 'TAA': '022', 'TAG': '002', 
'TGT': '002', 'TGC': '002', 'TGA': '020', 'TGG': '000', 
'CTT': '004', 'CTC': '004', 'CTA': '204', 'CTG': '204', 
'CCT': '004', 'CCC': '004', 'CCA': '004', 'CCG': '004', 
'CAT': '002', 'CAC': '002', 'CAA': '002', 'CAG': '002', 
'CGT': '004', 'CGC': '004', 'CGA': '204', 'CGG': '204', 
'ATT': '003', 'ATC': '003', 'ATA': '003', 'ATG': '000', 
'ACT': '004', 'ACC': '004', 'ACA': '004', 'ACG': '004', 
'AAT': '002', 'AAC': '002', 'AAA': '002', 'AAG': '002', 
'AGT': '002', 'AGC': '002', 'AGA': '202', 'AGG': '202', 
'GTT': '004', 'GTC': '004', 'GTA': '004', 'GTG': '004', 
'GCT': '004', 'GCC': '004', 'GCA': '004', 'GCG': '004', 
'GAT': '002', 'GAC': '002', 'GAA': '002', 'GAG': '002', 
'GGT': '004', 'GGC': '004', 'GGA': '004', 'GGG': '004'
} 

# http://arep.med.harvard.edu/labgc/adnan/projects/Utilities/revcomp.html
COMPLEMENT = {
    'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 
    'N': 'N', 'Y': 'R', 'R': 'Y', 
    'S': 'S', 'W': 'W', 'K': 'M', 'M': 'K', 
    'B': 'V', 'D': 'H', 'H': 'D', 'V': 'B'
    }

def revcom(sequence):
    revcom_seq = "".join([COMPLEMENT.get(nt.upper(), '') for nt in sequence[::-1]])
    return revcom_seq

def write_df(df, out_f='', sep='\t', header=True, status=True):
    if header == True:
        df.to_csv(out_f, index=False, sep=sep)
    else:
        df.to_csv(out_f, index=False, sep=sep, header=False)
    if status == True:
        print("[+] => Wrote %r" % str(out_f))

def parse_fasta(fasta_file):
    print("[+] Parsing FASTA file...")
    sequence_by_id = {}
    with open(fasta_file) as fh:
        header, seqs = '', []
        for l in fh:
            if l[0] == '>':
                if header:
                    sequence_by_id[header] = ''.join(seqs).upper()
                header, seqs = l[1:-1].split()[0], [] # Header is split at first whitespace
            else:
                seqs.append(l[:-1])
        sequence_by_id[header] = ''.join(seqs).upper()
    print("[+] Found %s sequences..." % (len(sequence_by_id)))
    return sequence_by_id

def parse_bed(bed_file):
    print("[+] Parsing BED file...")
    try:
        bed_df = pd.read_csv(bed_file, sep="\t", names=['sequence_id', 'start', 'end', 'transcript_id', 'score', 'orientation'], 
                                dtype={'sequence_id': str, 'start': np.int64, 'end': np.int64, 'transcript_id': str,  'orientation': str}).sort_values(['sequence_id', 'start'], ascending=[True, True])
    except ValueError:
        sys.exit("[X] BED file %r does not contain the following the columns: 'sequence_id', 'start', 'end', 'transcript_id', 'score', 'orientation'" % (bed_file))
    count_transcript = bed_df['transcript_id'].nunique()
    count_cds = len(bed_df.index)
    print("[+] Found %s CDSs in %s transcripts (%s CDS per transcript)..." % (
                                        lib.gimble.format_count(count_cds), 
                                        lib.gimble.format_count(count_transcript), 
                                        lib.gimble.format_fraction(count_cds / count_transcript)))
    return bed_df

def get_transcripts(parameterObj, sequence_by_id):
    bed_df = parse_bed(parameterObj.bed_file)
    transcriptObjs = []
    for transcript_id, transcript_df in tqdm(bed_df.groupby(['transcript_id']), desc="[%] Reading BED...", ncols=150):
        transcriptObj = TranscriptObj(transcript_id)
        transcriptObj.add_cds_from_df(transcript_df, sequence_by_id)
        if transcriptObj.positions.shape[0] > 6:
            transcriptObjs.append(transcriptObj)
    return transcriptObjs


class TranscriptObj(object):
    def __init__(self, transcript_id):
        self.transcript_id = transcript_id
        self.sequence = None
        self.positions = None
        self.degeneracy = None
        self.orientation = None
        self.sequence_id = None
        self.start = None
        self.end = None
        self.degeneracy_by_sample = collections.defaultdict(list)

    def add_cds_from_df(self, transcript_df, sequence_by_id):
        if not transcript_df['orientation'].nunique():
            sys.exit("[X] More than one orientation found in CDSs of transcript %s\n%s" % (self.transcript_id, transcript_df))
        if not transcript_df['sequence_id'].nunique():
            sys.exit("[X] More than one sequence_id found in CDSs of transcript %s\n%s" % (self.transcript_id, transcript_df))
        self.orientation = list(transcript_df['orientation'].unique())[0]
        self.sequence_id = list(transcript_df['sequence_id'].unique())[0]
        pos_arrays = []
        cds_list = []
        for sequence_id, start, end, _, score, orientation in transcript_df.values.tolist():
            pos_arrays.append(np.arange(start, end))
            cds = sequence_by_id[sequence_id][start:end]
            if orientation == '-':
                cds_list.insert(0, revcom(cds))
            else:
                cds_list.append(cds)
        sequence = "".join(cds_list)
        self.sequence = np.array(list(sequence))
        self.degeneracy = np.concatenate([degeneracy(["".join(sequence[i:i+3])]) for i in range(0, len(sequence), 3)])
        self.positions = np.concatenate(pos_arrays)
        if self.orientation == '-':
            self.positions = self.positions[::-1]
        self.start = np.min(self.positions)
        self.end = np.max(self.positions)

    def has_start(self):
        codon = "".join(self.sequence[0:3])
        if AMINOACID_BY_CODON[codon] == 'M':
            return True
        return False

    def has_stop(self):
        codon = "".join(self.sequence[-3:])
        if AMINOACID_BY_CODON[codon] == 'X':
            return True
        return False

    def is_divisible_by_three(self):
        if self.sequence.shape[0] % 3 == 0:
            return True
        return False

    def is_orf(self):
        if self.has_start() and self.has_stop() and self.is_divisible_by_three():
            return True
        return False

    def __str__(self):
        return ">%s [%s:%s-%s(%s)]\n%s" % (self.transcript_id, self.sequence_id, self.start, self.end, self.orientation, self.sequence)

class ParameterObj(RunObj):
    '''Sanitises command line arguments and stores parameters'''
    def __init__(self, params, args):
        super().__init__(params)
        self.fasta_file = self._get_path(args['--fasta_file'], path=True)
        self.vcf_file = self._get_path(args['--vcf_file'], path=True)
        self.bed_file = self._get_path(args['--bed_file'], path=True)
        self.outprefix = args['--outprefix']

def parse_vcf_file(vcf_file, sequence_ids, query_regions_by_sequence_id):
    samples = []
    variant_arrays_by_seq_id = collections.defaultdict(dict)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for seq_id in tqdm(sequence_ids, total=len(sequence_ids), desc="[%] Parsing VCF file... ", ncols=150):
            data_by_key = allel.read_vcf(str(vcf_file), region=seq_id, fields=['samples', 'calldata/GT', 'variants/POS', 'variants/REF', 'variants/ALT'])
            if data_by_key:
                pos = data_by_key['variants/POS']
                pos_in_query_mask = np.isin(pos, query_regions_by_sequence_id[seq_id])
                variant_arrays_by_seq_id[seq_id]['REF'] = np.array(data_by_key['variants/REF'])[pos_in_query_mask]
                variant_arrays_by_seq_id[seq_id]['ALT'] = np.array(data_by_key['variants/ALT'])[pos_in_query_mask]
                variant_arrays_by_seq_id[seq_id]['GT'] = np.array(data_by_key['calldata/GT'])[pos_in_query_mask]
                variant_arrays_by_seq_id[seq_id]['POS'] = np.array(data_by_key['variants/POS'])[pos_in_query_mask]
                samples = list(data_by_key['samples'])
    print("[+] Parsed %s variants." % sum([variant_arrays_by_seq_id[seq_id]['POS'].shape[0] for seq_id in variant_arrays_by_seq_id]))
    return (samples, variant_arrays_by_seq_id)

def degeneracy(array):
    '''
    Takes: list of codons 
    Returns: list of 'AA_DEGENERACY' ('AA_DEGENERACY|...') strings
    >>> degeneracy(['ATG'])
    ['M_0', 'M_0', 'M_0']
    >>> degeneracy(['AGT', 'ATG'])
    ['M_0|S_0', 'M_0|S_0', 'M_0|S_2']
    >>> degeneracy(['AGT', 'AGC']) 
    ['S_0', 'S_0', 'S_2'] 
    >>> degeneracy(['AGT', 'ATG', 'GGG', 'AAA', 'CCC', 'TTT', 'AGC'])
    ['F_0|G_0|K_0|M_0|P_0|S_0', 'F_0|G_0|K_0|M_0|P_0|S_0', 'F_2|G_4|K_2|M_0|P_4|S_2']
    >>> degeneracy(['AAA', 'XXX'])
    ['NA', 'NA', 'NA']
    '''
    if array:
        AA = [AMINOACID_BY_CODON.get(''.join(row), [])*3 for row in array]
        if not AA or None in AA:
            #print('Array', array)
            #print('AA', AA)
            
            return len(array[0]) * ['NA']
        DEG = [CODON_DEGENERACY_BY_CODON.get(''.join(row), "111") for row in array]
        if '111' in DEG:
            #print('Array', array)
            #print('AA', AA)
            #print('DEG', DEG)
            return len(array[0]) * ['NA']
        temp = []
        deg = "".join(DEG)
        for idx, a in enumerate("".join(AA)):
            if idx > 2:
                temp[idx % 3].append("_".join([a,deg[idx]]))
            else:
                temp.append(["_".join([a,deg[idx]])])
        result = ["|".join(sorted(set(_temp))) for _temp in temp] 
        return result

def infer_degeneracy(transcriptObjs, samples, variant_arrays_by_seq_id):
    #degeneracy_arrays_by_sample = collections.defaultdict(list)
    warnings = []
    # needs to know how many sites in output
    total_sites = 0
    transcriptObjs_blessed = []
    for transcriptObj in tqdm(transcriptObjs, total=len(transcriptObjs), desc="[%] Checking for ORFs... ", ncols=150, position=0, leave=True):
        if not transcriptObj.is_orf():
            warnings.append("[-] Transcript %s has no ORF: START=%s, STOP=%s, DIVISIBLE_BY_3=%s (will be skipped)" % (transcriptObj.transcript_id, transcriptObj.has_start(), transcriptObj.has_stop(), transcriptObj.is_divisible_by_three()))
        else:
            total_sites += transcriptObj.positions.shape[0]
            transcriptObjs_blessed.append(transcriptObj)
    data = np.zeros(total_sites, dtype={'names':('sequence_id', 'start', 'end', 'degeneracy', 'codon_pos', 'orientation'),
                                        'formats':('U16', 'i8', 'i8', 'U16', 'i1', 'U1')})
    
    degeneracy_array = np.zeros((total_sites, len(samples)), dtype='U16')
    print("\n".join(warnings))
    print(type(degeneracy_array), degeneracy_array.shape, degeneracy_array)
    
    offset = 0 
    #print(type(data), data.shape, data)
    for transcriptObj in tqdm(transcriptObjs_blessed, total=len(transcriptObjs_blessed), desc="[%] Inferring degeneracy... ", ncols=150, position=0, leave=True):
        start, end = offset, offset + transcriptObj.positions.shape[0]
        #print('start', start, 'end', end)
        if not transcriptObj.sequence_id in variant_arrays_by_seq_id:
            cds_pos_mask = np.array([])
        else:
            cds_pos_mask = np.isin(transcriptObj.positions, variant_arrays_by_seq_id[transcriptObj.sequence_id]['POS'], assume_unique=True) # will crash if non-unique pos

        data[start:end]['sequence_id'] = transcriptObj.sequence_id
        data[start:end]['start'] = transcriptObj.positions
        data[start:end]['end'] = transcriptObj.positions + 1
        data[start:end]['codon_pos'][0::3] = 1
        data[start:end]['codon_pos'][1::3] = 2
        data[start:end]['codon_pos'][2::3] = 3
        data[start:end]['orientation'] = transcriptObj.orientation
        if np.any(cds_pos_mask):
            gt_pos_mask = np.isin(variant_arrays_by_seq_id[transcriptObj.sequence_id]['POS'], transcriptObj.positions, assume_unique=True) # will crash if non-unique pos
            alleles = np.column_stack(
                    [
                    variant_arrays_by_seq_id[transcriptObj.sequence_id]['REF'][gt_pos_mask], 
                    variant_arrays_by_seq_id[transcriptObj.sequence_id]['ALT'][gt_pos_mask],
                    np.full(gt_pos_mask[gt_pos_mask==True].shape, '')
                    ])
            for idx, sample in enumerate(samples):
                # at least one variant pos in VCF (even if all variants are HOMREF)
                gt_sample = variant_arrays_by_seq_id[transcriptObj.sequence_id]['GT'][gt_pos_mask, idx]
                idx0 = np.array([np.arange(alleles.shape[0]), gt_sample[:,0]])
                idx1 = np.array([np.arange(alleles.shape[0]), gt_sample[:,1]])
                variants = np.vstack(
                    [
                    alleles[tuple(idx0)], 
                    alleles[tuple(idx1)]
                    ]).T
                cds_sample_array = np.full((transcriptObj.positions.shape[0], variants.shape[1]), '') # has shape (gt_sample, ploidy)
                cds_sample_array[cds_pos_mask] = variants
                cds_sample_array[~cds_pos_mask,0] = transcriptObj.sequence[~cds_pos_mask]
                degeneracies = []
                for i in range(0, len(transcriptObj.sequence), 3):
                    codon_list = list(filter(lambda codon: len (codon) == 3, ["".join(x) for x in itertools.product(*cds_sample_array[i:i+3])]))
                    _degeneracy = degeneracy(codon_list) if codon_list else 3 * ['NA']
                    degeneracies.append(_degeneracy)
                #_deg = np.concatenate(degeneracies)
                degeneracy_array[start:end, idx] = np.concatenate(degeneracies)
                #print("#", transcriptObj.transcript_id, transcriptObj.positions.shape, _deg.shape, _deg[0:3])
                #degeneracy_arrays_by_sample[sample].append(_deg)
        else:
            #sequence_id_arrays.append(np.full(site_count, transcriptObj.sequence_id))
            #start_arrays.append(transcriptObj.positions)
            #end_arrays.append(transcriptObj.positions + 1)
            for idx, sample in enumerate(samples):
                degeneracy_array[start:end, idx] = transcriptObj.degeneracy
                #degeneracy_arrays_by_sample[sample].append(transcriptObj.degeneracy)
        offset = end
    #print(data)
    for idx, sample in tqdm(enumerate(samples), total=len(samples), desc="[%] Writing output... ", ncols=150):
        df = pd.DataFrame(data=[], columns=['sequence_id', 'start', 'end', 'degeneracy', 'codon_pos', 'orientation'])
        #data[:]['degeneracy'] = degeneracy_array[:,idx]
        write_df(pd.DataFrame(data=data, columns=['sequence_id', 'start', 'end', 'degeneracy', 'codon_pos', 'orientation']).sort_values(['sequence_id', 'start'], ascending=[True, True]), out_f="%s.bed" % sample, sep='\t', header=False, status=False)

def get_query_regions(transcriptObjs):
    _query_regions_by_sequence_id = collections.defaultdict(list)
    for transcriptObj in tqdm(transcriptObjs, total=len(transcriptObjs), desc="[%] Inferring regions of interest in VCF file... ", ncols=150, position=0, leave=True):
        _query_regions_by_sequence_id[transcriptObj.sequence_id].append(transcriptObj.positions)
    query_regions_by_sequence_id = {}
    for sequence_id, regions in _query_regions_by_sequence_id.items():
        query_regions_by_sequence_id[sequence_id] = np.concatenate(regions)
    return query_regions_by_sequence_id

def main(params):
    try:
        start_time = timer()
        args = docopt(__doc__)
        #print(args)
        #log = lib.log.get_logger(run_params)
        parameterObj = ParameterObj(params, args)
        sequence_by_id = parse_fasta(parameterObj.fasta_file)
        transcriptObjs = get_transcripts(parameterObj, sequence_by_id)
        query_regions_by_sequence_id = get_query_regions(transcriptObjs)
        samples, variant_arrays_by_seq_id = parse_vcf_file(parameterObj.vcf_file, sequence_by_id, query_regions_by_sequence_id)
        infer_degeneracy(transcriptObjs, samples, variant_arrays_by_seq_id)
        print("[*] Total runtime: %.3fs" % (timer() - start_time))
    except KeyboardInterrupt:
        print("\n[X] Interrupted by user after %s seconds!\n" % (timer() - start_time))
        exit(-1)
