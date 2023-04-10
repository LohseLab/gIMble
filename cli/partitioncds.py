#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""usage: gimble partitioncds               -f FILE -b FILE -v FILE [-e STR -o STR] [-h|--help] [-V|--version]
                                            
    Options:
        -h --help                                   show this
        -f, --fasta_file FILE                       FASTA file
        -v, --vcf_file FILE                         VCF file (filtered)
        -b, --bed_file FILE                         BED file
        -o, --outprefix STR                         Outprefix [default: gimble]
        -e, --exclude STR                           Sample IDs to exclude : '-e sample_A,sample_B'
        -V, --version                               Print version
"""

'''
[To Do]
- Some GFF3 have stop-codon as part of CDS, some do not 
    - needs to have fallback parsing of stop_codon instances to infer stop_codon presence
- make standalone

'''


from timeit import default_timer as timer
from docopt import docopt
import warnings
import numpy as np
import sys
import tempfile
import allel
from tqdm import tqdm
import collections
import pandas as pd
import zarr
import itertools
import shutil
import pathlib

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

def format_proportion(fraction, precision=2):
    if fraction in set(['-', 'N/A']):
        return fraction
    return "{:.{}f}".format(fraction, precision)

def format_count(count):
    if count in set(['-', 'N/A']):
        return count
    return "%s" % str(format(count, ',d'))

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
                                        format_count(count_cds), 
                                        format_count(count_transcript), 
                                        format_proportion(count_cds / count_transcript)))
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
        self.bed = None

    def add_cds_from_df(self, transcript_df, sequence_by_id):
        if not transcript_df['orientation'].nunique():
            sys.exit("[X] More than one orientation found in CDSs of transcript %s\n%s" % (self.transcript_id, transcript_df))
        self.orientation = list(transcript_df['orientation'].unique())[0]
        if not transcript_df['sequence_id'].nunique():
            sys.exit("[X] More than one sequence_id found in CDSs of transcript %s\n%s" % (self.transcript_id, transcript_df))
        self.sequence_id = list(transcript_df['sequence_id'].unique())[0]
        pos_arrays = []
        cds_list = []
        beds = []
        for sequence_id, start, end, _, score, orientation in transcript_df.values.tolist():
            beds.append("\t".join([sequence_id, str(start), str(end), str(_), str(score), orientation]))
            pos_arrays.append(np.arange(start, end))
            cds = sequence_by_id[sequence_id][start:end]
            if orientation == '-':
                cds_list.insert(0, revcom(cds))
            else:
                cds_list.append(cds)
        self.bed = "\n".join(beds)
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
        return True # essentially means not filtering on stop codons

    def is_divisible_by_three(self):
        if self.sequence.shape[0] % 3 == 0:
            return True
        return False

    def is_orf(self):
        if self.has_start() and self.is_divisible_by_three():
            return True
        return False

    def __str__(self):
        return ">%s [%s:%s-%s(%s)]\n%s" % (self.transcript_id, self.sequence_id, self.start, self.end, self.orientation, self.sequence)

class PartitioncdsParameterObj():
    '''Sanitises command line arguments and stores parameters'''
    def __init__(self, params, args):
        self.fasta_file = self._get_path(args['--fasta_file'], path=True)
        self.vcf_file = str(self._get_path(args['--vcf_file'], path=True))
        self.bed_file = self._get_path(args['--bed_file'], path=True)
        self.outprefix = args['--outprefix']
        self.samples_to_exclude = set(args['--exclude'].split(",")) if args['--exclude'] is not None else set([])
        self.tmp_dir = str(tempfile.mkdtemp(prefix='.tmp_gimble_', dir="."))

    def _get_path(self, infile, path=False):
        if infile is None:
            return None
        _path = pathlib.Path(infile).resolve()
        if not _path.exists():
            sys.exit("[X] File not found: %r" % str(infile))
        if path:
            return _path
        return str(_path)

def parse_vcf_file(parameterObj, sequence_ids, query_regions_by_sequence_id):
    print("[+] Parsing VCF file...")
    zstore = zarr.open(parameterObj.tmp_dir, mode='w')
    variant_counts = []
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sample_key, gt_key, pos_key, ref_key, alt_key = 'samples', 'calldata/GT', 'variants/POS', 'variants/REF', 'variants/ALT'
        samples_vcf = list(allel.read_vcf(parameterObj.vcf_file, fields=[sample_key])[sample_key])
        samples_query = [vcf_sample for vcf_sample in samples_vcf if vcf_sample not in set(parameterObj.samples_to_exclude)] # order as they appear genotypes
        zstore.attrs['samples'] = samples_query
        print("[+] Will query %s samples from VCF: %r" % (len(samples_query), ", ".join(samples_query)))
        for idx, sequence_id in tqdm(enumerate(sequence_ids), total=len(sequence_ids), desc="[%] Reading variants...", ncols=100):
            vcf_data = allel.read_vcf(parameterObj.vcf_file, 
                region=sequence_id, 
                samples=samples_query, 
                fields=[gt_key, pos_key, ref_key, alt_key])
            if not vcf_data is None:
                if sequence_id in query_regions_by_sequence_id:
                    pos_array = np.array(vcf_data[pos_key]) - 1 # port to BED (0-based) coordinates
                    cds_mask = np.isin(pos_array, query_regions_by_sequence_id[sequence_id])
                    pos = pos_array[cds_mask]
                    zstore.create_dataset("seqs/%s/variants/pos" % sequence_id, data=pos)
                    ref = vcf_data[ref_key][cds_mask]
                    zstore.create_dataset("seqs/%s/variants/ref" % sequence_id, data=ref, dtype='str')
                    alt = vcf_data[alt_key][cds_mask]
                    zstore.create_dataset("seqs/%s/variants/alt" % sequence_id, data=alt, dtype='str')
                    gts = vcf_data[gt_key][cds_mask]
                    zstore.create_dataset("seqs/%s/variants/gts" % sequence_id, data=gts)
                    variant_counts.append(pos_array.shape[0])
    print("[+] Parsed %s variants." % sum(variant_counts))
    return zstore

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
        #print("Calculating degeneracy:")
        #print(array)
        AA = [AMINOACID_BY_CODON.get(''.join(row), [])*3 for row in array]
        if not AA or None in AA:
            return len(array[0]) * ['NA']
        DEG = [CODON_DEGENERACY_BY_CODON.get(''.join(row), "111") for row in array]
        if '111' in DEG:
            return len(array[0]) * ['NA']
        temp = []
        deg = "".join(DEG)
        for idx, a in enumerate("".join(AA)):
            if idx > 2:
                temp[idx % 3].append("_".join([a,deg[idx]]))
            else:
                temp.append(["_".join([a,deg[idx]])])
        result = ["|".join(sorted(set(_temp))) for _temp in temp] 
        #print("[+] DEG Result:", result)
        return result

def infer_degeneracy(parameterObj, transcriptObjs, zstore):
    #degeneracy_arrays_by_sample = collections.defaultdict(list)
    warnings = []
    beds_rejected = []
    # needs to know how many sites in output
    total_sites = 0
    transcriptObjs_by_sequence_id = collections.defaultdict(list)
    transcriptObjs_valid = 0
    length_by_sequence_id = collections.Counter()

    for transcriptObj in tqdm(transcriptObjs, total=len(transcriptObjs), desc="[%] Checking for ORFs... ", ncols=150, position=0, leave=True):
        if not transcriptObj.is_orf():
            #warnings.append("[-] Transcript %s has no ORF: START=%s, STOP=%s, DIVISIBLE_BY_3=%s (will be skipped)" % (transcriptObj.transcript_id, transcriptObj.has_start(), transcriptObj.has_stop(), transcriptObj.is_divisible_by_three()))
            warnings.append("[-] Transcript %s has no ORF: START=%s, DIVISIBLE_BY_3=%s (will be skipped)" % (transcriptObj.transcript_id, transcriptObj.has_start(), transcriptObj.is_divisible_by_three()))
            beds_rejected.append(transcriptObj.bed)
        else:
            total_sites += transcriptObj.positions.shape[0]
            transcriptObjs_by_sequence_id[transcriptObj.sequence_id].append(transcriptObj)
            length_by_sequence_id[transcriptObj.sequence_id] += transcriptObj.positions.shape[0]
            transcriptObjs_valid += 1
    samples = zstore.attrs['samples']
    degeneracy_chars = "U%s" % (len(samples) * 6 + len(samples) * 1) # could be improved with ploidy?
    chrom_chars = "U%s" % (max([len(sequence_id) for sequence_id in length_by_sequence_id.keys()]) + 1)
    #data = np.zeros(total_sites, dtype={'names':('sequence_id', 'start', 'end', 'degeneracy', 'codon_pos', 'orientation'),'formats':('U16', 'i8', 'i8', degeneracy_chars, 'i1', 'U1')})
    if warnings:
        with open("%s.cds.rejected_transcripts.bed" % parameterObj.outprefix, 'w') as fh:
            fh.write("\n".join(beds_rejected) + "\n")
        print("\n".join(warnings))
    dfs = []
    with tqdm(total=transcriptObjs_valid, ncols=150, desc="[%] Inferring degeneracy... ", position=0, leave=True) as pbar:
        for sequence_id, transcriptObjs in transcriptObjs_by_sequence_id.items():
            offset = 0 
            data = np.zeros(
                length_by_sequence_id[sequence_id], 
                dtype={'names':('sequence_id', 'start', 'end', 'degeneracy', 'codon_pos', 'orientation'),'formats':(chrom_chars, 'i8', 'i8', degeneracy_chars, 'i1', 'U1')})
            if sequence_id in zstore['seqs']: # sequence has variants
                pos = np.array(zstore["seqs/%s/variants/pos" % sequence_id]) 
                gts = np.array(zstore["seqs/%s/variants/gts" % sequence_id])
                alt = np.array(zstore["seqs/%s/variants/alt" % sequence_id])
                ref = np.array(zstore["seqs/%s/variants/ref" % sequence_id])
                alleles_raw = np.column_stack([ref, alt])
                # initiate boolean mask with False
                mask = np.zeros(alleles_raw.shape, dtype=bool)
            if gts.size:
                acs = allel.GenotypeArray(gts).count_alleles()
                # overwrite with True those alleles_raw that occur in gts 
                mask[:,0:acs.shape[1]] = acs
                alleles = np.where(mask, alleles_raw, '')
            for transcriptObj in transcriptObjs:
                start, end = offset, offset + transcriptObj.positions.shape[0]
                data[start:end]['sequence_id'] = transcriptObj.sequence_id
                data[start:end]['start'] = transcriptObj.positions
                data[start:end]['end'] = transcriptObj.positions + 1
                data[start:end]['codon_pos'][0::3] = 1
                data[start:end]['codon_pos'][1::3] = 2
                data[start:end]['codon_pos'][2::3] = 3
                data[start:end]['orientation'] = transcriptObj.orientation
                if not sequence_id in zstore['seqs']:
                    #print("transcriptObj.degeneracy", type(transcriptObj.degeneracy), transcriptObj.degeneracy.shape)
                    data[start:end]['degeneracy'] = transcriptObj.degeneracy
                else:
                    #pos_in_cds_mask = np.isin(pos, transcriptObj.positions, assume_unique=True) # will crash if non-unique pos
                    cds_in_pos_mask = np.isin(transcriptObj.positions, pos, assume_unique=True) # will crash if non-unique pos
                    for i in range(0, len(transcriptObj.sequence), 3):
                        codon_start = start+i
                        if not np.any(cds_in_pos_mask[codon_start:codon_start+3]):
                            data[codon_start:codon_start+3]['degeneracy'] = transcriptObj.degeneracy[i:i+3]
                        else:
                            codon_list = list(filter(lambda codon: len (codon) == 3, ["".join(x) for x in itertools.product(*alleles[codon_start:codon_start+3])]))
                            data[codon_start:codon_start+3]['degeneracy'] = degeneracy(codon_list) if codon_list else 3 * ['NA']
                offset = end
                pbar.update()
            dfs.append(pd.DataFrame(data=data, columns=['sequence_id', 'start', 'end', 'degeneracy', 'codon_pos', 'orientation']))
    shutil.rmtree(parameterObj.tmp_dir)
    #df = pd.DataFrame(data=data, columns=['sequence_id', 'start', 'end', 'degeneracy', 'codon_pos', 'orientation'])
    #write_df(df.sort_values(['sequence_id', 'start'], ascending=[True, True]), out_f="%s.cds.bed" % (parameterObj.outprefix), sep='\t', header=False, status=False)
    write_df(pd.concat(dfs), out_f="%s.cds.bed" % (parameterObj.outprefix), sep='\t', header=False, status=False)


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
        #log = lib.log.get_logger(run_params)
        parameterObj = PartitioncdsParameterObj(params, args)
        sequence_by_id = parse_fasta(parameterObj.fasta_file)
        transcriptObjs = get_transcripts(parameterObj, sequence_by_id)
        query_regions_by_sequence_id = get_query_regions(transcriptObjs)
        zstore = parse_vcf_file(parameterObj, sequence_by_id, query_regions_by_sequence_id)
        infer_degeneracy(parameterObj, transcriptObjs, zstore)
        print("[*] Total runtime: %.3fs" % (timer() - start_time))
    except KeyboardInterrupt:
        print("\n[X] Interrupted by user after %s seconds!\n" % (timer() - start_time))
        exit(-1)