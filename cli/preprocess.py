#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""usage: gimble preprocess                 -f FILE -v FILE -b DIR [-g INT -m INT -M INT -t INT -o STR -k] [-h|--help]
                                            
    Options:
        -h --help                                   show this
        -f, --fasta_file FILE                       FASTA file
        -v, --vcf_file FILE                         VCF file (raw)
        -b, --bam_dir FILE                          Directory containing all BAM files
        -g, --snpgap INT                            SnpGap [default: 2]
        -m, --min_depth FILE                        Min read depth [default: 8]
        -M, --max_depth FILE                        Max read depth (as multiple of SD from mean) [default: 2]
        -t, --threads INT                           Threads [default: 1]
        -o, --outprefix STR                         Outprefix [default: gimble]
        -k, --keep_tmp                              Do not delete temporary files [default: False]
"""

from timeit import default_timer as timer
from docopt import docopt
#import lib.gimblelog
from lib.gimble import RunObj
import lib.gimble
import os

import numpy as np
import tempfile
import pathlib
import subprocess
import sys
import shutil
from tqdm import tqdm
import pysam
import pandas as pd
import tabulate
import oyaml as yaml

def print_df(df):
    table = tabulate.tabulate(df, headers=df.columns, tablefmt="orgtbl", floatfmt=".2f")
    border = "[=] %s" % str('=' * (len(table.split("\n")[0]) - 4))
    print(border)
    print(table)
    print(border)

def write_df(df, out_f='', sep='\t', header=True, status=True):
    if header == True:
        df.to_csv(out_f, index=False, sep=sep)
    else:
        df.to_csv(out_f, index=False, sep=sep, header=False)
    if status == True:
        print("[+] \t=> Wrote %r" % str(out_f))

class ParameterObj(RunObj):
    '''Sanitises command line arguments and stores parameters'''
    def __init__(self, params, args):
        super().__init__(params)
        self.fasta_file = self._get_path(args['--fasta_file'], path=True)
        self.vcf_file = self._get_path(args['--vcf_file'], path=True)
        self.bam_dir = self._get_path(args['--bam_dir'], path=True)
        self.bam_files = list(self.bam_dir.glob('*.bam'))
        self.snpgap = self._get_int(args['--snpgap'])
        self.min_depth = self._get_int(args['--min_depth'])
        self.max_depth = self._get_float(args['--max_depth'])
        self.tmp_dir = tempfile.mkdtemp(prefix='.tmp_gimble_', dir=".")
        self.threads = self._get_int(args['--threads'])
        self.outprefix = args['--outprefix']
        self.keep_tmp = args['--keep_tmp']
        self.gimble_genome_file = pathlib.Path('%s.genomefile' % self.outprefix)
        self.gimble_coverage_summary_file = pathlib.Path('%s.coverage_summary.csv' % self.outprefix)
        self.gimble_vcf_file = pathlib.Path('%s.vcf.gz' % self.outprefix)
        self.gimble_sample_file = pathlib.Path("%s.samples.csv" % self.outprefix)
        self.gimble_log_file = pathlib.Path("%s.log.txt" % self.outprefix)
        self.bed_multi_callable_file =  self.tmp_dir / pathlib.Path('%s.multi_callable.bed' % self.outprefix)
        self.coverage_data = []
        self.commands = []

    def clean_up(self):        
        if self.keep_tmp:
            print("[+] Not deleting Temporary folder %r since '--keep_tmp' was specified." % self.tmp_dir)
        else:
            print("[+] Deleting temporary folder %r ..." % self.tmp_dir)
            shutil.rmtree(self.tmp_dir)

    def write_log(self):
        with open(self.gimble_log_file, 'w') as fh:
            fh.write("%s\n" % "\n".join(self.commands))

def run_command(args):
    try:
        call = subprocess.run(args, text=True, capture_output=True, shell=True)
        call.check_returncode()
        return (call.stdout, call.stderr)
    except subprocess.CalledProcessError as cpe:
        if call.stdout:
            sys.exit('[X] The command following command failed with error code %r:\n[X] => %s\n[X] (STDOUT): %r\n[X] (STDERR): %r' % (cpe.returncode, cpe.cmd, cpe.stdout, cpe.stderr))
        sys.exit('[X] The command following command failed with error code %r:\n[X] => %s\n[X] (STDERR): %r' % (cpe.returncode, cpe.cmd, cpe.stderr.rstrip("\n")))

def generate_genome_file(parameterObj):
    print("[+] Parsing FASTA file...")
    cmd = f"""samtools dict {parameterObj.fasta_file} | grep '^@SQ' | cut -f2-3 | perl -lpe 's/SN:|LN://g' > {parameterObj.gimble_genome_file}"""
    _stdout, _stderr = run_command(cmd)
    parameterObj.commands.append(cmd)
    print("[+] \t=> Wrote %r" % str(parameterObj.gimble_genome_file))
    return parameterObj

def get_coverage_data_from_bam(parameterObj):
    print("[+] Parsing BAM files...")
    for bam_file in tqdm(parameterObj.bam_files, desc="[%] Infer coverage levels from BAM files...", ncols=100):
        try:
            with pysam.AlignmentFile(bam_file, require_index=True) as bam:
                readgroup_id = [line.split("\t")[1] for line in str(bam.header).split("\n") if line.startswith("@RG")][0].replace("ID:", "")
                if not readgroup_id:
                    sys.exit("[X] BAM file %r has no readgroup id set.")
                tmp_prefix = pathlib.Path(pathlib.Path(parameterObj.tmp_dir), pathlib.Path(bam_file.stem))
                cmd = f"mosdepth --fast-mode -t {parameterObj.threads} {tmp_prefix} {bam_file}" # --no-per-base
                _stdout, _stderr = run_command(cmd)
                parameterObj.commands.append(cmd)
                bed_f = '%s.per-base.bed.gz' % tmp_prefix
                bed_df = pd.read_csv(bed_f, compression='gzip', sep="\t", names=['sequence_id', 'start', 'end', 'depth'], 
                    dtype={'sequence_id': str, 'start': np.int, 'end': np.int, 'depth': np.int})
                bed_df = bed_df[bed_df['depth'] > 0] # ignore 0-depth sites
                bed_df['length'] = bed_df['end'] - bed_df['start']
                mean = round(np.average(bed_df['depth'], weights=bed_df['length']), 2)
                sd = round(np.sqrt(np.average((bed_df['depth'] - mean)**2, weights=bed_df['length'])), 2)
                parameterObj.coverage_data.append([bam_file.name, readgroup_id, mean, sd, parameterObj.min_depth, int(mean + parameterObj.max_depth * sd)])
        except IOError:
            sys.exit("[X] BAM file %r has no index. Please (sort and) index your BAM files." % bam_file)
    coverage_df = pd.DataFrame(parameterObj.coverage_data, columns=['bam_f', 'read_group_id', 'cov_mean', 'cov_sd', 'depth_min', 'depth_max'])
    print_df(coverage_df)
    write_df(coverage_df, out_f=parameterObj.gimble_coverage_summary_file, sep='\t')
    return parameterObj

def get_multi_callable_f(parameterObj):
    print("[+] Infer callable regions...")
    read_group_ids, bed_callable_fs = [], [] 
    for bam_f, read_group_id, mean, sd, min_depth, max_depth in tqdm(parameterObj.coverage_data, desc="[%] Inferring callable regions...", ncols=100):
        tmp_prefix = pathlib.Path(parameterObj.tmp_dir, pathlib.Path(bam_f).stem)
        bam_path = pathlib.Path(parameterObj.bam_dir, pathlib.Path(bam_f))
        cmd = f"export MOSDEPTH_Q0=NULL; export MOSDEPTH_Q1=LOW; export MOSDEPTH_Q2=CALLABLE; export MOSDEPTH_Q3=HIGH && mosdepth -t {parameterObj.threads} -n --quantize 0:1:{min_depth}:{max_depth}: {tmp_prefix} {bam_path}"
        _stdout, _stderr = run_command(cmd)
        parameterObj.commands.append(cmd)
        bed_path = tmp_prefix.parent / (tmp_prefix.name + '.quantized.bed.gz')
        bed_df = pd.read_csv(bed_path, compression='gzip', sep="\t", names=['sequence_id', 'start', 'end', 'callability'], 
                    dtype={'sequence_id': str, 'start': np.int, 'end': np.int, 'callability': str}).sort_values(['sequence_id', 'start'], ascending=[True, True])
        bed_callable_df = bed_df[bed_df['callability'] == "CALLABLE"] # subset CALLABLE regions
        bed_callable_f = tmp_prefix.parent / (tmp_prefix.name + '.callable.bed')
        write_df(bed_callable_df, out_f=bed_callable_f, sep='\t', header=False, status=False)
        read_group_ids.append(read_group_id)
        bed_callable_fs.append(str(bed_callable_f))
    cmd = f"bedtools multiinter -i {' '.join(bed_callable_fs)} -names {' '.join(read_group_ids)} | cut -f1-5 > {parameterObj.bed_multi_callable_file}"
    _stdout, _stderr = run_command(cmd)
    parameterObj.commands.append(cmd)
    return parameterObj

def process_vcf_f(parameterObj):
    depth_filter_by_sample_id = {read_group_id: (min_depth, max_depth) for (_, read_group_id, _, _, min_depth, max_depth) in parameterObj.coverage_data}
    bam_sample_ids = depth_filter_by_sample_id.keys()
    cmd = f"""bcftools query -l {parameterObj.vcf_file}"""
    _stdout, _stderr = run_command(cmd)
    vcf_sample_ids = [sample_id for sample_id in _stdout.split("\n") if sample_id]
    sample_id_intersection = set(vcf_sample_ids).intersection(bam_sample_ids)
    if not sample_id_intersection:
        sys.exit('[X] Sample ID mismatch\n[X] BAM files : %s\n[X] VCF files : %s\n' % (bam_sample_ids, vcf_sample_ids))
    sample_id_dict = {sample_id: "" for sample_id in sorted(sample_id_intersection)} 
    sample_df = pd.DataFrame(sorted(sample_id_intersection))
    write_df(sample_df, out_f=parameterObj.gimble_sample_file, sep=',', header=False)
    print("[+] Process variants (this might take a while)...")
    vcf_tmp = pathlib.Path(parameterObj.tmp_dir, '%s.filtered.vcf.gz' % parameterObj.outprefix)
    print("[+] Filtering VCF file ...")
    cmd = f"""bcftools norm --threads {parameterObj.threads} -Ov -f {parameterObj.fasta_file} {parameterObj.vcf_file} | \
            vcfallelicprimitives --keep-info --keep-geno -t decomposed | \
            bcftools plugin fill-AN-AC --threads {parameterObj.threads} -Oz | \
            bcftools filter --threads {parameterObj.threads} -Oz -s Balance -m+ -e 'RPL<1 | RPR<1 | SAF<1 | SAR<1' | \
            bcftools filter --threads {parameterObj.threads} -Oz -m+ -s+ --SnpGap 2 | \
            bcftools filter --threads {parameterObj.threads} -Oz -e 'TYPE!="snp"' -s NonSnp -m+ > {vcf_tmp}"""
    # other alternative is fixing AC/AN with vcffixup (hopefully bcftools is faster/better)
    _stdout, _stderr = run_command(cmd)
    parameterObj.commands.append(cmd)
    print("[+] Subsetting VCF file ...")
    bed_fail_f = pathlib.Path(parameterObj.tmp_dir, '%s.filtered.fail.bed' % parameterObj.outprefix)
    cmd = f"""bcftools view --threads {parameterObj.threads} -H -i "%FILTER!='PASS'" {vcf_tmp} | \
            perl -lane '$pad=0; print($F[0]."\t".($F[1]-1)."\t".(($F[1]-1)+length($F[3]))."\t".$F[6])' | bedtools merge -i - > {bed_fail_f}"""
    _stdout, _stderr = run_command(cmd)
    parameterObj.commands.append(cmd)
    vcf_tmp_pass = pathlib.Path(parameterObj.tmp_dir, '%s.filtered.pass.vcf.gz' % parameterObj.outprefix)
    cmd = f"""bcftools view --threads {parameterObj.threads} -Oz -f "PASS" {vcf_tmp} > {vcf_tmp_pass}"""
    _stdout, _stderr = run_command(cmd)
    print("[+] Adjusting callable regions (Filter FAILs are removed)...")
    gimble_bed_f = '%s.bed' % parameterObj.outprefix
    cmd = f"bedtools subtract -a {parameterObj.bed_multi_callable_file} -b {bed_fail_f} > {gimble_bed_f}"
    _stdout, _stderr = run_command(cmd)
    parameterObj.commands.append(cmd)
    print("[+]\t => Wrote %r ..." % str(gimble_bed_f))
    print("[+] Adjust genotypes in uncallable regions...")
    filter_string = " | ".join(["(FMT/DP[%s]<%s | FMT/DP[%s]>%s)" % (
        idx, depth_filter_by_sample_id[vcf_sample_id][0], 
        idx, depth_filter_by_sample_id[vcf_sample_id][1]) for idx, vcf_sample_id in enumerate(vcf_sample_ids)])
    cmd = f"""bcftools filter --threads {parameterObj.threads} -Oz -S . -e '{filter_string}' {vcf_tmp_pass} > {parameterObj.gimble_vcf_file} && bcftools index -t {parameterObj.gimble_vcf_file}"""
    _stdout, _stderr = run_command(cmd)
    parameterObj.commands.append(cmd)
    print("[+]\t => Wrote %r ..." % str(parameterObj.gimble_vcf_file))
    return parameterObj

def preprocess(parameterObj):
    print("[+] Running 'gimble preprocess'")
    parameterObj = generate_genome_file(parameterObj)
    parameterObj = get_coverage_data_from_bam(parameterObj)
    parameterObj = get_multi_callable_f(parameterObj)
    parameterObj = process_vcf_f(parameterObj)
    parameterObj.write_log()
    parameterObj.clean_up()

def main(params):
    try:
        start_time = timer()
        args = docopt(__doc__)
        #print(args)
        #log = lib.log.get_logger(run_params)
        parameterObj = ParameterObj(params, args)
        #print(parameterObj.__dict__)
        preprocess(parameterObj)
        print("[*] Total runtime: %.3fs" % (timer() - start_time))
    except KeyboardInterrupt:
        print("\n[X] Interrupted by user after %s seconds!\n" % (timer() - start_time))
        exit(-1)
