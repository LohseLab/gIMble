#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""usage: gIMble preprocess                 -f FILE -v FILE -b DIR [-g INT -m INT -M INT] [-h|--help]
                                            
    Options:
        -h --help                                   show this
        -f, --fasta_file FILE                       FASTA file
        -v, --vcf_file FILE                         VCF file (raw)
        -b, --bam_dir FILE                          Directory containing all BAM files
        -g, --snpgap INT                            SnpGap [default: 2]
        -m, --min_depth FILE                        Min read depth [default: 8]
        -M, --max_depth FILE                        Max read depth (as multiple of SD from mean) [default: 3]
"""

'''
preprocess 
- arguments
    + BAM_DIR
    + FASTA
    + SNP_GAP [2]
    + MIN_DEPTH [8]
    + MAX_DEPTH_AS_SD_MULTIPLICATOR [3]

- output
    + genome_file
    + samples.csv (not including POPULATION_IDs)
    + filtered VCF file
        * PASS
        * DP-Fails to missing GTs
    + BED file
        * Multiinter CALLABLE w/ filter-FAILs subtracted
- steps
    + 1. Parse FASTA, write genomefile          
        * Validates FASTA
    + 2. Parse BAMs

- BAM
    + mosdepth -> min (8) / max (MEAN+3SD) depth for filtering
    + modepth quantize -> CALLABLE BED
- CALLABLE BEDs
    + BEDTOOLS multiinter    
- VCF
    + VT normalize
    + BCFTOOLS view | vcfallelicprimitives | regex | BCFTOOLS SORT
    + 


# Test Files
## BAMs
parallel -j4 'samtools view -h {} | head -100000 | samtools view -b > ~/heliconius/{.}.10K.bam && samtools index ~/heliconius/{.}.10K.bam' ::: *sorted.MD.bam
## VCF
head -75000 heliconius.freebayes.raw.no_pop_priors.vcf > ~/heliconius/heliconius.freebayes.raw.no_pop_priors.75K.vcf
bgzip -c heliconius.freebayes.raw.no_pop_priors.75K.vcf > heliconius.freebayes.raw.no_pop_priors.75K.vcf.gz && bcftools index -t heliconius.freebayes.raw.no_pop_priors.75K.vcf.gz

mosdepth -b 9_bed_files/hmel2_5.chromosomes.autosomes.intergenic_and_nonrepeats.bed -t 10 chi.CAM25091.f 6_mapping/chi.CAM25091.f.vs.hmel2_5.sorted.MD.bam
mosdepth -b 9_bed_files/hmel2_5.chromosomes.autosomes.intergenic_and_nonrepeats.bed -t 10 chi.CAM25137.f 6_mapping/chi.CAM25137.f.vs.hmel2_5.sorted.MD.bam
mosdepth -b 9_bed_files/hmel2_5.chromosomes.autosomes.intergenic_and_nonrepeats.bed -t 10 chi.CAM580.m 6_mapping/chi.CAM580.m.vs.hmel2_5.sorted.MD.bam
mosdepth -b 9_bed_files/hmel2_5.chromosomes.autosomes.intergenic_and_nonrepeats.bed -t 10 chi.CAM582.m 6_mapping/chi.CAM582.m.vs.hmel2_5.sorted.MD.bam
mosdepth -b 9_bed_files/hmel2_5.chromosomes.autosomes.intergenic_and_nonrepeats.bed -t 10 chi.CAM585.m 6_mapping/chi.CAM585.m.vs.hmel2_5.sorted.MD.bam
mosdepth -b 9_bed_files/hmel2_5.chromosomes.autosomes.intergenic_and_nonrepeats.bed -t 10 chi.CAM586.m 6_mapping/chi.CAM586.m.vs.hmel2_5.sorted.MD.bam
mosdepth -b 9_bed_files/hmel2_5.chromosomes.autosomes.intergenic_and_nonrepeats.bed -t 10 chi.CJ553.m 6_mapping/chi.CJ553.m.vs.hmel2_5.sorted.MD.bam
mosdepth -b 9_bed_files/hmel2_5.chromosomes.autosomes.intergenic_and_nonrepeats.bed -t 10 chi.CJ560.m 6_mapping/chi.CJ560.m.vs.hmel2_5.sorted.MD.bam
mosdepth -b 9_bed_files/hmel2_5.chromosomes.autosomes.intergenic_and_nonrepeats.bed -t 10 chi.CJ564.m 6_mapping/chi.CJ564.m.vs.hmel2_5.sorted.MD.bam
mosdepth -b 9_bed_files/hmel2_5.chromosomes.autosomes.intergenic_and_nonrepeats.bed -t 10 chi.CJ565.m 6_mapping/chi.CJ565.m.vs.hmel2_5.sorted.MD.bam
mosdepth -b 9_bed_files/hmel2_5.chromosomes.autosomes.intergenic_and_nonrepeats.bed -t 10 ros.CAM1841.m 6_mapping/ros.CAM1841.m.vs.hmel2_5.sorted.MD.bam
mosdepth -b 9_bed_files/hmel2_5.chromosomes.autosomes.intergenic_and_nonrepeats.bed -t 10 ros.CAM1880.m 6_mapping/ros.CAM1880.m.vs.hmel2_5.sorted.MD.bam
mosdepth -b 9_bed_files/hmel2_5.chromosomes.autosomes.intergenic_and_nonrepeats.bed -t 10 ros.CAM2045.m 6_mapping/ros.CAM2045.m.vs.hmel2_5.sorted.MD.bam
mosdepth -b 9_bed_files/hmel2_5.chromosomes.autosomes.intergenic_and_nonrepeats.bed -t 10 ros.CAM2059.m 6_mapping/ros.CAM2059.m.vs.hmel2_5.sorted.MD.bam
mosdepth -b 9_bed_files/hmel2_5.chromosomes.autosomes.intergenic_and_nonrepeats.bed -t 10 ros.CAM2519.m 6_mapping/ros.CAM2519.m.vs.hmel2_5.sorted.MD.bam
mosdepth -b 9_bed_files/hmel2_5.chromosomes.autosomes.intergenic_and_nonrepeats.bed -t 10 ros.CAM2552.m 6_mapping/ros.CAM2552.m.vs.hmel2_5.sorted.MD.bam
mosdepth -b 9_bed_files/hmel2_5.chromosomes.autosomes.intergenic_and_nonrepeats.bed -t 10 ros.CJ2071.m 6_mapping/ros.CJ2071.m.vs.hmel2_5.sorted.MD.bam
mosdepth -b 9_bed_files/hmel2_5.chromosomes.autosomes.intergenic_and_nonrepeats.bed -t 10 ros.CJ531.m 6_mapping/ros.CJ531.m.vs.hmel2_5.sorted.MD.bam
mosdepth -b 9_bed_files/hmel2_5.chromosomes.autosomes.intergenic_and_nonrepeats.bed -t 10 ros.CJ533.m 6_mapping/ros.CJ533.m.vs.hmel2_5.sorted.MD.bam
mosdepth -b 9_bed_files/hmel2_5.chromosomes.autosomes.intergenic_and_nonrepeats.bed -t 10 ros.CJ546.m 6_mapping/ros.CJ546.m.vs.hmel2_5.sorted.MD.bam

'''

from timeit import default_timer as timer
from docopt import docopt
#import lib.gimblelog
from lib.gimble import RunObj
import lib.gimble
import os

import tempfile
import pathlib
import subprocess
import sys
import shutil
from tqdm import tqdm
import pysam
import re

class ParameterObj(RunObj):
    '''Sanitises command line arguments and stores parameters'''
    def __init__(self, params, args):
        super().__init__(params)
        self.fasta_file = self._get_path(args['--fasta_file'])
        self.vcf_file = self._get_path(args['--vcf_file'])
        self.bam_files = list(self._get_path(args['--bam_dir']).glob('*.bam'))
        self.snpgap = self._get_int(args['--snpgap'])
        self.min_depth = self._get_int(args['--min_depth'])
        self.max_depth = self._get_int(args['--max_depth'])
        self.tmp_dir = tempfile.mkdtemp()

    def clean_up(self):
        shutil.rmtree(self.tmp_dir)

def run_command(args):
    try:
        call = subprocess.run(args, text=True, capture_output=True, shell=True)
        call.check_returncode()
        return (call.stdout, call.stderr)
    except subprocess.CalledProcessError as cpe:
        sys.exit('[X] The command %r failed with error code %r\n[X] STDOUT: %r\n[X] STDERR: %r' % (" ".join(cpe.cmd), cpe.returncode, cpe.stdout, cpe.stderr))

def generate_genome_file(parameterObj):
    tmp_fai_file = pathlib.Path(pathlib.Path(parameterObj.tmp_dir), pathlib.Path(parameterObj.fasta_file.with_suffix('.fai').name))
    genome_file = pathlib.Path(parameterObj.fasta_file.with_suffix('.genomefile').name)
    fasta_file = parameterObj.fasta_file
    command = [f"samtools faidx {parameterObj.fasta_file} && cut -f1,2 {tmp_fai_file} > {genome_file}"]
    return genome_file

def get_cov_mean_sd_by_readgroup_id(parameterObj):
    RG_ID_REGEX = re.compile('^@RG\tID:\w+\t\s+')
    for bam_file in tqdm(parameterObj.bam_files, desc="[%] Infer coverage levels from BAM files...", ncols=100):
        try:
            with pysam.AlignmentFile(bam_file, require_index=True) as bam:
                if not bam.has_tag('RG'):
                    sys.exit("[X] BAM file %r has no READGROUP set.")
                readgroup_id = bam.get_tag('RG')
                tmp_prefix = pathlib.Path(pathlib.Path(tmp_dir), pathlib.Path(bam_file.stem))
                mosdepth_call = [f"mosdepth --fast-mode -t {parameterObj.threads} {tmp_prefix} {bam_file}"] # --no-per-base
                # generates 
                # {tmp_prefix}.per-base.bed.gz          => Relevant for coverage 
                # {tmp_prefix}.mosdepth.global.dist.txt => Relevant for plotting 
                _stdout, _stderr = run_command(mosdepth_call)
                
                bed_df = pd.read_csv(tmp_prefix.with_suffix('per-base.bed.gz'), compression='gzip')
                print(bed_df)
        except IOError:
            sys.exit("[X] BAM file %r has no index.")
    return cov_mean_sd_by_readgroup_id

def preprocess(parameterObj):
    genome_file = generate_genome_file(parameterObj)
    cov_mean_sd_by_readgroup_id = get_cov_mean_sd_by_readgroup_id(parameterObj)

def main(params):
    try:
        start_time = timer()
        args = docopt(__doc__)
        #print(args)
        #log = lib.log.get_logger(run_params)
        parameterObj = ParameterObj(params, args)
        print(parameterObj.__dict__)
        preprocess(parameterObj)
        print("[*] Total runtime: %.3fs" % (timer() - start_time))
    except KeyboardInterrupt:
        print("\n[X] Interrupted by user after %s seconds!\n" % (timer() - start_time))
        exit(-1)
