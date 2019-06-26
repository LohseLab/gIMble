gIMble
=========

Dependencies (via [conda](https://conda.io/miniconda.html))
-------

```
# clone repository
git clone https://github.com/DRL/gimble.git

# create conda enviroment with dependencies
conda create -n gimble && \
source activate gimble && \
conda install -c conda-forge more-itertools tqdm scipy numpy matplotlib sympy giac networkx psutil pandas docopt pytables tabulate git htop && \
conda install -c bioconda pysam 
```

Usage
-----

```
/gIMble -h
Usage: gIMble <module> [<args>...]

Modules:
    blocks                Makes blocks
    variants              Fetches and analyses variants for blocks
    modify                Modifies/filters blocks/variants
    windows               Constructs windows of blocks
    likelihood            Infer likelihood for data given model
    gridsearch            Perform gridsearch on precomputed grid

Options:
    -h, --help                         Show this screen.
    -v, --version                      Show version.
```

Documentation
-------------

For a short introduction into gIMble, please see:
- [SMBE Speciation Workshop 2019 (1)](https://github.com/DRL/SMBE-SGE-2019/blob/master/Session_3/README.md)
- [SMBE Speciation Workshop 2019 (2)](https://github.com/DRL/SMBE-SGE-2019/blob/master/Session_4/README.md)

Preparation of input files
--------------------------

### required software (in no particular order)
- samtools (http://www.htslib.org/)
- GATK CallableLoci: https://software.broadinstitute.org/gatk/documentation/tooldocs/3.8-0/org_broadinstitute_gatk_tools_walkers_coverage_CallableLoci.php
- genomeTools (http://genometools.org/tools/gt_gff3.html)
- bedtools (https://bedtools.readthedocs.io/en/latest/)
- bedops (https://bedops.readthedocs.io/en/latest/content/reference/file-management/conversion/gff2bed.html)
- awk
- GNU parallel

### 1. Make genomefile from assembly FASTA (SEQID\tLENGTH)
```
samtools faidx assembly.fa
awk -v OFS='\t' {'print $1,$2'} assembly.fa.fai > assembly.genomefile
```

### 2. Generate Multintersect BED file
#### 2.1 Find BAM files used for variant calling
Assumptions
- BAM files are sorted
- Duplicted read pairs are marked in BAM files (using PicardTools MarkDuplicates)

### 2.2 Run GATK callableLoci on each BAM file
Assumptions
- BAM files are called sample_X.vs.assembly.sorted.deduped.bam, where X is [1..10] 
- MINCOV, MAXCOV are know for each BAM file (informed based on coverage profile), e.g.:
```
gatk -T CallableLoci -mmq 1 -minDepth $MINCOV -maxDepth $MAXCOV -R assembly.fa -I sample_1.vs.assembly.sorted.deduped.bam --summary sample_1.vs.assembly.sorted.deduped.callable_loci.summary.txt --format BED -o sample_1.vs.assembly.sorted.deduped.callable_loci.bed -U ALLOW_N_CIGAR_READS
```

#### 2.3 Run BedTools multintersect of all CallableLoci BED files
Caution
- important to specify sample IDs (-names) the same way they are named in the VCF file!!!
```
bedtools multiinter -i \
	sample_1.vs.assembly.sorted.deduped.callable_loci.CALLABLE.bed \
	sample_2.vs.assembly.sorted.deduped.callable_loci.CALLABLE.bed \
	sample_3.vs.assembly.sorted.deduped.callable_loci.CALLABLE.bed \
	sample_4.vs.assembly.sorted.deduped.callable_loci.CALLABLE.bed \
	sample_5.vs.assembly.sorted.deduped.callable_loci.CALLABLE.bed \
	sample_6.vs.assembly.sorted.deduped.callable_loci.CALLABLE.bed \
	sample_7.vs.assembly.sorted.deduped.callable_loci.CALLABLE.bed \
	sample_8.vs.assembly.sorted.deduped.callable_loci.CALLABLE.bed \
	sample_9.vs.assembly.sorted.deduped.callable_loci.CALLABLE.bed \
	sample_10.vs.assembly.sorted.deduped.callable_loci.CALLABLE.bed \
	-names 
	sample_1 \
	sample_2 \
	sample_3 \
	sample_4 \
	sample_5 \
	sample_6 \
	sample_7 \
	sample_8 \
	sample_9 \
	sample_10 \
	> samples.vs.assembly.multiinter.callableLoci.bed
```
#### 2.4 Only keep CALLABLE intervals
```
grep CALLABLE samples.vs.assembly.multiinter.callableLoci.bed > samples.vs.assembly.multiinter.only_callable.bed
```
#### 2.5 sort BED file
```
sort -k1,1V -k2,2n -k3,3n samples.vs.assembly.multiinter.only_callable.bed > samples.vs.assembly.multiinter.only_callable.sorted.bed
```
### 3. Generate genic+non_repeat / intergenic+non_repeat BED files
#### 3.1 Tidy BRAKER2 GFF3 and add introns using genometools 
```
gt gff3 -sort -tidy -retainids -fixregionboundaries -addintrons annotation.braker.gff3 > annotation.braker.gt.gff3
```
#### 3.2 Convert GFF3 to BED file using bedops 
```
gff2bed < annotation.braker.gt.gff3 > annotation.braker.gt.bed
```
#### 3.3 Convert RepeatMasker GFF3 to BED file
```
gff2bed < repeatmasker.gff > repeatmasker.bed
```
#### 3.4 Divide braker BED
```
awk '$8=="exon"' annotation.braker.gt.bed > annotation.braker.gt.exons.bed
awk '$8=="intron"' annotation.braker.gt.bed > annotation.braker.gt.introns.bed
awk '$8=="gene"' annotation.braker.gt.bed > annotation.braker.gt.genes.bed
```
#### 3.5 Merge intervals within all BED files (so that overlapping intervals get put together)
```
parallel -j1 'bedtools merge -i {} > {.}.merged.bed' ::: *braker.gt*.bed 
```
#### 3.6 non-repeat + genic
##### 3.6.1 Complement repeat BED
```
bedtools complement -i repeatmasker.merged.bed -g assembly.genomefile > assembly.non_repeats.bed
```
##### 3.6.2 Intersect with genic
```
bedtools intersect -a annotation.braker.gt.genes.merged.bed -b assembly.non_repeats.bed > assembly.non_repeats.genic.bed
```
#### 3.7 non-repeat + intergenic
##### 3.7.1 create intergenic BED
```
bedtools complement -i annotation.braker.gt.genes.merged.bed -g assembly.genomefile > annotation.braker.gt.intergenic.bed
```
##### 3.7.2 non-repeats + intergenic
```
bedtools intersect -a annotation.braker.gt.intergenic.bed -b assembly.non_repeats.bed > assembly.non_repeats.intergenic.bed
```
#### 3.8 Sort BED files 
```
sort -k1,1V -k2,2n -k3,3n assembly.non_repeats.genic.bed > assembly.non_repeats.genic.sorted.bed
sort -k1,1V -k2,2n -k3,3n assembly.non_repeats.intergenic.bed > assembly.non_repeats.intergenic.sorted.bed
```
