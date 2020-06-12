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
conda install -c conda-forge oyaml zarr scikit-allel pandas numpy tqdm docopt parallel more-itertools networkx scipy sagelib networkx pygraphviz msprime 
```

Usage
-----

```
Usage: ./gimble <module> [<args>...] [-D -V -h]

  [Modules]
    setup                 Setup DataStore
    info                  Print information about DataStore
    blocks                Generate blocks from data in DataStore
    windows               Generate windows from blocks in DataStore
    model                 Build new model
    inference             Make inference [TBI]
    simulate              Simulate [TBI]
    
  [Options]
    -h, --help                         Show this screen.
    -D, --debug                        Print debug information.
    -V, --version                      Show version.
```

Running gimble 
--------------
 
```
# 1. setup
# - will store input data in ZARR DataStore
./gimble setup -v normalised.decomposed.Balance.SnpGap.NonSnp.PASS.vcf.gz -b master.bed -g assembly.genomefile -s samples.csv -o analysis

# 2. blocks
# - infers bSFs for a given block-size
./gimble blocks -z analysis.z -l 64

# 3. windows
# - constructs windows of blocks along the genome
./gimble windows -z analysis.z -w 500 -s 100 -z analysis.z
```

Preparation of input files
--------------------------

### required software (in no particular order)
- bcftools/samtools (http://www.htslib.org/)
- mosdepth (https://github.com/brentp/mosdepth/)
- genomeTools (http://genometools.org/tools/gt_gff3.html)
- bedtools (https://bedtools.readthedocs.io/en/latest/)
- bedops (https://bedops.readthedocs.io/en/latest/content/reference/file-management/conversion/gff2bed.html)
- GNU parallel
- awk/sed/perl

### Necessary files to run gimble (setup/blocks/windows)

1. Genomefile (sequence_id, length) of assembly
    + Only listed sequence IDs will be considered when parsing VCF 
    + Example:
    ```
        chr01   19000101
        chr02    8880001
        chr03    1280001
    ```

2. Sample file (sample_id, population_id)
    + Only listed sample IDs will be considered when parsing VCF
    + Example:
    ```
        sample_1,population_A
        sample_2,population_A
        sample_3,population_B
        sample_4,population_B
    ```

3. VCF file (filtered, bgzip'ed, indexed)
    + Should only list high quality SNPs. 
    + Multiallelicity of SNPs is not a problem.

4. BED file
    + as produced by BEDTOOLS multiinter, based on output of mosdepth's Callable loci
        ```
        chr01   79      93      1       sample_1
        chr01   93      105     2       sample_1,sample_4
        chr01   105     106     3       sample_1,sample_2,sample_4
        chr01   106     111     4       sample_1,sample_2, sample_3, sample_4
        [...]
        ```
    + lists the regions that are accessible to analysis for each of the samples. 
        => Only variants within those regions are considered
        => All sites within those regions that are not variant (listed in VCF), are considered invariant
    + can be intersected with BED regions of certain genomic features (i.e. intergenic regions) ...
    
### 1. Make genomefile from assembly FASTA 

```
samtools faidx assembly.fa
awk -v OFS='\t' {'print $1,$2'} assembly.fa.fai > assembly.genomefile
```

### 2. Make sample file

```
bcftools query -l heliconius.freebayes.normalised.decomposed.Balance.SnpGap.vcf.gz | sort > samples.txt
# manually add population IDs...
```

### 3. Filter VCF file

#### 3.1 Normalise variants 

- See [Tan et. al. 2015](https://academic.oup.com/bioinformatics/article/31/13/2202/196142)

```
vt normalize -r assembly.fasta raw.vcf.gz | bgzip -c > normalised.vcf.gz
```

#### 3.2 Decompose variants into allelic primitives 

- breaks MNPs into SNPs

```
bcftools view normalised.vcf.gz | vcfallelicprimitives --keep-info --keep-geno -t decomposed | sed '/^##/! s/|/\//g' | sed 's/\.:\.:\.:\.:\.:\.:\.:\./\.\/\.:\.:\.:\.:\.:\.:\.:\./g' | bcftools sort -O z > normalised.decomposed.vcf.gz
```

#### 3.3 Filter variants based on read balance/SNPGap/NonSnps

- Read Balance: `-e RPL<1 | RPR<1 | SAF<1 | SAR<1`
- SNPGap: filter SNPs within `$SNPGAP` base pairs of an indel or other other variant type (`INT[:'indel',mnp,bnd,other,overlap]`)
- NonSnps: `TYPE!="snp"`

```
bcftools filter -Oz -s Balance -m+ -e 'RPL<1 | RPR<1 | SAF<1 | SAR<1' normalised.decomposed.vcf.gz | bcftools filter -Oz -m+ -s+ --SnpGap $SNPGAP:indel,other | bcftools filter -Oz -e 'TYPE!="snp"' -s NonSnp -m+ > normalised.decomposed.Balance.SnpGap.NonSnp.vcf.gz
```

#### 3.4 Subset variants that PASS filters 

- These are the variants gimble will analyse

```
bcftools view -O z -f PASS normalised.decomposed.Balance.SnpGap.NonSnp.vcf.gz > normalised.decomposed.Balance.SnpGap.NonSnp.PASS.vcf.gz
bcftools index -t normalised.decomposed.Balance.SnpGap.NonSnp.PASS.vcf.gz 
```

#### 3.5 Generate BED file of variant positions that FAIL'ed filters 

- this is based on the assumption that genotypes at positions of FAIL'ed variants have unknown genotypes
- if this step is not done, the assumption would be that FAIL'ed variants are invariant, i.e. are HOMREF
- the resulting BED file is later intersected with the BED file of CALLABLE regions below

```
bcftools view -H -f SnpGap,NonSNP,Balance normalised.decomposed.Balance.SnpGap.vcf.gz | perl -lane '$pad=0; print($F[0]."\t".($F[1]-1)."\t".(($F[1]-1)+length($F[3]))."\t".$F[6])' > normalised.decomposed.Balance.SnpGap.NonSnp.FAIL.bed
```

### 4. Generate Multintersect BED file of Callable sites across all samples

#### 4.1 Find BAM files used for variant calling

Assumptions:
- BAM files are sorted
- Duplicted read pairs are marked in BAM files (using PicardTools MarkDuplicates)

#### 4.2 Run mosdepth on each BAM file

Assumptions:
- `FILTER_MINCOVs`, `FILTER_MAXCOVs` are known for each BAM file (informed based on coverage profiles), e.g.:
    + `FILTER_MINCOV`: minimum read depth at a site necessary to call variant|invariant (probably uniform across all BAMs)
    + `FILTER_MAXCOV`: maximum read depth at a site necessary to call variant|invariant (probably a function of mean/median-coverage in each BAM)
- See [link](https://github.com/brentp/mosdepth/#callable-regions-example) for details

```
export MOSDEPTH_Q0=NO_COVERAGE; export MOSDEPTH_Q1=LOW_COVERAGE; export MOSDEPTH_Q2=CALLABLE; export MOSDEPTH_Q3=HIGH_COVERAGE
mosdepth -t 8 -n --quantize 0:1:$FILTER_MINCOV_SAMPLE_1:$FILTER_MAXCOV_SAMPLE_1: sample_1.vs.genome.sorted.MD sample_1.vs.genome.sorted.MD.bam
mosdepth -t 8 -n --quantize 0:1:$FILTER_MINCOV_SAMPLE_2:$FILTER_MAXCOV_SAMPLE_2: sample_2.vs.genome.sorted.MD sample_2.vs.genome.sorted.MD.bam
mosdepth -t 8 -n --quantize 0:1:$FILTER_MINCOV_SAMPLE_3:$FILTER_MAXCOV_SAMPLE_3: sample_3.vs.genome.sorted.MD sample_3.vs.genome.sorted.MD.bam
mosdepth -t 8 -n --quantize 0:1:$FILTER_MINCOV_SAMPLE_4:$FILTER_MAXCOV_SAMPLE_4: sample_4.vs.genome.sorted.MD sample_4.vs.genome.sorted.MD.bam
[...]
```

#### 4.3 Subset CALLABLE regions for each BED file

```
# decompress
gunzip *quantized.bed.gz
# Subset CALLABLE regions
parallel -j1 'grep CALLABLE {} > {.}.callable.bed' ::: *quantized.bed
```

#### 4.4 Run BedTools multintersect with all CALLABLE BED files
Caution
- important to specify sample IDs (-names) in the correct order and in the same way they are named in the VCF file.

```
bedtools multiinter -i \
	sample_1.vs.genome.sorted.MD.quantized.callable.bed \
	sample_2.vs.genome.sorted.MD.quantized.callable.bed \
	sample_3.vs.genome.sorted.MD.quantized.callable.bed \
	sample_4.vs.genome.sorted.MD.quantized.callable.bed \
	-names \
	sample_1 \
	sample_2 \
	sample_3 \
	sample_4 \
	| sort -k1,1V -k2,2n -k3,3n > samples.vs.genome.sorted.MD.quantized.callable.bed
```

### 4.5. Subset CALLABLE BED regions with genome features from annotation

- One might want to restrict the regions which are accessible to gimble analysis further, since different genomic regions are subject to different evolutionary forces. 
- For example, one might want to focus on *intergenic regions* which are expected to be affected more by demographic processes than by selection.

For curation/checking of GFF/GFF3/GTF files, we recommend:

- [Genometools](https://github.com/genometools/genometools)' gff3
- [AGAT](https://github.com/NBISweden/AGAT)

Transfer of GFF3 to BED is often easiest with:

- [BedOps](https://bedops.readthedocs.io/en/latest/content/reference/file-management/conversion/gff2bed.html#usage)' gff2bed

For example:

```
# Tidy BRAKER2 GFF3 and add introns using genometools 
gt gff3 -sort -tidy -retainids -fixregionboundaries -addintrons annotation.braker.gff3 > annotation.braker.gt.gff3

# Convert GFF3 to BED file using bedops 
gff2bed < annotation.braker.gt.gff3 > annotation.braker.gt.bed

# Subset genic features from annotation BED 
awk '$8=="gene"' annotation.braker.gt.bed | bedtools merge -i - | sort -k1,1V -k2,2n -k3,3n > annotation.braker.gt.genic.bed

# Generate intergenic BED (complement of genic features)
bedtools complement -i annotation.braker.gt.genic.bed -g assembly.genomefile | sort -k1,1V -k2,2n -k3,3n > annotation.braker.gt.intergenic.bed

# Intersect intergenic BED with CALLABLE BED
bedtools intersect \
    -a samples.vs.genome.sorted.MD.quantized.callable.bed \
    -b annotation.braker.gt.intergenic.bed | sort -k1,1V -k2,2n -k3,3n > samples.intergenic.callable.bed
```

### 4.6. Subtract FAIL BED file (see 3.5.) from `samples.intergenic.callable.bed`

```
bedtools subtract -a samples.intergenic.callable.bed -b normalised.decomposed.Balance.SnpGap.NonSnp.FAIL.bed | sort -k1,1V -k2,2n -k3,3n > master.bed   
```

Documentation
-------------

Further reading, please see (some parts might be outdated):

- [SMBE Speciation Workshop 2019 (1)](https://github.com/DRL/SMBE-SGE-2019/blob/master/Session_3/README.md)
- [SMBE Speciation Workshop 2019 (2)](https://github.com/DRL/SMBE-SGE-2019/blob/master/Session_4/README.md)