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
conda install bedtools bcftools samtools vcflib mosdepth pysam numpy docopt tqdm pandas tabulate oyaml zarr scikit-allel parallel more-itertools networkx sagelib matplotlib msprime networkx pygraphviz -c conda-forge -c bioconda 
```

Usage
-----

```
Usage: gimble <module> [<args>...] [-D -V -h]
  [Modules]
    preprocess            Preprocess input files
    setup                 Setup data store (GStore)
    info                  Print information about GStore [TBI]
    blocks                Generate blocks from data in GStore 
    windows               Generate windows from blocks in GStore (requires blocks)
    model                 Build demographic model
    simulate              Simulate data [TBI] 
    inference             Make inference [TBI] (requires blocks)
    grid                  Make grid [TBI]
    scan                  Scan using grid [TBI] (requires windows)
    
    partitioncds          Partition CDS sites in BED file by degeneracy in sample GTs 
    plotbed               Plot BED file [TBR]

```

Running gimble 
--------------
 
# 0. preprocess

[Steps]
    a generates Genomefile (sequence_id, length) based FASTA file: 
        => `gimble.genomefile`
    b generates Sample file (sample_id) based on ReadGroupIDs in BAM files: 
        => `gimble.samples.csv`
    c generates coverage thresholds (min/max) for each BAM file: 
        => `gimble.coverage_summary.csv`
    d processes VCF file: 
        + decomposition of MNPs into SNPs
        + `{RAW_VARIANTS}` = all variants in VCF
        + `{NONSNP}`: non-SNP variants 
        + `{SNPGAP}`: all variants within +/- X b of {NONSNP} variants
        + `{BALANCE}`: all variabts with any un-balanced allele observation (`-e 'RPL<1 | RPR<1 | SAF<1 | SAR<1'`) 
        + `{FAIL} = {{NONSNP} U {SNPGAP} U {BALANCE}} 
        + `{VARIANTS} = {RAW_VARIANTS} - {FAIL}`
        + sample genotypes in `{VARIANTS}` with read depths outside of coverage thresholds are set to missing (`./.`)
        => `gimble.vcf.gz`
    e processes BAM files
        + {RAW_INVARIANT} = union of all sites with read depths within coverage thresholds in their respective sample (bedtools multiinter)
        + {INVARIANTS} = {SITES} - {RAW_INVARIANT}
        => `gimble.bed`
    f document all executed commands
        => `gimble.log.txt`

[Command]
- Make sure `RAW.vcf.gz` is indexed
`~/gIMble/gIMble preprocess -f FASTA -b BAM_DIR/ -v RAW.vcf.gz -k`

# 1. Adjust input files
a `gimble.genomefile`:
    + [OPTIONAL] remove sequence IDs to ignore them in the analyses
b `gimble.samples.csv` 
    + [REQUIRED] add population IDs to the sample IDs (must be exactly 2)
    + [OPTIONAL] remove sample IDs to ignore them in the analyses
c `gimble.bed`
    + [RECOMMENDED] intersect with BED regions of interest to analyse particular genomic regions

    e.g `bedtools intersect -a gimble.bed -b my_intergenic_regions.bed > gimble.intergenic.bed` 

# 2. setup
+ will extract input data into DataStore 
`./gimble setup -v gimble.vcf.gz -b gimble.intergenic.bed -g gimble.genomefile -s gimble.samples.csv -o analysis`

# 3. blocks
+ infers bSFs for a given block length `'-l'` 
+ block span (`end - start`) can be adjusted (default is `2 * '-l'`)
`./gimble blocks -z analysis.z -l 64`

# 4. windows
+ constructs windows of blocks along the genome
./gimble windows -z analysis.z -w 500 -s 100 -z analysis.z