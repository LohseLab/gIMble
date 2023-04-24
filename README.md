gimble
=========

[![DOI](https://zenodo.org/badge/176883840.svg)](https://zenodo.org/badge/latestdoi/176883840)

# Table of contents

- [Installation](#installation)
- [Usage](#usage)
- [Workflow](#workflow)
- [Gimble modules](#gimble-modules)

# Installation
```
# clone repository
>>> git clone https://github.com/DRL/gimble.git
# Install miniconda from https://conda.io/miniconda.html
# ...
# Create the following conda environment 
>>> conda create -n gimble python=3.7.12 bedtools bcftools samtools vcflib mosdepth pysam numpy docopt tqdm pandas tabulate zarr scikit-allel parallel matplotlib msprime demes dask numcodecs python-newick nlopt -c conda-forge -c bioconda -y
# Load the environment (needs to be activated when using gimble)
>>> conda activate gimble
# Install agemo (make sure you have the conda environment activated)
>>> (gimble) pip install agemo
# Start gimble'ing ...
>>> (gimble) gIMble/gimble --help
```

# Usage
```
usage: gimble <module> [<args>...] [-V -h]

  [Input]
    preprocess            Preprocess input files
    parse                 Parse files into GimbleStore
    blocks                Generate blocks from parsed data in GimbleStore (requires 'parse')
    windows               Generate windows from blocks in GimbleStore (requires 'blocks')
    tally                 Tally variation for inference (requires 'blocks' or 'windows')

  [Simulation]
    simulate              Simulate data based on specific parameters or gridsearch results  
    
  [Inference]
    optimize              Perform global parameter optimisation on tally/simulation
    makegrid              Precalculate grid of parameters
    gridsearch            Evaluate tally/simulation against a precomputed grid (requires 'makegrid')

  [Info]
    info                  Print metrics about data in GimbleStore
    list                  List information saved in GimbleStore
    query                 Extract information from GimbleStore
    delete                Delete information in GimbleStore

  [Experimental]
    partitioncds          Partition CDS sites in BED file by degeneracy in sample GTs

  [Options]
    -h, --help            Show this screen
    -V, --version         Show version
```

# Workflow

# Gimble modules
  * [preprocess](#preprocess)
    + [VCF processing details](#vcf-processing-details)
    + [BAM processing details](#bam-processing-details)
    + [Manually modify preprocessed files](#manually-modify-preprocessed-files)
  * [parse](#parse)
  * [blocks](#blocks)
  * [windows](#windows)
  * [info](#info)
  * [tally](#tally)
  * [optimize](#optimize)
  * [makegrid](#makegrid)
  * [gridsearch](#gridsearch)
  * [simulate](#simulate)


## preprocess

The `preprocess` module assures that input files are adequately filtered and processed so that the `gimble` workflow can be completed successfully. 
While this processing of input files could be done more efficiently with other means, it has the advantage of generating a VCF file complies with `gimble` data requirements but which can also be used in alternative downstream analyses.

```
./gimble preprocess -f FASTA -b BAM_DIR/ -v RAW.vcf.gz -k
```

Based on the supplied input files:
- `-f` : FASTA of reference
- `-b`: directory of BAM files, composed of readgroup-labelled reads mapped to reference 
- `-v`: compressed+indexed Freebayes VCF file 

the module produces the following output files:

- **genome file** (sequence_id, length) based on FASTA file
- **sample file** (sample_id) based on ReadGroupIDs in BAM files
- **coverage threshold report** for each BAM file
- **gimble VCF file** (see [VCF processing details](#vcf-processing-details))
- **gimble BED file** (see [BAM processing details](#bam-processing-details))
- log of executed commands

After running, output files require manual user input. See [Manually modify files](#manually-modify-preprocessed-files))

---

### VCF processing details 
- MNPs are decomposed into SNPs
- variant sets are defined as
  ```
  - {RAW_VARIANTS} := all variants in VCF
  - {NONSNP}       := non-SNP variants 
  - {SNPGAP}       := all variants within +/- X b of {NONSNP} variants
  - {QUAL}         := all variants with QUAL below --min_qual
  - {BALANCE}      := all variants with any un-balanced allele observation (-e 'RPL<1 | RPR<1 | SAF<1 | SAR<1') 
  - {FAIL}         := {{NONSNP} U {SNPGAP} U {QUAL} U {BALANCE}} 
  - {VARIANTS}     := {RAW_VARIANTS} - {FAIL}```
The processed VCF file 
  - only contains variants from set `{VARIANTS}` 
  - only contains sample genotypes with sample read depths within coverage thresholds (others are set to missing, i.e. `./.`)

---

### BAM processing details 
- definition of BED region sets
  ```
  - {CALLABLE_SITES}   := for each BAM, regions along which sample read depths lie within coverage thresholds.
  - {CALLABLES}        := bedtools multiintersect of all {CALLABLE_SITES} regions across all BAMs/samples
  - {FAIL_SITES}       := sites excluded due to {FAIL} variant set (during VCF processing)
  - {SITES}            := {CALLABLES} - {FAIL_SITES}
  ```
Resulting BED file 
- only contains BED regions from set `{SITES}` 
- lists which samples are adequately covered along each region

---

#### Manually modify preprocessed files
+ `gimble.genomefile`:
    + **[Optional]** remove sequence IDs to ignore them in the analyses
+ `gimble.samples.csv` 
    + **[Required]** add population IDs the second column. Must be exactly 2 populations
    + **[Optional]** remove sample IDs to ignore them in the analyses 
+ `gimble.bed`
    + **[Recommended]** intersect with BED regions of interest to analyse particular genomic regions, e.g:  
    ```
    bedtools intersect -a gimble.bed -b my_intergenic_regions.bed > gimble.intergenic.bed
    ``` 
---

## parse
+ will parse input data into GimbleStore 
```
./gimble parse -v gimble.vcf.gz -b gimble.intergenic.bed -g gimble.genomefile -s gimble.samples.csv -o analysis
```

## blocks
--------------
+ infers bSFs for a given block length `'-l'` 
+ block span (`end - start`) can be adjusted (default is `2 * '-l'`)
```
./gimble blocks -z analysis.z -l 64
```

## windows
--------------
+ constructs windows of blocks along the genome
```
./gimble windows -z analysis.z -w 500 -s 100 -z analysis.z
```

## info
--------------
+ lists basic metrics about the GimbleStore
```
./gimble info -z analysis.z
```

## tally
--------------
+ tallies variation in blocks
+ parameter k-max, max count per mutation type beyond which counts are treated 
    as marginals. Order of mutation types is (hetB, hetA, hetAB, fixed) [default: 2,2,2,2]
```
./gimble tally -z analysis.z -k 2,2,2,2 -l blocks_kmax2 -t blocks
./gimble tally -z analysis.z -k 2,2,2,2 -l windows_kmax2 -t windows
```

## optimize
------------
+ Searches parameter space for a tally or simulation under a given model

## makegrid
------------
+ Computes a grid of parameter combinations for a list of parameters under a given model

## gridsearch
--------------
+ Searches a tally or simulation against a grid made with makegrid

## simulate
------------
+ Simulates tallies based on fixed parameters under a model
+ Simulates tallies based on results of gridsearches
