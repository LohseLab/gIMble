gimble
=========

[![DOI](https://zenodo.org/badge/176883840.svg)](https://zenodo.org/badge/latestdoi/176883840)

# Installation
```
0. clone repository
>>> git clone https://github.com/DRL/gimble.git

1. Install miniconda from https://conda.io/miniconda.html

2. Create the following conda environment 
>>> conda create -n gimble python=3.7.12 bedtools bcftools samtools vcflib mosdepth pysam numpy docopt tqdm pandas tabulate zarr scikit-allel parallel matplotlib msprime demes dask numcodecs python-newick nlopt -c conda-forge -c bioconda -y

3. Load the environment (needs to be activated when using gimble)
>>> conda activate gimble

4. Install agemo (make sure you have the conda environment activated)
>>> (gimble) pip install agemo

5. Rock'n'roll ...
>>> (gimble) ./gimble --help
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
## 0. `gimble preprocess`

A. generates **genome file** (sequence_id, length) based on FASTA file

B. generates **sample file** (sample_id) based on ReadGroupIDs in BAM files

C. generates **coverage threshold report** for each BAM file

D. processes **VCF file** (best to use a freebayes VCF)

+ decomposition of MNPs into SNPs
+ `{RAW_VARIANTS}` = all variants in VCF
+ `{NONSNP}`: non-SNP variants 
+ `{SNPGAP}`: all variants within +/- X b of {NONSNP} variants
+ `{BALANCE}`: all variabts with any un-balanced allele observation (`-e 'RPL<1 | RPR<1 | SAF<1 | SAR<1'`) 
+ `{FAIL} = {{NONSNP} U {SNPGAP} U {BALANCE}} 
+ `{VARIANTS} = {RAW_VARIANTS} - {FAIL}`
+ sample genotypes in `{VARIANTS}` with read depths outside of coverage thresholds are set to missing (`./.`)
    
E. processes **BAM files**:

+ `{RAW_INVARIANT}` = union of all sites with read depths within coverage thresholds in their respective sample (bedtools multiinter)
+ `{INVARIANTS} = {SITES} - {RAW_INVARIANT}`
    
F. logs all executed commands

```
./gimble preprocess -f FASTA -b BAM_DIR/ -v RAW.vcf.gz -k

# output files 
`gimble.samples.csv`            # A)
`gimble.genomefile`             # B)
`gimble.coverage_summary.csv`   # C)
`gimble.vcf.gz`                 # D)
`gimble.bed`                    # E)
`gimble.log.txt`                # F)
```

### 0.1 Modify input files (manually)
--------------
+ `gimble.genomefile`:
    + [OPTIONAL] remove sequence IDs to ignore them in the analyses
+ `gimble.samples.csv` 
    + [REQUIRED] add population IDs to the sample IDs (must be exactly 2)
    + [OPTIONAL] remove sample IDs to ignore them in the analyses
+ `gimble.bed`
    + [RECOMMENDED] intersect with BED regions of interest to analyse particular genomic regions e.g:  
    ```
    bedtools intersect -a gimble.bed -b my_intergenic_regions.bed > gimble.intergenic.bed
    ``` 

## 1. `gimble parse`
--------------
+ will parse input data into GimbleStore 
```
./gimble parse -v gimble.vcf.gz -b gimble.intergenic.bed -g gimble.genomefile -s gimble.samples.csv -o analysis
```

## 2. `gimble blocks`
--------------
+ infers bSFs for a given block length `'-l'` 
+ block span (`end - start`) can be adjusted (default is `2 * '-l'`)
```
./gimble blocks -z analysis.z -l 64
```

## 3. `gimble windows`
--------------
+ constructs windows of blocks along the genome
```
./gimble windows -z analysis.z -w 500 -s 100 -z analysis.z
```

## 4. `gimble info`
--------------
+ lists basic metrics about the GimbleStore
```
./gimble info -z analysis.z
```

## 5. `gimble tally`
--------------
+ tallies variation in blocks
+ parameter k-max, max count per mutation type beyond which counts are treated 
    as marginals. Order of mutation types is (hetB, hetA, hetAB, fixed) [default: 2,2,2,2]
```
./gimble tally -z analysis.z -k 2,2,2,2 -l blocks_kmax2 -t blocks
./gimble tally -z analysis.z -k 2,2,2,2 -l windows_kmax2 -t windows
```

## 6. `gimble optimize`
------------
+ Searches parameter space for a tally or simulation under a given model

## 7. `gimble makegrid`
------------
+ Computes a grid of parameter combinations for a list of parameters under a given model

## 8. `gimble gridsearch`
--------------
+ Searches a tally or simulation against a grid made with makegrid

## 9. `gimble simulate`
------------
+ Simulates tallies based on fixed parameters under a model
+ Simulates tallies based on results of gridsearches
