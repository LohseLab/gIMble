Installation
================

clone repository::
    >>> git clone https://github.com/DRL/gimble.git

create conda enviroment with dependencies::
    >>> conda create -n gimble
    >>> source activate gimble && \
    >>> conda install bedtools bcftools samtools vcflib mosdepth pysam numpy docopt tqdm pandas tabulate zarr nlopt scikit-allel parallel more-itertools networkx giac sagelib matplotlib msprime networkx pygraphviz nlopt sympy cerberus maxima -c conda-forge -c bioconda