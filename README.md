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
    gridsearch            TBE

Options:
    -h, --help                         Show this screen.
    -v, --version                      Show version.
```

Documentation
-------------

For a short introduction into gIMble, please see:
- [SMBE Speciation Workshop 2019 (1)](https://github.com/DRL/SMBE-SGE-2019/blob/master/Session_3/README.md)
- [SMBE Speciation Workshop 2019 (2)](https://github.com/DRL/SMBE-SGE-2019/blob/master/Session_4/README.md)
