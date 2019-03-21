gIMble
=========

Dependencies (via [conda](https://conda.io/miniconda.html))
-------

```
# clone repository
git clone https://github.com/DRL/gimble.git

# create conda enviroment with dependencies
conda env create -f=gimble.conda.yaml -n gimble

# Activate blocktools conda environment
conda activate gimble
```

Usage
-----

```
# Two populations of diploid samples + Migration ({} -> [])
./gIMble graph -s "[{'a'}, ['b']]" -p 2 -o output/test -m

# Two populations of diploid samples + Migration ({} -> []) + Exodus ({} -> [])
./gIMble graph -s "[{'a'}, ['b']]" -p 2 -o output/test -m -e

# Two populations of diploid samples + Migration ({} -> []) + RExodus ({} <- [])
./gIMble graph -s "[{'a'}, ['b']]" -p 2 -o output/test -m -E

```
