# Changelog

## [1.1.0] - 2023-XX-XX
**Breaking changes**:
- preprocess module has been replaced by a separate bioconda package gimble-preprocess.

## [1.1.0b] - 2023-11-06
**Breaking changes**:
- Removal of parallelisation for `makegrid` and `optimize` since implementation is not viable anymore in python 3.8+.

## [1.1.1a] - 2024-06-24
- Added relational constraints argument to query of gridsearch results, e.g. "Ne_A==Ne_B"
