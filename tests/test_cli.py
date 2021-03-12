import pytest
import lib.gimble
import numpy as np

@pytest.mark.cli

'''
[to test]

# general
./gimble --help

# parse
./gimble parse -g ~/data/gimble/test_data/gimble.test.2chr_1M.genomefile -b ~/data/gimble/test_data/gimble.test.2chr_1M.bed -v ~/data/gimble/test_data/gimble.test.2chr_1M.vcf.gz -s ~/data/gimble/test_data/gimble.test.samples.csv -o ~/data/gimble/test_master/gimble.v0.7.1.master.parse -f
# blocks 
./gIMble blocks -z /Users/dlaetsch/data/gimble/test_master/gimble.v0.7.1.master.parse.z --force
'''

def test_entrypoint():
    exit_status = os.system('gimble --help')
    assert exit_status == 0
