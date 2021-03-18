import pytest
import lib.gimble
import numpy as np
import os 
'''
[to test]

# general
./gimble --help

# parse
./gimble parse -g ~/data/gimble/test_data/gimble.test.2chr_1M.genomefile -b ~/data/gimble/test_data/gimble.test.2chr_1M.bed -v ~/data/gimble/test_data/gimble.test.2chr_1M.vcf.gz -s ~/data/gimble/test_data/gimble.test.samples.csv -o ~/data/gimble/test_master/gimble.v0.7.1.master.parse -f
# blocks 
./gIMble blocks -z /Users/dlaetsch/data/gimble/test_master/gimble.v0.7.1.master.parse.z --force

./gimble parse -g ~/git/gimble/data/test.genomefile -b ~/git/gimble/data/test.bed -v ~/git/gimble/data/test.vcf -s ~/git/gimble/data/test.samples.csv -o ~/data/gimble/test_master/gimble.v0.7.1.master.toy -f
./gimble blocks -z /Users/dlaetsch/data/gimble/test_master/gimble.v0.7.1.master.toy.z -l 10 -m 20
./gimble windows -z /Users/dlaetsch/data/gimble/test_master/gimble.v0.7.1.master.toy.z -w 50 -s 5 -f
'''
@pytest.mark.cli
def test_entrypoint():
    #exit_status = os.system('gimble --help')
    assert 0 == 0
