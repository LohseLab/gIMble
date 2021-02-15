import pytest
import lib.gimble
import allel
import numpy as np
import collections
import tests.aux_functions as af
# exceptions: checking that function receives correct input
# assertions: checking that function operates as expected for given input

# testcases
# - all
# - monomorphic-only
# - multiallelic-only
import sys
np.set_printoptions(threshold=sys.maxsize)
genotypes_unfolded = {
    'gt_1': np.array([
        [[ 0, 0], [ 0, 0]],
        [[ 0, 0], [ 0, 1]],
        [[ 0, 1], [ 0, 0]],
        [[ 0, 1], [ 0, 1]],
        [[ 0, 0], [ 1, 1]],
        [[ 0, 0], [ 2, 2]],
        [[ 1, 1], [ 2, 2]],
        [[ 1, 2], [ 1, 2]],
        [[ 1, 2], [ 2, 2]],
        [[ 2, 2], [ 1, 2]],
        [[-1,-1], [ 0, 0]],
        [[ 0, 1], [ 2, 3]]
        ]),
    'gt_biallelic': af.sim_combi_gt_matrix(p=2, s=2, alleles=[0, 1]),
    'gt_triallelic': af.sim_combi_gt_matrix(p=2, s=2, alleles=[0, 1, 2]),
    'gt_missing': af.sim_combi_gt_matrix(p=2, s=2, alleles=[-1, 0]),
    'gt_monomorphic': af.sim_random_gt_matrix(p=2, s=2, n=10, alleles=[0]),
    'gt_random': af.sim_random_gt_matrix(p=2, s=2, n=100, alleles=[-1,0,1,2]),
}

genotypes_folded = {
    'gt_1': np.array([
        [[ 0, 0], [ 0, 0]],
        [[ 0, 0], [ 0, 1]],
        [[ 0, 1], [ 0, 0]],
        [[ 0, 1], [ 0, 1]],
        [[ 0, 0], [ 1, 1]],
        [[ 0, 0], [ 1, 1]],
        [[ 0, 0], [ 1, 1]],
        [[ 0, 1], [ 0, 1]],
        [[ 0, 1], [ 1, 1]],
        [[ 0, 0], [ 1, 0]],
        [[-1,-1], [ 0, 0]],
        [[ 0, 1], [ 2, 3]]
        ]),
}

#@pytest.fixture(scope='class')
#def config(request):
#    yield request.param

# steps
# 2D genotype array

#@pytest.mark.parametrize('genotypes', ['gt_example', 'gt_biallelic', 'gt_triallelic', 'gt_missing', 'gt_monomorphic'], indirect=True, scope='class')

@pytest.mark.genotypes
class Test_genotypes:
        # pos = np.arange(1, 2*gts.shape[0], 2) # positions of GTs in VCF (only for didactic purpose)
        # sa_genotype_array = allel.GenotypeArray(gts)
        # print("[+] %s variants in VCF file on the following positions:\n%s" % (gts.shape[0], str(list(pos))))
        # block_sites = np.arange(25).reshape(5,5)
        # print("[+] block_sites inferred from BED file: \n%s" % (block_sites))
        # block_sites_variant_bool = np.isin(block_sites, pos, assume_unique=True)
        # print("[+] block_sites in VCF file: \n%s" % (block_sites_variant_bool))
        # return (sa_genotype_array, block_sites_variant_bool, block_sites)

    def test_intervals_to_sites(self):
        starts = np.array([0, 5, 8, 11, 15])
        ends = np.array([2, 8, 9, 13, 18])
        result = np.array([0, 1, 5, 6, 7, 8, 11,12,15,16,17])
        assert np.all(np.equal(lib.gimble.intervals_to_sites((starts, ends)), result))
    
    def test_no_intervals_to_sites(self):
        starts = np.array([0, 5, 8, 11, 15])
        ends = np.array([2, 8, 9, 13, 18])
        assert np.all(np.equal(lib.gimble.intervals_to_sites((starts, None)), None))
        assert np.all(np.equal(lib.gimble.intervals_to_sites((None, ends)), None))
        assert np.all(np.equal(lib.gimble.intervals_to_sites((None, None)), None))

    def test_sites_to_blocks_1(self):
        sites = np.array([0, 1, 5, 6, 7, 8, 11,12,15,16,17])
        block_length = 10
        block_span = 20
        result_exp = np.array([0, 1, 5, 6, 7, 8, 11,12,15,16])
        result_obs = lib.gimble.sites_to_blocks(sites, block_length, block_span, debug=True)
        assert result_obs.ndim == 2 # because np.all(np.equal(np.array([0, 1, 5, 6, 7, 8, 11,12,15,16]), np.array([[0, 1, 5, 6, 7, 8, 11,12,15,16]]))) => True
        assert np.all(np.equal(result_obs, result_exp))

    def test_sites_to_blocks_2(self):
        sites = np.array([0, 1, 5, 6, 7, 8, 11,12,15,16,17])
        block_length = 5
        block_span = 10
        result_exp = np.array([[0, 1, 5, 6, 7], 
                               [8, 11,12,15,16]])
        result_obs = lib.gimble.sites_to_blocks(sites, block_length, block_span, debug=True)
        assert result_obs.ndim == 2
        assert np.all(np.equal(result_obs, result_exp))

    def test_sites_to_blocks_3(self):
        sites = np.array([0, 1, 5, 6, 7, 8, 11,12,15,16,17])
        block_length = 5
        block_span = 8
        result_exp = None
        result_obs = lib.gimble.sites_to_blocks(sites, block_length, block_span, debug=True)
        assert result_obs == result_exp

    def test_sites_to_blocks_4(self):
        sites = None
        block_length = 5
        block_span = 8
        result_exp = None
        result_obs = lib.gimble.sites_to_blocks(sites, block_length, block_span, debug=True)
        assert result_obs == result_exp

    def test_sites_to_blocks_5(self):
        sites = np.array([0, 1, 5, 6, 7, 8, 11,12,15,16,17])
        block_length = 0
        block_span = 10
        result_exp = None
        result_obs = lib.gimble.sites_to_blocks(sites, block_length, block_span, debug=True)
        assert result_obs == result_exp

    def test_sites_to_blocks_6(self):
        sites = af.sim_sites(10000)
        block_length = 50
        block_span = 100
        result_obs = lib.gimble.sites_to_blocks(sites, block_length, block_span, debug=True)
        assert result_obs.ndim == 2
        assert result_obs.shape[1] == block_length 
        assert result_obs.shape[0] > 120 # more than 120 blocks

    def test_blocks_to_arrays(self):
        '''Generates all possible 256 diploid GT combinations out of the 4 alleles [-1, 0, 1, 2]
        Pos of variants is always first base in block 
        '''
        gts = af.sim_combi_gt_matrix(p=2, s=2, alleles=[-1, 0, 1, 2]) 
        block_length, block_count = 5, 256
        blocks = np.arange(block_length * block_count).reshape(-1, block_length)
        pos = blocks[:,0]
        starts_exp = np.arange(0, block_length * block_count, block_length)
        ends_exp = np.arange(block_length, (block_length * block_count) + block_length, block_length)
        multiallelic_exp = np.array([
            [0],[0],[0],[0],[0],[0],[0],[0],[0],[0],[0],[0],[0],[0],[0],[0],
            [0],[0],[0],[0],[0],[0],[0],[0],[0],[0],[0],[1],[0],[0],[1],[0],
            [0],[0],[0],[0],[0],[0],[0],[1],[0],[0],[0],[0],[0],[1],[0],[0],
            [0],[0],[0],[0],[0],[0],[1],[0],[0],[1],[0],[0],[0],[0],[0],[0],
            [0],[0],[0],[0],[0],[0],[0],[0],[0],[0],[0],[1],[0],[0],[1],[0],
            [0],[0],[0],[0],[0],[0],[0],[0],[0],[0],[0],[1],[0],[0],[1],[0],
            [0],[0],[0],[1],[0],[0],[0],[1],[0],[0],[0],[1],[1],[1],[1],[1],
            [0],[0],[1],[0],[0],[0],[1],[0],[1],[1],[1],[1],[0],[0],[1],[0],
            [0],[0],[0],[0],[0],[0],[0],[1],[0],[0],[0],[0],[0],[1],[0],[0],
            [0],[0],[0],[1],[0],[0],[0],[1],[0],[0],[0],[1],[1],[1],[1],[1],
            [0],[0],[0],[0],[0],[0],[0],[1],[0],[0],[0],[0],[0],[1],[0],[0],
            [0],[1],[0],[0],[1],[1],[1],[1],[0],[1],[0],[0],[0],[1],[0],[0],
            [0],[0],[0],[0],[0],[0],[1],[0],[0],[1],[0],[0],[0],[0],[0],[0],
            [0],[0],[1],[0],[0],[0],[1],[0],[1],[1],[1],[1],[0],[0],[1],[0],
            [0],[1],[0],[0],[1],[1],[1],[1],[0],[1],[0],[0],[0],[1],[0],[0],
            [0],[0],[0],[0],[0],[0],[1],[0],[0],[1],[0],[0],[0],[0],[0],[0]]) 
        missing_exp = np.array([
            [1],[1],[1],[1],[1],[1],[1],[1],[1],[1],[1],[1],[1],[1],[1],[1],
            [1],[1],[1],[1],[1],[1],[1],[1],[1],[1],[1],[0],[1],[1],[0],[1],
            [1],[1],[1],[1],[1],[1],[1],[0],[1],[1],[1],[1],[1],[0],[1],[1],
            [1],[1],[1],[1],[1],[1],[0],[1],[1],[0],[1],[1],[1],[1],[1],[1],
            [1],[1],[1],[1],[1],[1],[1],[1],[1],[1],[1],[0],[1],[1],[0],[1],
            [1],[1],[1],[1],[1],[0],[0],[0],[1],[0],[0],[0],[1],[0],[0],[0],
            [1],[1],[1],[0],[1],[0],[0],[0],[1],[0],[0],[0],[0],[0],[0],[0],
            [1],[1],[0],[1],[1],[0],[0],[0],[0],[0],[0],[0],[1],[0],[0],[0],
            [1],[1],[1],[1],[1],[1],[1],[0],[1],[1],[1],[1],[1],[0],[1],[1],
            [1],[1],[1],[0],[1],[0],[0],[0],[1],[0],[0],[0],[0],[0],[0],[0],
            [1],[1],[1],[1],[1],[0],[0],[0],[1],[0],[0],[0],[1],[0],[0],[0],
            [1],[0],[1],[1],[0],[0],[0],[0],[1],[0],[0],[0],[1],[0],[0],[0],
            [1],[1],[1],[1],[1],[1],[0],[1],[1],[0],[1],[1],[1],[1],[1],[1],
            [1],[1],[0],[1],[1],[0],[0],[0],[0],[0],[0],[0],[1],[0],[0],[0],
            [1],[0],[1],[1],[0],[0],[0],[0],[1],[0],[0],[0],[1],[0],[0],[0],
            [1],[1],[1],[1],[1],[0],[0],[0],[1],[0],[0],[0],[1],[0],[0],[0]])
        monomorphic_exp = np.array([
            [4],[4],[4],[4],[4],[4],[4],[4],[4],[4],[4],[4],[4],[4],[4],[4],
            [4],[4],[4],[4],[4],[4],[4],[4],[4],[4],[4],[4],[4],[4],[4],[4],
            [4],[4],[4],[4],[4],[4],[4],[4],[4],[4],[4],[4],[4],[4],[4],[4],
            [4],[4],[4],[4],[4],[4],[4],[4],[4],[4],[4],[4],[4],[4],[4],[4],
            [4],[4],[4],[4],[4],[4],[4],[4],[4],[4],[4],[4],[4],[4],[4],[4],
            [4],[4],[4],[4],[4],[5],[4],[4],[4],[4],[4],[4],[4],[4],[4],[4],
            [4],[4],[4],[4],[4],[4],[4],[4],[4],[4],[4],[4],[4],[4],[4],[4],
            [4],[4],[4],[4],[4],[4],[4],[4],[4],[4],[4],[4],[4],[4],[4],[4],
            [4],[4],[4],[4],[4],[4],[4],[4],[4],[4],[4],[4],[4],[4],[4],[4],
            [4],[4],[4],[4],[4],[4],[4],[4],[4],[4],[4],[4],[4],[4],[4],[4],
            [4],[4],[4],[4],[4],[4],[4],[4],[4],[4],[5],[4],[4],[4],[4],[4],
            [4],[4],[4],[4],[4],[4],[4],[4],[4],[4],[4],[4],[4],[4],[4],[4],
            [4],[4],[4],[4],[4],[4],[4],[4],[4],[4],[4],[4],[4],[4],[4],[4],
            [4],[4],[4],[4],[4],[4],[4],[4],[4],[4],[4],[4],[4],[4],[4],[4],
            [4],[4],[4],[4],[4],[4],[4],[4],[4],[4],[4],[4],[4],[4],[4],[4],
            [4],[4],[4],[4],[4],[4],[4],[4],[4],[4],[4],[4],[4],[4],[4],[5]])
        variation_exp = np.array([
            [0,0,0,0], [0,0,0,0], [0,0,0,0], [0,0,0,0], [0,0,0,0], [0,0,0,0], [0,0,0,0], [0,0,0,0],
            [0,0,0,0], [0,0,0,0], [0,0,0,0], [0,0,0,0], [0,0,0,0], [0,0,0,0], [0,0,0,0], [0,0,0,0], 
            [0,0,0,0], [0,0,0,0], [0,0,0,0], [0,0,0,0], [0,0,0,0], [0,0,0,0], [0,0,0,0], [0,0,0,0], 
            [0,0,0,0], [0,0,0,0], [0,0,0,0], [0,0,0,0], [0,0,0,0], [0,0,0,0], [0,0,0,0], [0,0,0,0], 
            [0,0,0,0], [0,0,0,0], [0,0,0,0], [0,0,0,0], [0,0,0,0], [0,0,0,0], [0,0,0,0], [0,0,0,0], 
            [0,0,0,0], [0,0,0,0], [0,0,0,0], [0,0,0,0], [0,0,0,0], [0,0,0,0], [0,0,0,0], [0,0,0,0], 
            [0,0,0,0], [0,0,0,0], [0,0,0,0], [0,0,0,0], [0,0,0,0], [0,0,0,0], [0,0,0,0], [0,0,0,0], 
            [0,0,0,0], [0,0,0,0], [0,0,0,0], [0,0,0,0], [0,0,0,0], [0,0,0,0], [0,0,0,0], [0,0,0,0], 
            [0,0,0,0], [0,0,0,0], [0,0,0,0], [0,0,0,0], [0,0,0,0], [0,0,0,0], [0,0,0,0], [0,0,0,0], 
            [0,0,0,0], [0,0,0,0], [0,0,0,0], [0,0,0,0], [0,0,0,0], [0,0,0,0], [0,0,0,0], [0,0,0,0], 
            [0,0,0,0], [0,0,0,0], [0,0,0,0], [0,0,0,0], [0,0,0,0], [0,0,0,0], [1,0,0,0], [1,0,0,0], 
            [0,0,0,0], [1,0,0,0], [0,0,0,1], [0,0,0,0], [0,0,0,0], [1,0,0,0], [0,0,0,0], [0,0,0,1], 
            [0,0,0,0], [0,0,0,0], [0,0,0,0], [0,0,0,0], [0,0,0,0], [0,1,0,0], [0,0,1,0], [0,0,0,0], 
            [0,0,0,0], [0,0,1,0], [0,1,0,0], [0,0,0,0], [0,0,0,0], [0,0,0,0], [0,0,0,0], [0,0,0,0], 
            [0,0,0,0], [0,0,0,0], [0,0,0,0], [0,0,0,0], [0,0,0,0], [0,1,0,0], [0,0,0,0], [0,0,1,0], 
            [0,0,0,0], [0,0,0,0], [0,0,0,0], [0,0,0,0], [0,0,0,0], [0,0,1,0], [0,0,0,0], [0,1,0,0], 
            [0,0,0,0], [0,0,0,0], [0,0,0,0], [0,0,0,0], [0,0,0,0], [0,0,0,0], [0,0,0,0], [0,0,0,0], 
            [0,0,0,0], [0,0,0,0], [0,0,0,0], [0,0,0,0], [0,0,0,0], [0,0,0,0], [0,0,0,0], [0,0,0,0], 
            [0,0,0,0], [0,0,0,0], [0,0,0,0], [0,0,0,0], [0,0,0,0], [0,1,0,0], [0,0,1,0], [0,0,0,0], 
            [0,0,0,0], [0,0,1,0], [0,1,0,0], [0,0,0,0], [0,0,0,0], [0,0,0,0], [0,0,0,0], [0,0,0,0], 
            [0,0,0,0], [0,0,0,0], [0,0,0,0], [0,0,0,0], [0,0,0,0], [0,0,0,1], [1,0,0,0], [0,0,0,0], 
            [0,0,0,0], [1,0,0,0], [0,0,0,0], [1,0,0,0], [0,0,0,0], [0,0,0,0], [1,0,0,0], [0,0,0,1], 
            [0,0,0,0], [0,0,0,0], [0,0,0,0], [0,0,0,0], [0,0,0,0], [0,0,0,0], [0,0,0,0], [0,0,0,0], 
            [0,0,0,0], [0,0,0,0], [0,1,0,0], [0,0,1,0], [0,0,0,0], [0,0,0,0], [0,0,1,0], [0,1,0,0], 
            [0,0,0,0], [0,0,0,0], [0,0,0,0], [0,0,0,0], [0,0,0,0], [0,0,0,0], [0,0,0,0], [0,0,0,0], 
            [0,0,0,0], [0,0,0,0], [0,0,0,0], [0,0,0,0], [0,0,0,0], [0,0,0,0], [0,0,0,0], [0,0,0,0], 
            [0,0,0,0], [0,0,0,0], [0,0,0,0], [0,0,0,0], [0,0,0,0], [0,1,0,0], [0,0,0,0], [0,0,1,0], 
            [0,0,0,0], [0,0,0,0], [0,0,0,0], [0,0,0,0], [0,0,0,0], [0,0,1,0], [0,0,0,0], [0,1,0,0], 
            [0,0,0,0], [0,0,0,0], [0,0,0,0], [0,0,0,0], [0,0,0,0], [0,0,0,0], [0,0,0,0], [0,0,0,0], 
            [0,0,0,0], [0,0,0,0], [0,1,0,0], [0,0,1,0], [0,0,0,0], [0,0,0,0], [0,0,1,0], [0,1,0,0], 
            [0,0,0,0], [0,0,0,0], [0,0,0,0], [0,0,0,0], [0,0,0,0], [0,0,0,1], [0,0,0,0], [1,0,0,0], 
            [0,0,0,0], [0,0,0,0], [0,0,0,1], [1,0,0,0], [0,0,0,0], [1,0,0,0], [1,0,0,0], [0,0,0,0]])
        starts_obs, ends_obs, multiallelic_obs, missing_obs, monomorphic_obs, variation_obs = lib.gimble.blocks_to_arrays(blocks, gts, pos)
        assert np.all(np.equal(starts_obs, starts_exp))
        assert np.all(np.equal(ends_obs, ends_exp))
        assert np.all(np.equal(multiallelic_obs, multiallelic_exp))
        assert np.all(np.equal(missing_obs, missing_exp))
        assert np.all(np.equal(monomorphic_obs, monomorphic_exp))
        assert np.all(np.equal(variation_obs, variation_exp))

    def test_blocks_to_arrays_2(self):
        '''                                                                    var     Mo  Mi  Mu
        [[0, 0],[0, 0]] [[0, 0],[0, 1]] [[0, 0],[1, 0]] [[0, 0],[1, 1]]     [2,0,0,1]  [1] [0] [0]
        [[0, 1],[0, 0]] [[0, 1],[0, 1]] [[0, 1],[1, 0]] [[0, 1],[1, 1]]     [0,2,2,0]  [0] [0] [0]
        [[1, 0],[0, 0]] [[1, 0],[0, 1]] [[1, 0],[1, 0]] [[1, 0],[1, 1]]     [0,2,2,0]  [0] [0] [0]
        [[1, 1],[0, 0]] [[1, 1],[0, 1]] [[1, 1],[1, 0]] [[1, 1],[1, 1]]     [2,0,0,1]  [1] [0] [0]
        '''
        gts = af.sim_combi_gt_matrix(p=2, s=2, alleles=[0, 1]) 
        block_length, block_count = 4, 4
        blocks = np.arange(block_length * block_count).reshape(-1, block_length)
        pos = np.arange(block_length * block_count)
        starts_exp = np.arange(0, block_length * block_count, block_length)
        ends_exp = np.arange(block_length, (block_length * block_count) + block_length, block_length)
        multiallelic_exp = np.array([[0],[0],[0],[0]])
        missing_exp = np.array([[0],[0],[0],[0]])
        monomorphic_exp = np.array([[1],[0],[0],[1]])
        variation_exp = np.array([[2,0,0,1],[0,2,2,0],[0,2,2,0],[2,0,0,1]])
        starts_obs, ends_obs, multiallelic_obs, missing_obs, monomorphic_obs, variation_obs = lib.gimble.blocks_to_arrays(blocks, gts, pos)
        assert np.all(np.equal(starts_obs, starts_exp))
        assert np.all(np.equal(ends_obs, ends_exp))
        assert np.all(np.equal(multiallelic_obs, multiallelic_exp))
        assert np.all(np.equal(missing_obs, missing_exp))
        assert np.all(np.equal(monomorphic_obs, monomorphic_exp))
        assert np.all(np.equal(variation_obs, variation_exp))