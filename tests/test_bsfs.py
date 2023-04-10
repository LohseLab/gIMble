import pytest
import lib.gimble
import numpy as np
import tests.aux_functions as af

@pytest.mark.bsfs
class Test_bsfs:

    def test_tally_variation_bsfs(self):
        variation_blocks = np.array([[0,0,0,0], [0,0,0,0], [0,0,0,0], 
                                     [0,0,0,1], [0,0,0,1], 
                                     [0,0,1,0],
                                   # [0,0,1,1], # FGV
                                     [0,1,0,0], [0,1,0,0], [0,1,0,0], 
                                     [0,1,0,1], [0,1,0,1], 
                                     [0,1,1,0], 
                                   # [0,1,1,1], # FGV
                                     [1,0,0,0], [1,0,0,0], [1,0,0,0], 
                                     [1,0,0,1], [1,0,0,1], 
                                     [1,0,1,0],  
                                   # [1,0,1,1], # FGV
                                     [1,1,0,0], [1,1,0,0], [1,1,0,0], 
                                     [1,1,0,1], [1,1,0,1], 
                                     [1,1,1,0], 
                                   # [1,1,1,1]] # FGV
                                     ])
        result_exp = np.array([[[[3, 2],[1, 0]],
                                [[3, 2],[1, 0]]],
                               [[[3, 2],[1, 0]],
                                [[3, 2],[1, 0]]]])
        result_obs = lib.gimble.tally_variation(variation_blocks, form='bsfs')
        print('result_exp', result_exp)
        print('result_obs', result_obs)
        assert np.all(np.equal(result_exp, result_obs))

    def test_tally_variation_blocks_bsfs_vs_tally(self):
        variation_blocks = af.sim_block_variation(b=10, kmax=3)
        max_k = np.array([2,2,2,2])
        variation_bsfs = lib.gimble.tally_variation(variation_blocks, form='bsfs', max_k=max_k)
        variation_tally = lib.gimble.tally_variation(variation_blocks, form='tally', max_k=max_k)
        assert np.all(np.equal(af.bsfs_to_2d(variation_bsfs), variation_tally))

    def test_tally_variation_windows_bsfs_vs_tally(self):
        variation_blocks = af.sim_window_variation(b=10, w=10, kmax=3)
        max_k = np.array([2,2,2,2])
        variation_bsfs = lib.gimble.tally_variation(variation_blocks, form='bsfs', max_k=max_k)
        variation_tally = lib.gimble.tally_variation(variation_blocks, form='tally', max_k=max_k)
        assert np.all(np.equal(af.bsfs_to_2d(variation_bsfs), variation_tally))