import itertools
import numpy as np
import pytest
import lib.math
import lib.gimble
import lib.simulate
import tests.aux_functions as af

from lib.GeneratingFunction.gf import togimble

@pytest.mark.ETPs
@pytest.mark.parametrize('config_file', [
    af.get_test_config_file(task='simulate', model='DIV', label='test_ETPs_1'),
    af.get_test_config_file(task='simulate', model='MIG_BA', label='test_ETPs_2'),
    af.get_test_config_file(task='simulate', model='IM_AB', label='test_ETPs_3')
    ], scope='class')
class Test_ETPs:

    @pytest.fixture(scope='class')
    def gimble_ETPs(self):
        def _gimble_ETPs(config):
            gf = lib.math.config_to_gf(config)
            gfEvaluatorObj = togimble.gfEvaluator(gf, config['max_k'], lib.gimble.MUTYPES, config['gimble']['precision'], exclude=[(2,3),])
            ETPs = lib.math.new_calculate_all_ETPs(
                gfEvaluatorObj, 
                config['parameters_expanded'], 
                config['simulate']['reference_pop_id'], 
                config['simulate']['block_length'], 
                config['mu']['mu'], 
                processes=1, 
                verbose=False
                )
            return ETPs
        return _gimble_ETPs

    def test_ETPs_model(self, gimble_ETPs, config_file, cmd_plot=True):
        config = lib.gimble.load_config(config_file)
        simmed_ETPs = lib.simulate.run_sims(config['demographies'], config['simulate']['recombination_rate'], config, config['replicates'], config['seeds'], 1)
        gimbled_ETPs = gimble_ETPs(config)
        for idx, (gimble_combo, sim_combo) in enumerate(zip(gimbled_ETPs, simmed_ETPs)):
            summed_combo=np.sum(sim_combo, axis=0)
            total_reps = np.sum(summed_combo)
            freqs_simmed_ETPs = summed_combo / total_reps
            counts_gimble_ETPs = gimble_combo * total_reps
            assert summed_combo.shape == gimble_combo.shape
            if cmd_plot:
                af.scatter_loglog(freqs_simmed_ETPs, gimble_combo, -np.log(1 / total_reps), '%s_rep%s' % (config['gimble']['model'], idx))
            af.chisquare(summed_combo, counts_gimble_ETPs)

# @pytest.mark.chunks
# class Test_n_blocks_chunks:
    
#   @pytest.fixture(scope='class')
#   def sim_ETPs_basic(self):
#       n, blocks, chunks = 1, 1, 1
#       replicates = 1000
#       yield af.sim_ETPs(**all_configs['IM_AB'], n=n, blocks=blocks, chunks=chunks, replicates=replicates)
    
#   def test_sample_size(self, sim_ETPs_basic, cmd_plot):
#       n, blocks, chunks = 4, 1, 1
#       replicates = 1000
#       test_replicates = af.sim_ETPs(
#           **all_configs['IM_AB'], 
#           n=n, 
#           blocks=blocks, 
#           chunks=chunks, 
#           replicates=replicates
#           )
#       self.run_actualtest(
#           sim_ETPs_basic, 
#           test_replicates, 
#           'sim_increased_sample_size', 
#           correction=True, 
#           plot=cmd_plot
#           )
    
#   def test_chunking(self, sim_ETPs_basic):
#       chunks = 0 
#       for combo in sim_ETPs_basic:
#           combo_sum = np.reshape(np.sum(combo, axis=0), -1)
#           chunks = np.sum(combo_sum)
#           reshaped = np.reshape(lib.simulate._combine_chunks(combo, chunks),-1)
#           assert np.all(np.equal(reshaped, combo_sum))
#   """
#   def test_sim_blocking(self, sim_ETPs_basic):
#       #needs to be replaced by a test that actually verifies the blocking rather than perform a sim-sim comparison
#       #block making needs to be carved out to achieve this.
#       n, blocks, chunks = 1, 10, 1
#       replicates = 1000
#       test_replicates = af.sim_ETPs(**all_configs['IM_AB'], n=n, blocks=blocks, chunks=chunks, replicates=replicates) 
#       self.run_actualtest(sim_ETPs_basic, test_replicates, 'sim_blocking')
#   """
#   def run_actualtest(self, sim_basic, sim_test, name, correction=False, plot=False):
#       for idx, (basic_combo, test_combo) in enumerate(zip(sim_basic, sim_test)):
#           basic_counts=np.sum(basic_combo,axis=0)
#           basic_total_reps = np.sum(basic_counts)
#           test_counts=np.sum(test_combo,axis=0)
#           test_total_reps = np.sum(test_counts)   
#           freqs_basic_ETPs = basic_counts/basic_total_reps
#           freqs_test_ETPs = test_counts/test_total_reps
#           assert basic_counts.shape == test_counts.shape
#           if plot:
#               af.scatter_loglog(freqs_test_ETPs, freqs_basic_ETPs, -np.log(1/basic_total_reps), f'{name}_rep{idx}', 'n_1', 'n_4')
#           if correction:
#               af.chisquare_contingency(np.round(test_counts/(test_total_reps/basic_total_reps)), basic_counts)
#           else:
#               af.chisquare_contingency(test_counts, basic_counts)

@pytest.mark.optimize
@pytest.mark.parametrize('opt_config_file', [
    # ~/git/gIMble/gIMble makeconfig -l test_optimize_1 -t optimize -m 1
    #af.get_test_config_file(task='optimize', model='DIV', label='test_optimize_1'),
    # ~/git/gIMble/gIMble makeconfig -l test_optimize_2 -t optimize -m 3
    af.get_test_config_file(task='optimize', model='MIG_BA', label='test_optimize_2'),
    # ~/git/gIMble/gIMble makeconfig -l test_optimize_3 -t optimize -m 4
    #af.get_test_config_file(task='optimize', model='IM_AB', label='test_optimize_3')
    ], scope='class')
@pytest.mark.parametrize('etp_config_file', [
    # ~/git/gIMble/gIMble makeconfig -l test_optimize_1 -t optimize -m 1
    #af.get_test_config_file(task='optimize', model='DIV', label='test_etp_1'),
    # ~/git/gIMble/gIMble makeconfig -l test_optimize_2 -t optimize -m 3
    af.get_test_config_file(task='optimize', model='MIG_BA', label='test_etp_2'),
    # ~/git/gIMble/gIMble makeconfig -l test_optimize_3 -t optimize -m 4
    #af.get_test_config_file(task='optimize', model='IM_AB', label='test_etp_3')
    ], scope='class')
class Test_optimize:

    @pytest.fixture(scope='class')
    def gimble_ETPs(self):
        def _gimble_ETPs(config):
            gf = lib.math.config_to_gf(config)
            gfEvaluatorObj = togimble.gfEvaluator(gf, config['max_k'], lib.gimble.MUTYPES, config['gimble']['precision'], exclude=[(2,3),])
            ETPs = lib.math.new_calculate_all_ETPs(
                gfEvaluatorObj, 
                config['parameters_expanded'], 
                config['populations']['reference_pop_id'], 
                config['block_length'], 
                config['mu']['mu'], 
                processes=1, 
                verbose=False
                )
            print(ETPs.shape)
            return ETPs
        return _gimble_ETPs              

    def test_known_truth(self, gimble_ETPs, opt_config_file, etp_config_file):
        # ETPs
        config_ETPs = lib.gimble.load_config(etp_config_file)
        config_ETPs['block_length'] = 64
        dataset = gimble_ETPs(config_ETPs)
        # optimize
        config_optimize = lib.gimble.load_config(opt_config_file)
        config_optimize['start_point'] = np.array([6e5, 7e-7])
        config_optimize['block_length'] = 64
        config_optimize['num_cores'] = 1
        config_optimize['max_iterations'] = 10
        config_optimize['xtol_rel'] = -1
        config_optimize['ftol_rel'] = 0.1
        gf = lib.math.config_to_gf(config_optimize)
        gfEvaluatorObj = togimble.gfEvaluator(gf, config_optimize['max_k'], lib.gimble.MUTYPES, config_optimize['gimble']['precision'], exclude=[(2,3),])
        # truth
        truth_values = {'Ne_B': 6e5, 'me': 7e-7, 'Ne_A': 1.3e6}
        for data_idx, data in ((0, dataset) for _ in (0,)):
            optimize_result = lib.math.optimize(gfEvaluatorObj, data_idx, data, config_optimize)
            print(optimize_result)
            truth_comparison = np.array([optimize_result['nlopt_values_by_dataset_idx'][0][k] - v for k, v in truth_values.items()])
            print(truth_comparison)
            assert np.all(np.isclose(truth_comparison, np.zeros((len(truth_values))), atol=1e-3))