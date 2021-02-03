import pytest
import numpy as np 
import lib.math
import lib.gimble
import lib.simulate
import tests.aux_functions as af

#ini files translated to the necessary parameters to test ETPs
all_configs = {
	'IM_AB': {
		"global_info" : {
					'model_file':"models/IM_AB.tsv",
					'mu':3e-9, 
					'ploidy': 2, 
					'sample_pop_ids': ['A','B'], 
					'blocklength': 64,
					'k_max': {'m_1':2, 'm_2':2, 'm_3':2, 'm_4':2},
					'reference_pop': 'A_B'
					},
		"sim_configs": [{'Ne_A': 1.3e6 , 'Ne_B': 6e5, 'Ne_A_B': 1.5e6, 'T': 1e7, 'me_A_B':7e-7, 'recombination':0}]
		},
	'DIV' : {
		"global_info" : {
					'model_file':"models/DIV.tsv",
					'mu':3e-9, 
					'ploidy': 2, 
					'sample_pop_ids': ['A','B'], 
					'blocklength': 64,
					'k_max': {'m_1':2, 'm_2':2, 'm_3':2, 'm_4':2},
					'reference_pop': 'A_B'
					},
		"sim_configs": [{'Ne_A': 1.3e6 , 'Ne_B': 6e5, 'Ne_A_B': 1.5e6, 'T': 1e7, 'recombination':0}]
		},
	'MIG_BA' : {
		"global_info" : {
					'model_file':"models/MIG_BA.tsv",
					'mu':3e-9, 
					'ploidy': 2, 
					'sample_pop_ids': ['A','B'], 
					'blocklength': 64,
					'k_max': {'m_1':2, 'm_2':2, 'm_3':2, 'm_4':2},
					'reference_pop': 'A'
					},
		"sim_configs": [{'Ne_A': 1.3e6 , 'Ne_B': 6e5, 'me_B_A':7e-7, 'recombination':0}]
		}
	}

@pytest.fixture(scope='class')
def config(request):
	yield request.param

@pytest.mark.ETPs
@pytest.mark.parametrize('config', ['IM_AB', 'DIV', 'MIG_BA'], indirect=True, scope='class')
class Test_ETPs:

	@pytest.fixture(scope='class')
	def gimble_ETPs(self):
		def _gimble_ETPs(global_info, sim_configs):
			equationSob = lib.math.EquationSystemObj(
				model_file=global_info['model_file'], 
				reference_pop=global_info['reference_pop'], 
				k_max_by_mutype=global_info['k_max'], 
				block_length=global_info['blocklength'], 
				mu=global_info['mu']
				)
			equationSob.initiate_model(
				sync_ref=None, 
				sync_targets=None
				)
			#parameter_combinations: DOL
			parameter_combinations = lib.gimble.LOD_to_DOL(sim_configs)
			equationSob.ETPs = equationSob.calculate_all_ETPs(
				parameter_combinations, 
				threads=2, 
				gridThreads=2
				)
	
			return equationSob.ETPs
		return _gimble_ETPs

	def test_ETPs_model(self, gimble_ETPs, config, cmd_plot):
		#simmed_ETPs = sim_ETPs(**all_configs[config])
		simmed_ETPs = af.sim_ETPs(**all_configs[config], n=1, blocks=1, chunks=1, replicates=1000)
		gimbled_ETPs = gimble_ETPs(**all_configs[config])
		for idx, (gimble_combo, sim_combo) in enumerate(zip(gimbled_ETPs, simmed_ETPs)):
			summed_combo=np.sum(sim_combo,axis=0)
			total_reps = np.sum(summed_combo)
			freqs_simmed_ETPs = summed_combo/total_reps
			counts_gimble_ETPs = gimble_combo*total_reps
			assert summed_combo.shape == gimble_combo.shape
			if cmd_plot:
				af.scatter_loglog(freqs_simmed_ETPs, gimble_combo, -np.log(1/total_reps), f'{config}_rep{idx}')
			af.chisquare(summed_combo, counts_gimble_ETPs)

@pytest.mark.chunks
class Test_n_blocks_chunks:
	
	@pytest.fixture(scope='class')
	def sim_ETPs_basic(self):
		n, blocks, chunks = 1, 1, 1
		replicates = 1000
		yield af.sim_ETPs(**all_configs['IM_AB'], n=n, blocks=blocks, chunks=chunks, replicates=replicates)
	
	def test_sample_size(self, sim_ETPs_basic, cmd_plot):
		n, blocks, chunks = 4, 1, 1
		replicates = 1000
		test_replicates = af.sim_ETPs(**all_configs['IM_AB'], n=n, blocks=blocks, chunks=chunks, replicates=replicates)
		self.run_actualtest(sim_ETPs_basic, test_replicates, 'sim_increased_sample_size', correction=True, plot=cmd_plot)
	
	def test_chunking(self, sim_ETPs_basic):
		chunks = 0 
		for combo in sim_ETPs_basic:
			combo_sum = np.reshape(np.sum(combo, axis=0), -1)
			chunks = np.sum(combo_sum)
			reshaped = np.reshape(lib.simulate._combine_chunks(combo, chunks),-1)
			assert np.all(np.equal(reshaped, combo_sum))
	"""
	def test_sim_blocking(self, sim_ETPs_basic):
		#needs to be replaced by a test that actually verifies the blocking rather than perform a sim-sim comparison
		#block making needs to be carved out to achieve this.
		n, blocks, chunks = 1, 10, 1
		replicates = 1000
		test_replicates = af.sim_ETPs(**all_configs['IM_AB'], n=n, blocks=blocks, chunks=chunks, replicates=replicates) 
		self.run_actualtest(sim_ETPs_basic, test_replicates, 'sim_blocking')
	"""
	def run_actualtest(self, sim_basic, sim_test, name, correction=False, plot=False):
		for idx, (basic_combo, test_combo) in enumerate(zip(sim_basic, sim_test)):
			basic_counts=np.sum(basic_combo,axis=0)
			basic_total_reps = np.sum(basic_counts)
			test_counts=np.sum(test_combo,axis=0)
			test_total_reps = np.sum(test_counts)	
			freqs_basic_ETPs = basic_counts/basic_total_reps
			freqs_test_ETPs = test_counts/test_total_reps
			assert basic_counts.shape == test_counts.shape
			af.scatter_loglog(freqs_test_ETPs, freqs_basic_ETPs, -np.log(1/basic_total_reps), f'{name}_rep{idx}', 'n_1', 'n_4')
			if correction:
				af.chisquare_contingency(np.round(test_counts/(test_total_reps/basic_total_reps)), basic_counts)
			else:
				af.chisquare_contingency(test_counts, basic_counts)