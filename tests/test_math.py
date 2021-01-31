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

@pytest.mark.ETPs_new
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

	@pytest.fixture(scope='function')
	def sim_ETPs(self):
		def _sim_ETPs(global_info, sim_configs):
			threads = 1
			blocks, chunks, n = 1, 1, 1
			add_global_info = {
								'blocks':blocks, 
								'chunks':chunks, 
								'sample_pop_sizes':[n for _ in global_info['sample_pop_ids']],
								'replicates':5000
								}
			all_interpop_comparisons = lib.simulate.all_interpopulation_comparisons(*add_global_info['sample_pop_sizes'])
			sim_global_info = {**global_info, **add_global_info}
			simmed_ETPs = lib.simulate.run_sims(sim_configs, sim_global_info, all_interpop_comparisons, chunks, threads)
			return simmed_ETPs
	
		return _sim_ETPs

	def test_ETPs_IM(self, gimble_ETPs, sim_ETPs, config):
		simmed_ETPs = sim_ETPs(**all_configs[config])
		gimbled_ETPs = gimble_ETPs(**all_configs[config])
		for idx, (gimble_combo, sim_combo) in enumerate(zip(gimbled_ETPs, simmed_ETPs)):
			summed_combo=np.sum(sim_combo,axis=0)
			total_reps = np.sum(summed_combo)
			freqs_simmed_ETPs = summed_combo/total_reps
			counts_gimble_ETPs = gimble_combo*total_reps
			assert summed_combo.shape == gimble_combo.shape
			af.scatter_loglog(freqs_simmed_ETPs, gimble_combo, -np.log(1/total_reps), f'{config}_rep{idx}')
			af.chisquare(summed_combo, counts_gimble_ETPs)

@pytest.mark.chunks
class Test_n_blocks_chunks:
	
	@pytest.fixture(scope='class')
	def sim_ETPs_basic(self, sim_ETPs):
		n, blocks, chunks = 1, 1, 1
		replicates = 5000
		yield sim_ETPs(**all_configs['IM_AB'], n=n, blocks=blocks, chunks=chunks, replicates=replicates)

	@pytest.fixture(scope='function')
	def sim_ETPs_n4(self, sim_ETPs):
		n, blocks, chunks = 5, 1, 1
		replicates = 1000
		yield sim_ETPs(**all_configs['IM_AB'], n=n, blocks=blocks, chunks=chunks, replicates=replicates)

	@pytest.fixture(scope='function')
	def sim_ETPs_chunked(self, sim_ETPs):
		n, blocks, chunks = 1, 10, 10
		replicates = 500
		yield sim_ETPs(**all_configs['IM_AB'], n=n, blocks=blocks, chunks=chunks, replicates=replicates)

	@pytest.fixture(scope='function')
	def sim_ETPs_blocking(self, sim_ETPs):
		n, blocks, chunks = 1, 50, 1
		replicates = 100
		yield sim_ETPs(**all_configs['IM_AB'], n=n, blocks=blocks, chunks=chunks, replicates=replicates) 
	
	@pytest.fixture(scope='function')
	def sim_ETPs(self):
		def _sim_ETPs(global_info, sim_configs, n, blocks, chunks, replicates):
			threads = 1
			add_global_info = {
								'blocks':blocks, 
								'chunks':chunks, 
								'sample_pop_sizes':[n for _ in global_info['sample_pop_ids']],
								'replicates':replicates
								}
			all_interpop_comparisons = lib.simulate.all_interpopulation_comparisons(*add_global_info['sample_pop_sizes'])
			sim_global_info = {**global_info, **add_global_info}
			simmed_ETPs = lib.simulate.run_sims(sim_configs, sim_global_info, all_interpop_comparisons, chunks, threads)
			return simmed_ETPs
	
		return _sim_ETPs

	def run_actualtest(self, sim_basic, sim_test, name):
		for idx, (sim_combo, test_combo) in enumerate(zip(sim_basic, sim_test)):
			summed_combo=np.sum(sim_combo,axis=0)
			total_reps = np.sum(summed_combo)
			freqs_simmed_ETPs = summed_combo/total_reps
			assert summed_combo.shape == gimble_combo.shape
			af.scatter_loglog(freqs_simmed_ETPs, gimble_combo, -np.log(1/total_reps), f'{name}_rep{idx}')
			af.chisquare(freqs_simmed_ETPs, gimble_combo)

	def test_chunking(self, sim_ETPs_basic, sim_ETPs_chunked):
		run_actualtest(self, sim_ETPs_basic, sim_ETPs_chunked, 'chunking')

	def test_sample_size(self, sim_ETPs_basic, sim_ETPs_n4):
		run_actualtest(self, sim_ETPs_basic, sim_ETPs_n4, 'increased_sample_size')

	def test_blocking(self, sim_ETP_basic, sim_ETPs_blocking):
		run_actualtest(self, sim_ETPs_basic, sim_ETPs_blocking, 'blocking')

"""
@pytest.mark.simple
class Test_simple:
	global_info = {'mu':3e-9, 
						'ploidy': 2, 
						'sample_pop_ids': ['A','B'], 
						'blocklength': 64,
						'k_max': {'m_1':2, 'm_2':2, 'm_3':2, 'm_4':2},
						'reference_pop': 'A_B'
						}
	#note: k_max dict keys should always be m_x (mind underscore)
	sim_configs = [{'Ne_A': 1.3e6 , 'Ne_B': 6e5, 'Ne_A_B': 1.5e6, 'T': 1e7, 'me_A_B':7e-7, 'recombination':0}]

	@pytest.fixture(scope='function') #need to pass multiple parameters to this one
	def sim_ETPs(self):
	
		def _sim_ETPs(n, blocks, chunks, threads=1):
			#add sample_pop_sizes, blocks, replicates to global_info dict
			add_global_info = {
								'blocks':blocks, 
								'chunks':chunks, 
								'sample_pop_sizes':[n for _ in Test_simple.global_info['sample_pop_ids']],
								'replicates':5000
								}
			all_interpop_comparisons = lib.simulate.all_interpopulation_comparisons(*add_global_info['sample_pop_sizes'])
			sim_global_info = {**Test_simple.global_info, **add_global_info}
			simmed_ETPs = lib.simulate.run_sims(Test_simple.sim_configs, sim_global_info, all_interpop_comparisons, chunks, threads)
			return simmed_ETPs
	
		return _sim_ETPs
	
	@pytest.mark.parametrize("n, blocks, chunks",
								[(1,1,1),
								#(1,10,10),
								#(4,1,1)
								]
							)
	def test_ETPs_IM(self, sim_ETPs, n, blocks, chunks):
		simmed_ETPs = sim_ETPs(n, blocks, chunks)
		print(simmed_ETPs.shape)
		for combo in simmed_ETPs:
			summed_combo=np.sum(combo, axis=0)
			total_reps = np.sum(summed_combo)
			freqs_simmed_ETPs = summed_combo/total_reps
			print(freqs_simmed_ETPs)
			assert False

@pytest.mark.ETPs
class Test_IM:
	global_info = {'mu':3e-9, 
						'ploidy': 2, 
						'sample_pop_ids': ['A','B'], 
						'blocklength': 64,
						'k_max': {'m_1':2, 'm_2':2, 'm_3':2, 'm_4':2},
						'reference_pop': 'A_B'
						}
	#note: k_max dict keys should always be m_x (mind underscore)
	sim_configs = [{'Ne_A': 1.3e6 , 'Ne_B': 6e5, 'Ne_A_B': 1.5e6, 'T': 1e7, 'me_A_B':7e-7, 'recombination':0}]

	@pytest.fixture(scope='class')
	def gimble_ETPs(self):
		equationSob = lib.math.EquationSystemObj(
			model_file="models/IM_AB.tsv", 
			reference_pop=Test_IM.global_info['reference_pop'], 
			k_max_by_mutype=Test_IM.global_info['k_max'], 
			block_length=Test_IM.global_info['blocklength'], 
			mu=Test_IM.global_info['mu']
			)
		equationSob.initiate_model(
			sync_ref=None, 
			sync_targets=None
			)
		#parameter_combinations: DOL
		parameter_combinations = lib.gimble.LOD_to_DOL(Test_IM.sim_configs)
		equationSob.ETPs = equationSob.calculate_all_ETPs(
			parameter_combinations, 
			threads=2, 
			gridThreads=2
			)

		yield equationSob.ETPs

	@pytest.fixture(scope='function') #need to pass multiple parameters to this one
	def sim_ETPs(self):

		def _sim_ETPs(n, blocks, chunks, threads=1):
			#add sample_pop_sizes, blocks, replicates to global_info dict
			add_global_info = {
								'blocks':blocks, 
								'chunks':chunks, 
								'sample_pop_sizes':[n for _ in Test_IM.global_info['sample_pop_ids']],
								'replicates':5000
								}
			all_interpop_comparisons = lib.simulate.all_interpopulation_comparisons(*add_global_info['sample_pop_sizes'])
			sim_global_info = {**Test_IM.global_info, **add_global_info}
			simmed_ETPs = lib.simulate.run_sims(Test_IM.sim_configs, sim_global_info, all_interpop_comparisons, chunks, threads)
			return simmed_ETPs

		return _sim_ETPs

	@pytest.mark.parametrize("n, blocks, chunks",
								[(1,1,1),
								#(1,10,10),
								#(4,1,1)
								]
							)
	def test_ETPs_IM(self, gimble_ETPs, sim_ETPs, n, blocks, chunks):
		simmed_ETPs = sim_ETPs(n, blocks, chunks)
		for gimble_combo, sim_combo in zip(gimble_ETPs, simmed_ETPs):
			summed_combo=np.sum(sim_combo,axis=0)
			total_reps = np.sum(summed_combo)
			freqs_simmed_ETPs = summed_combo/total_reps
			assert summed_combo.shape == gimble_combo.shape
			af.scatter_loglog(freqs_simmed_ETPs, gimble_combo, -np.log(1/total_reps), f'IM_{n}_{blocks}_{chunks}')
			af.chisquare(freqs_simmed_ETPs, gimble_combo)

@pytest.mark.ETPs
class Test_MIG:
	global_info = {'mu':3e-9, 
						'ploidy': 2, 
						'sample_pop_ids': ['A','B'], 
						'blocklength': 64,
						'k_max': {'m_1':2, 'm_2':2, 'm_3':2, 'm_4':2},
						'reference_pop': 'A_B'
						}
	#note: k_max dict keys should always be m_x (mind underscore)
	sim_configs = [{'Ne_A': 1.3e6 , 'Ne_B': 6e5, 'me_A_B':7e-7, 'recombination':0}]
	
	@pytest.fixture(scope='class')
	def gimble_ETPs(self):
		pass

	@pytest.fixture(scope='class')
	def sim_ETPs(self):
		pass

	def test_ETPs_MIG(self, gimble_ETPs, sim_ETPs):
		pass

@pytest.mark.ETPs
class Test_DIV:
	global_info = {'mu':3e-9, 
						'ploidy': 2, 
						'sample_pop_ids': ['A','B'], 
						'blocklength': 64,
						'k_max': {'m_1':2, 'm_2':2, 'm_3':2, 'm_4':2},
						'reference_pop': 'A_B'
						}
	#note: k_max dict keys should always be m_x (mind underscore)
	sim_configs = [{'Ne_A': 1.3e6 , 'Ne_B': 6e5, 'Ne_A_B': 1.5e6, 'T': 1e7, 'recombination':0}]

	@pytest.fixture(scope='class')
	def gimble_ETPs(self):
		equationSob = lib.math.EquationSystemObj(
			model_file="models/DIV.tsv", 
			reference_pop=Test_.global_info['reference_pop'], 
			k_max_by_mutype=Test_DIV.global_info['k_max'], 
			block_length=Test_DIV.global_info['blocklength'], 
			mu=Test_DIV.global_info['mu']
			)
		equationSob.initiate_model(
			sync_ref=None, 
			sync_targets=None
			)
		#parameter_combinations: DOL
		parameter_combinations = lib.gimble.LOD_to_DOL(Test_DIV.sim_configs)
		equationSob.ETPs = equationSob.calculate_all_ETPs(
			parameter_combinations, 
			threads=2, 
			gridThreads=2
			)

		yield equationSob.ETPs

	@pytest.fixture(scope='class')
	def sim_ETPs(self):
			threads = 1
			add_global_info = {
								'blocks':1, 
								'chunks':1, 
								'sample_pop_sizes':[1 for _ in Test_DIV.global_info['sample_pop_ids']],
								'replicates':5000
								}
			all_interpop_comparisons = lib.simulate.all_interpopulation_comparisons(*add_global_info['sample_pop_sizes'])
			sim_global_info = {**Test_DIV.global_info, **add_global_info}
			simmed_ETPs = lib.simulate.run_sims(Test_DIV.sim_configs, sim_global_info, all_interpop_comparisons, 1, threads)
			return simmed_ETPs

	def test_EPTs_DIV(self, gimble_ETPs, sim_ETPs):
		for gimble_combo, sim_combo in zip(gimble_ETPs, sim_ETPs):
			summed_combo=np.sum(sim_combo, axis=0)
			total_reps = np.sum(summed_combo)
			freqs_simmed_ETPs = summed_combo/total_reps
			assert summed_combo.shape == gimble_combo.shape
			af.scatter_loglog(freqs_simmed_ETPs, gimble_combo, -np.log(1/total_reps), f'IM_{n}_{blocks}_{chunks}')
			af.chisquare(freqs_simmed_ETPs, gimble_combo)
"""