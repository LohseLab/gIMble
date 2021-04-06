import itertools
import numpy as np
import pytest
import sage.all

import lib.math
import lib.gimble
import lib.simulate
import tests.aux_functions as af

from lib.GeneratingFunction.gf import togimble

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
			model = lib.gimble.get_model_name(global_info['model_file'])
			mutype_labels, max_k = zip(*sorted(global_info['k_max'].items()))
			gf = lib.math.config_to_gf(model, mutype_labels)
			gfEvaluatorObj = togimble.gfEvaluator(gf, max_k, mutype_labels)
			parameter_combinations = lib.gimble.LOD_to_DOL(sim_configs)
			ETPs = lib.math.new_calculate_all_ETPs(
				gfEvaluatorObj,
				parameter_combinations,
				global_info['reference_pop'],
				global_info['blocklength'],
				global_info['mu']
				)
			return ETPs
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
		test_replicates = af.sim_ETPs(
			**all_configs['IM_AB'], 
			n=n, 
			blocks=blocks, 
			chunks=chunks, 
			replicates=replicates
			)
		self.run_actualtest(
			sim_ETPs_basic, 
			test_replicates, 
			'sim_increased_sample_size', 
			correction=True, 
			plot=cmd_plot
			)
	
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
			if plot:
				af.scatter_loglog(freqs_test_ETPs, freqs_basic_ETPs, -np.log(1/basic_total_reps), f'{name}_rep{idx}', 'n_1', 'n_4')
			if correction:
				af.chisquare_contingency(np.round(test_counts/(test_total_reps/basic_total_reps)), basic_counts)
			else:
				af.chisquare_contingency(test_counts, basic_counts)

@pytest.mark.optimize
class Test_optimize:
	model_to_test = 'MIG_BA'

	@pytest.fixture(scope='class')
	def get_equationSystemObj(self):
		model_file = all_configs[Test_optimize.model_to_test]['global_info']['model_file']  
		k_max = all_configs[Test_optimize.model_to_test]['global_info']['k_max']
		model = lib.gimble.get_model_name(model_file)
		mutype_labels, max_k = zip(*sorted(k_max.items()))
		gf = lib.math.config_to_gf(model, mutype_labels)
		gfEvaluatorObj = togimble.gfEvaluator(gf, max_k, mutype_labels)		
		yield gfEvaluatorObj

	@pytest.fixture(scope='class')
	def gimble_ETPs(self, get_equationSystemObj):
		parameter_combinations = lib.gimble.LOD_to_DOL(all_configs[Test_optimize.model_to_test]['sim_configs'])
		global_info = all_configs[Test_optimize.model_to_test]['global_info']		
		ETPs = lib.math.new_calculate_all_ETPs(
			get_equationSystemObj,
			parameter_combinations,
			global_info['reference_pop'],
			global_info['blocklength'],
			global_info['mu']
			)
		return ETPs

	def optimize_parameters(self, parameter_combinations, numPoints, gfEvaluatorObj, ETPs, truth=None):
		"""
		method is stripped down version of lib.math.optimize_parameters: 
		ideally same function, but dependency of parameterObj and printing logfiles --> can it be avoided?
		"""
		max_eval, xtol_rel, ftol_rel = 10, -1, 0.1
		threads=1
		verbose=False
		reference_pop = global_info = all_configs[Test_optimize.model_to_test]['global_info']['reference_pop']
		#fixedParams list of parameters that are fixed including ,mu
		fixedParams = {k:v for k,v in all_configs[Test_optimize.model_to_test]['sim_configs'][0].items() if k not in parameter_combinations.keys()}
		mu = all_configs[Test_optimize.model_to_test]['global_info']['mu']
		block_length = all_configs[Test_optimize.model_to_test]['global_info']['blocklength'] 
		#boundaryNames: list of parameters that are not fixed
		boundaryNames = list(parameter_combinations.keys())
		starting_points, parameter_combinations_lowest, parameter_combinations_highest = lib.math._optimize_get_boundaries(parameter_combinations, boundaryNames, numPoints)
		if truth is not None:
			starting_points = truth
		trackHistoryPath = [[] for _ in range(numPoints)]
		#specify the objective function to be optimized
		specified_objective_function_list = lib.math._optimize_specify_objective_function(
			gfEvaluatorObj, 
			'simulate', 
			trackHistoryPath, 
			ETPs, 
			boundaryNames,
			fixedParams, 
			mu, 
			block_length, 
			reference_pop, 
			verbose
			)
		allResults=[]
		for startPos, specified_objective_function in zip(itertools.cycle(starting_points), specified_objective_function_list):
			single_run_result = lib.math.run_single_optimize(
				startPos, 
				parameter_combinations_lowest, 
				parameter_combinations_highest, 
				specified_objective_function, 
				max_eval, 
				xtol_rel, 
				ftol_rel
				) 
			allResults.append(single_run_result)
		return allResults                    

	def test_known_truth(self, get_equationSystemObj, gimble_ETPs):
		known_truth = np.array([6e5, 7e-7])
		parameter_combinations = {'Ne_B':[3e5, 9e5], 'me_B_A':[7e-8,7e-6]}
		numPoints=1
		#gimble_ETPs is data
		optimized_results = self.optimize_parameters(
			parameter_combinations, 
			numPoints, 
			get_equationSystemObj, 
			gimble_ETPs, 
			truth=[known_truth,]
			)
		print(optimized_results)
		compare_truth_optimize = known_truth-optimized_results[0]['optimum']
		assert np.all(np.isclose(compare_truth_optimize, np.zeros((len(parameter_combinations))), atol=1e-3))