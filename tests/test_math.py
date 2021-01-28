import pytest
import numpy as np 
import lib.math
import lib.gimble
import lib.simulate
import tests.aux_functions as af

@pytest.mark.simple
class Test_simple:

	@pytest.fixture()
	def generate_numbers(self):
		yield np.random.rand(2,20)

	def test_plot(self, generate_numbers):
		af.scatter_xy(generate_numbers[0], generate_numbers[1], 0.0)
		assert 1==1

@pytest.mark.ETPs
class Test_IM:
	global_info = {'mu':1e-8, 
						'ploidy': 2, 
						'sample_pop_ids': ['A','B'], 
						'blocklength': 100,
						'k_max': {'m1':2, 'm2':2, 'm3':2, 'm4':2},
						'reference_pop': 'A_B'
						}
	sim_configs = [{'Ne_A': 10000 , 'Ne_B': 20000, 'Ne_A_B': 10000, 'T': 1e6, 'me_A_B':1e-7, 'recombination':0}]

	@pytest.fixture(scope='class')
	def gimble_ETPs(self):
		EquationSystemObj = lib.math.EquationSystemObj(
			model="models/IM_AB.tsv", 
			reference_pop="A_B", 
			k_max_by_mutype={'m1':2, 'm2':2, 'm3':2, 'm4':2}, 
			block_length=100, 
			mu=1e-7
			)
		equationSystem.initiate_model(
			sync_ref=None, 
			sync_targets=None
			)
		#parameter_combinations: DOL
		equationSystem.ETPs = equationSystem.calculate_all_ETPs(
			parameter_combinations, 
			threads=1, 
			gridThreads=1
			)

		yield equationSystem.ETPs

	@pytest.fixture(scope='function') #need to pass multiple parameters to this one
	def sim_ETPs(self):

		def _sim_ETPs(n, blocks, chunks, threads=1):
			#add sample_pop_sizes, blocks, replicates to global_info dict
			add_global_info = {
								'blocks':blocks, 
								'chunks':chunks, 
								'sample_pop_sizes':[n for _ in len(Test_IM.global_info['sample_pop_ids'])],
								'replicates':2500
								}
			all_interpop_comparisons = lib.simulate.all_interpopulation_comparisons(add_global_info['sample_pop_sizes'])
			sim_global_info = {**Test_IM.global_info, **add_global_info}
			simmed_ETPs = lib.simulate.run_sims(Test_IM.sim_configs, sim_global_info, all_interpop_comparisons, chunks, threads)
			sum_simmed_ETPs = np.sum(simmed_ETPs, axis=0)
			return sum_simmed_ETPs/np.sum(sum_simmed_ETPs)

		return _sim_ETPs

	@pytest.mark.parametrize("n, blocks, chunks",
								[(1,1,1),
								(1,10,10),
								(4,1,1)
								]
							)
	def test_ETPs_IM(self, gimble_ETPs, sim_ETPs, n, blocks, chunks):
		simmed_ETPs = Test_IM.sim_configs(n, blocks, chunks)
		assert simmed_ETPs.shape == gimble_ETPs.shape
		af.chisquare(simmed_ETPs, gimble_ETPs)

@pytest.mark.ETPs
class Test_MIG:
	@pytest.fixture(scope='class')
	def gimble_ETPs(self):
		pass

	@pytest.fixture(scope='class') #need to pass multiple parameters to this one
	def sim_ETPs(self):
		pass

	def test_ETPs_MIG(self, gimble_ETPs, sim_ETPs):
		pass

@pytest.mark.ETPs
class Test_DIV:
	@pytest.fixture(scope='class')
	def gimble_ETPs(self):
		pass

	@pytest.fixture(scope='class')
	def sim_ETPs(self):
		pass

	def test_EPTs_DIV(self, gimble_ETPs, sim_ETPs):
		pass