import pytest
import numpy as np
import pandas as pd 
import lib.gimble
import tests.aux_functions as af
import itertools

@pytest.mark.gridsearch
class Test_gridsearch:
	"""
	test1: 
	* generate ETPs: simulate or exp
	* populate other matrices with similar values (shuffle)
	* build grid from those arrays
	* generate data from those ETPs -> multiple windows
	* return best gridpoint

	gridsearch can involve: blocks, windows, simulate
	* get_bsfs
	* gridsearch_np
	* find local/global winning parameter

	sims:
	*what happens to fixed parameter?
	*gridsearch sims_single: for single parameter
		* multiplies replicates with bSFS -> result replicates*num_grid_points
		* optimal index for each replicate
		* if fixed_param_grid is specified: returns df with distribution of lncls
			for each value of the fixed parameter

	test2: test simulate fixed grid

	"""
	def test_gridsearch(self):
		pass

@pytest.mark.simgrid
class Test_simgrid:
	"""
	test simulate --grid

	"""
	fixed_param_grid_value=1
	@pytest.fixture(scope='class')
	def params(self):
		grid_meta_dict = {'Ne_A':np.linspace(10,100,10), 'Ne_B':np.linspace(10,100,10)[::-1],'m_e':np.array([Test_simgrid.fixed_param_grid_value]*5+[5]*5, dtype=np.int16)}
		lncls_global = np.array([-1]+9*[-100])
		lncls_windows = np.full((5,10), -100) #lncls_windows dimensions: windows*grid_size
		local_winning_grid_points = (0,1,2,1,0)
		lncls_windows[(0,1,2,3,4),local_winning_grid_points] = -1
		lncls_windows[0,6] = -0.5
		fixed_param_grid = 'm_e'
		yield (grid_meta_dict, lncls_global, lncls_windows, fixed_param_grid, local_winning_grid_points)

	def test_simgrid_fixed_rec(self, params):
		grid_meta_dict, lncls_global, lncls_windows, fixed_param_grid, local_winning_grid_points = params
		rec_rate = 0.0
		num_windows = lncls_windows.shape[0]
		entries = {
			'sequence':[0]*num_windows, 
			'start':[i*100 for i in range(num_windows)], 
			'end':[i*100+100 for i in range(num_windows)]
				}
		window_coordinates = pd.DataFrame(entries)
		sim_grid, window_df = lib.gimble._get_sim_grid(
			grid_meta_dict, 
			lncls_global, 
			lncls_windows, 
			fixed_param_grid, 
			rec_rate,
			window_coordinates
			)
		print('sim_grid:',sim_grid)
		print('window_df:',window_df)
		fixed_param_grid_values = [combo[fixed_param_grid] for combo in sim_grid]
		assert len(set(fixed_param_grid_values))==1
		assert fixed_param_grid_values[0] == Test_simgrid.fixed_param_grid_value
		assert len(sim_grid) == window_df['param_idx'].nunique()
		assert max(window_df['param_idx'])==len(sim_grid)-1
		#for numbers with same index in tuple, param_idx should be equal as well
		window_np = window_df['param_idx'].to_numpy()
		local_winning_grid_points_np = np.array(local_winning_grid_points, dtype=np.int8)
		for idx in range(window_df['param_idx'].nunique()):
			corresponding_idxs = local_winning_grid_points_np==idx
			assert np.unique(window_np[corresponding_idxs]).size==1
		
	def test_simgrid_recmap(self, params):
		grid_meta_dict, lncls_global, lncls_windows, fixed_param_grid, local_winning_grid_points = params
		num_windows = lncls_windows.shape[0]
		entries = {
			'sequence':[0]*num_windows,
			'start':[i*100 for i in range(num_windows)],
			'end':[i*100+100 for i in range(num_windows)],
			'rec_bins':[0.001, 0.01, 0.1, 0.01, 0.1]
		}
		rec_map = pd.DataFrame(entries)
		sim_grid, window_df = lib.gimble._get_sim_grid(
			grid_meta_dict, 
			lncls_global, 
			lncls_windows, 
			fixed_param_grid, 
			rec_map
			)
		print('sim_grid:',sim_grid)
		print('window_df:',window_df)
		fixed_param_grid_values = [combo[fixed_param_grid] for combo in sim_grid] 
		assert len(set(fixed_param_grid_values))==1
		assert fixed_param_grid_values[0] == Test_simgrid.fixed_param_grid_value
		assert len(sim_grid) == window_df['param_with_rec_idx'].nunique()
		assert max(window_df['param_with_rec_idx']) == len(sim_grid) -1
		window_np = window_df['param_with_rec_idx'].to_numpy()
		local_winning_grid_points_np = np.array(local_winning_grid_points, dtype=np.int8)
		rec_bins_np = rec_map['rec_bins'].to_numpy()
		for idx in range(window_df['param_with_rec_idx'].nunique()):
			corresponding_idxs = local_winning_grid_points_np==idx
			rec_bins_selection = rec_bins_np[corresponding_idxs]
			assert np.unique(window_np[corresponding_idxs]).size==np.unique(rec_bins_selection).size
		
	def test_get_slice_grid_meta_single_value(self, params):
		#case of single gridpoint corresponding to the global optimum value for the fixed parameter
		grid_meta_dict, lncls_global, lncls_windows, fixed_param_grid, local_winning_grid_points = params
		modified_grid_meta_dict = {k:v for k,v in grid_meta_dict.items()}
		modified_grid_meta_dict['m_e'] = np.array([Test_simgrid.fixed_param_grid_value]*1+[5]*9, dtype=np.int16)
		rec_rate = 0.0
		num_windows = lncls_windows.shape[0]
		entries = {
			'sequence':[0]*num_windows, 
			'start':[i*100 for i in range(num_windows)], 
			'end':[i*100+100 for i in range(num_windows)]
				}
		window_coordinates = pd.DataFrame(entries)
		sim_grid, window_df = lib.gimble._get_sim_grid(
			modified_grid_meta_dict, 
			lncls_global, 
			lncls_windows, 
			fixed_param_grid, 
			rec_rate,
			window_coordinates
			)
		print('sim_grid:',sim_grid)
		print('window_df:',window_df)
		
		assert len(sim_grid) == window_df['param_idx'].nunique()
		assert len(sim_grid) == 1

	def test_simgrid_no_fixed_param(self, params):
		grid_meta_dict, lncls_global, lncls_windows, fixed_param_grid, local_winning_grid_points = params
		fixed_param_grid = None
		rec_rate = 0.0
		num_windows = lncls_windows.shape[0]
		entries = {
			'sequence':[0]*num_windows, 
			'start':[i*100 for i in range(num_windows)], 
			'end':[i*100+100 for i in range(num_windows)]
				}
		window_coordinates = pd.DataFrame(entries)
		sim_grid, window_df = lib.gimble._get_sim_grid(
			grid_meta_dict, 
			lncls_global, 
			lncls_windows, 
			fixed_param_grid, 
			rec_rate,
			window_coordinates
			)
		print('sim_grid:',sim_grid)
		print('window_df:',window_df)
		#add additional assertions
		fixed_param_grid_values = [combo['m_e'] for combo in sim_grid] 
		assert len(set(fixed_param_grid_values))>1