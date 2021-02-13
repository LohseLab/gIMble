import pytest
import numpy as np
import pandas as pd 
import lib.gimble
import tests.aux_functions as af
import itertools

@pytest.mark.gridsearch
class Test_gridsearch:
	def test_gridsearch_sims(self):
		num_replicates, num_gridpoints = 5, 10
		pos_max = tuple(i for i in range(num_replicates))
		fixed_param_grid, fixed_param_grid_value = 'm_e', 1
		gridded_params = ['Ne_A', 'Ne_B', 'm_e']
		grid_meta_dict = {
			'Ne_A':np.array(list(pos_max)+[100]*5), 
			'Ne_B':np.array(list(pos_max)+[100]*5), 
			'm_e':np.array([1]*(num_gridpoints//2)+[10]*(num_gridpoints//2))
			}
		unique_values_fixed_param = np.unique(grid_meta_dict[fixed_param_grid])
		fixed_param_grid_value_idx = np.where(unique_values_fixed_param==fixed_param_grid_value)[0][0]
		bsfs_replicates = np.random.randint(low=0, high=1000, size=(num_replicates, 4,4,4,4), dtype=np.int16)
		num_bsfs_entries = np.sum(bsfs_replicates, axis=(1,2,3,4))
		etps = bsfs_replicates/num_bsfs_entries[:,None,None,None,None]
		grid = np.zeros((num_gridpoints,4,4,4,4),dtype=np.float64)
		grid[pos_max,:] = etps
		etps_flat = etps[0].reshape(-1)
		for i in range(num_replicates,num_gridpoints):
			new_entry = np.random.permutation(etps_flat)
			while np.all(new_entry == etps_flat):
				new_entry = np.random.permutation(etps_flat)
			grid[i] = new_entry.reshape((4,4,4,4))
		df, df_fixed_param = lib.gimble._gridsearch_sims_single(
					bsfs_replicates, 
					grid, 
					fixed_param_grid, 
					gridded_params,
					grid_meta_dict, 
					label=None,
					name=None, 
					fixed_param_grid_value_idx=0, 
					output_df=False
					)
		assert df['m_e'].nunique()==1
		assert np.all(df['Ne_A'].to_numpy()==pos_max)
		assert np.all(df['Ne_B'].to_numpy()==pos_max)
		assert df_fixed_param.shape == (num_replicates, np.unique(grid_meta_dict['m_e']).size) 

	def test_gridsearch_blocks(self):
		num_gridpoints = 10
		pos_max = 0
		bsfs = np.random.randint(low=0, high=1000, size=(4,4,4,4), dtype=np.int16)
		num_bsfs_entries = np.sum(bsfs)
		etps = bsfs/num_bsfs_entries
		grid = np.zeros((num_gridpoints,4,4,4,4),dtype=np.float64)
		grid[pos_max] = etps
		etps_flat = etps.reshape(-1)
		for i in range(1,num_gridpoints):
			new_entry = np.random.permutation(etps_flat)
			while np.all(new_entry == etps_flat):
				new_entry = np.random.permutation(etps_flat)
			grid[i] = new_entry.reshape((4,4,4,4))
		result = lib.gimble.gridsearch_np(bsfs=bsfs, grids=grid)
		assert result.shape==(num_gridpoints,)
		assert np.argmax(result)==pos_max

	def test_gridsearch_windows(self):
		num_windows, num_gridpoints = 3, 10
		pos_max = tuple(i for i in range(num_windows))
		bsfs = np.random.randint(low=0, high=1000, size=(num_windows,4,4,4,4), dtype=np.int16)
		num_bsfs_entries = np.sum(bsfs, axis=(1,2,3,4))
		etps = bsfs/num_bsfs_entries[:,None,None,None,None]
		grid = np.zeros((num_gridpoints,4,4,4,4),dtype=np.float64)
		grid[pos_max,:] = etps
		etps_flat = etps[0].reshape(-1)
		for i in range(num_windows,num_gridpoints):
			new_entry = np.random.permutation(etps_flat)
			while np.all(new_entry == etps_flat):
				new_entry = np.random.permutation(etps_flat)
			grid[i] = new_entry.reshape((4,4,4,4))
		result = lib.gimble.gridsearch_np(bsfs=bsfs, grids=grid)
		assert result.shape==(num_windows,num_gridpoints)
		assert np.all(np.argmax(result, axis=1)==pos_max)

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