import pytest
import lib.gimble
import lib.simulate
import numpy as np
import tests.aux_functions as af
import msprime
import scipy.stats

@pytest.mark.simulate
class Test_sim:

    def test_IM(self):
        Ne = 10000
        t = 2*Ne
        me = 1e-5
        M = 4*Ne*me
        T = t/(2*Ne)
        config = {
            'populations' : {
                'pop_ids': ['A', 'B', 'A_B']
                },
            'parameters_expanded' : {
                'T' : [t,],
                'Ne_A' : [Ne,],
                'Ne_B' : [Ne,],
                'Ne_A_B' : [Ne,],
                'me' : [me,]
                },
            'events' : {
                'exodus' : [(0,1,2),],
                'migration' : [(1,0),],
            }
        }
        demography = next(lib.gimble.config_to_demes_graph(config, idxs=[0,]))
        demography = msprime.Demography.from_demes(demography)
        replicates = 500
        seeds = np.random.randint(1, 2 ** 32, (replicates, 2))
        recombination_rate = 0.0
        samples = {'A':1, 'B':1}
        ploidy = 2
        sequence_length = 1000
        blocks = 1
        comparisons = lib.simulate.all_interpopulation_comparisons(*list(samples.values()))
        max_k = np.array([2,2,2,2])
        mutation_rate = 1e-6
        mutation_model = 'infinite_alleles'
        discrete = True
        result = lib.simulate.run_sim_serial(seeds, recombination_rate, samples, demography, ploidy, sequence_length, comparisons, max_k, mutation_rate, mutation_model, discrete, blocks, True, 0)
        result = np.sum(result, axis=0)
        #extract info on topologies
        #extract those containing ab = m_3 mutation
        obs_freq_ab = np.sum(result[:,:,1:,:])
        obs_freq_aabb = np.sum(result[:,:,:,1:])
        total_obs = obs_freq_ab + obs_freq_aabb
        obs_prob_ab = obs_freq_ab/total_obs
        #prob of aa,bb topology
        exp_prob_ab = (4*np.exp(-(2+M)*T)+2*M)/(3*(2+M))
        print(obs_prob_ab)
        print(exp_prob_ab)
        af.chisquare([obs_freq_ab, total_obs-obs_freq_ab], [exp_prob_ab*total_obs,(1-exp_prob_ab)*total_obs])