#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import allel
import pandas as pd
from timeit import default_timer as timer

def block_sites_to_variation_arrays(block_sites, max_type_count=4):
    # to make it general we must be able to calculate number of mutypes based on sample-set-length and ploidy
    block_count = block_sites.shape[0]
    max_type_count += 3 # ideally this should be the maximum amount of mutypes + 2 + 1 
    temp_sites = block_sites + (max_type_count * np.arange(block_count).reshape(block_count, 1))
    # return multiallelic, missing, monomorphic, variation
    return np.hsplit(np.bincount(temp_sites.ravel(), minlength=(block_count * max_type_count)).reshape(-1, max_type_count), [1, 2, 3])

def szudzik_pairing(folded_minor_allele_counts):
    '''1=HetB, 2=HetA, 3=HetAB, 4=Fixed'''
    # adapted from: https://drhagen.com/blog/superior-pairing-function/
    # for unpairing and multidimensional pairing functions, see: https://drhagen.com/blog/multidimensional-pairing-functions/
    if isinstance(folded_minor_allele_counts, np.ndarray):
        return np.where(
            (folded_minor_allele_counts[:,0] >= folded_minor_allele_counts[:,1]),
            np.square(folded_minor_allele_counts[:,0]) + folded_minor_allele_counts[:,0] + folded_minor_allele_counts[:,1],
            folded_minor_allele_counts[:,0] + np.square(folded_minor_allele_counts[:,1])
            )
    elif isinstance(folded_minor_allele_counts, tuple):
        a, b = folded_minor_allele_counts
        if a >= b:
            return (a**2) + a + b 
        return a + (b**2)
    else:
        pass

def genotype_to_mutype_array(sa_genotype_array, block_sites_variant_bool, block_sites, debug=False):

    np_genotype_array = np.array(sa_genotype_array)
    np_allele_count_array = np.ma.masked_equal(sa_genotype_array.count_alleles(), 0, copy=False)    
    allele_map = np.ones((np_allele_count_array.shape), dtype='int8') * np.arange(np_allele_count_array.shape[-1])
    idx_max_global_allele_count = np.nanargmax(np_allele_count_array, axis=1)
    idx_min_global_allele_count = np.nanargmin(np_allele_count_array, axis=1)
    has_major_allele = (idx_max_global_allele_count != idx_min_global_allele_count)
    idx_min_prime_allele = np.amin(np_genotype_array[:,0], axis=1)
    idx_min_global_allele = np.amin(np.amin(np_genotype_array, axis=1), axis=1)
    idx_max_global_allele = np.amax(np.amax(np_genotype_array, axis=1), axis=1)
    idx_major_allele = np.where(
        has_major_allele, 
        idx_max_global_allele_count, 
        idx_min_prime_allele)
    idx_minor_allele = np.where(
        has_major_allele, 
        idx_min_global_allele_count, 
        np.where((
            idx_min_global_allele == idx_min_prime_allele),
            np.max((idx_min_global_allele, idx_max_global_allele), axis=0), 
            np.min((idx_min_global_allele, idx_max_global_allele), axis=0)))
    # for each genotype (np.arange(allele_map.shape[0])), set minor allele to 1 (1st do minor, so that overwritten if monomorphic)
    allele_map[np.arange(allele_map.shape[0]), idx_minor_allele] = 1 
    # for each genotype (np.arange(allele_map.shape[0])), set major allele to 0
    allele_map[np.arange(allele_map.shape[0]), idx_major_allele] = 0
    folded_minor_allele_counts = sa_genotype_array.map_alleles(allele_map).to_n_alt(fill=-1)
    folded_minor_allele_counts[np.any(sa_genotype_array.is_missing(), axis=1)] = np.ones(2) * -1        # -1, -1 for missing => -1
    folded_minor_allele_counts[(np_allele_count_array.count(axis=1) > 2)] = np.ones(2) * (-1, -2)       # -1, -2 for multiallelic => -2
    block_sites_pos = block_sites.flatten()
    block_sites[block_sites_variant_bool] = szudzik_pairing(folded_minor_allele_counts) + 2               # add 2 so that not negative for bincount
    block_sites[~block_sites_variant_bool] = 2                                                            # monomorphic = 2 (0 = multiallelic, 1 = missing)
    if debug == True:
        pos_df = pd.DataFrame(block_sites_pos[block_sites_variant_bool.flatten()], dtype='int', columns=['pos'])
        genotypes_df = pd.DataFrame(np_genotype_array.reshape(np_genotype_array.shape[0], 4), dtype='i4', columns=['a1', 'a2', 'b1', 'b2'])        
        block_sites_df = pos_df.join(genotypes_df)
        folded_minor_allele_count_df = pd.DataFrame(folded_minor_allele_counts, dtype='int8', columns=['fmAC_a', 'fmAC_b'])
        block_sites_df = block_sites_df.join(folded_minor_allele_count_df)
        variants = pd.DataFrame(block_sites[block_sites_variant_bool], dtype='int', columns=['SVar'])
        block_sites_df = block_sites_df.join(variants)
        print('[#] Debugging info for variant sites:\n [#] Mutypes: 0=MULTI, 1=MISS, 2=MONO, 3=HetB, 4=HetA, 5=HetAB, 6=Fixed')
        print(block_sites_df)
    return block_sites

if __name__ == '__main__':
    '''
    gts                         := numpy array of genotypes
                                    - shape = (n=sites, s=samples, p=ploidy)
    block_sites                 := numpy array of blocked sites 
                                    - shape = (n=sites, b=block_length) 
                                    - this is based on BED file of intervals
    block_sites_variant_bool    := numpy array (boolean) of whether site in block_sites exists in VCF file 
                                    - shape = block_sites
                                    - accounts for monomorphic sites are not recorded in VCF
    '''
    run = 'basic'

    if run == 'basic':
        gts = np.array(
                        [[[0,0], [0,0]],
                        [[0,0], [0,1]],
                        [[0,1], [0,0]],
                        [[0,1], [0,1]],
                        [[0,0], [1,1]],
                        [[0,0], [2,2]],
                        [[1,1], [2,2]],
                        [[1,2], [1,2]],
                        [[1,2], [2,2]],
                        [[2,2], [1,2]],
                        [[-1,-1], [0,0]],
                        [[0,1], [2,3]]]
                        )
    
        pos = np.arange(1, 2*gts.shape[0], 2) # positions of GTs in VCF (only for didactic purpose)
        sa_genotype_array = allel.GenotypeArray(gts)
        print("[+] %s variants in VCF file on the following positions:\n%s" % (gts.shape[0], str(list(pos))))
        
        block_sites = np.arange(25).reshape(5,5)
        print("[+] block_sites inferred from BED file: \n%s" % (block_sites))
        
        block_sites_variant_bool = np.isin(block_sites, pos, assume_unique=True)
        print("[+] block_sites in VCF file: \n%s" % (block_sites_variant_bool))
    
        block_sites = genotype_to_mutype_array(sa_genotype_array, block_sites_variant_bool, block_sites, debug=True)
        print("[+] Actual data used in downstream analyses of gimble: \n%s" % (block_sites))    
    
        multiallelic, missing, monomorphic, variation = block_sites_to_variation_arrays(block_sites)
        print("[+] Variation (bSFS) array [IDX: 1=hetB, 2=hetA, 3=hetAB, 4=fixed]: \n%s" % (variation))    
    else:
        parameters = [
            (64, 5*10**5, 1/100),
            (64, 5*10**5, 1/100),
            (64, 5*10**5, 1/100),
            (64, 5*10**5, 1/100),
        ]
        for block_size, block_count, snp_density in parameters:
            # setup
            genome_sites = np.arange(block_size * block_count)
            block_sites = genome_sites.reshape(-1, block_size)
            genotypes = np.random.randint(low=0, high=2, size=(int(genome_sites.shape[0] * snp_density), 2, 2))
            pos = np.random.choice(genome_sites, genotypes.shape[0], replace=False)
            block_sites_variant_bool = np.isin(block_sites, pos, assume_unique=True)
            # work
            start_time = timer()
            sa_genotype_array = allel.GenotypeArray(genotypes)
            block_sites = genotype_to_mutype_array(sa_genotype_array, block_sites_variant_bool, block_sites, debug=False)
            multiallelic, missing, monomorphic, variation = block_sites_to_variation_arrays(block_sites)
            print("[+] block_size=%s, block_count=%s, SNP-density=%s => %.3fs" % (block_size, block_count, snp_density, timer() - start_time))
            
        