from p_from_scaled_jaccard import compute_confidence_interval_one_step, compute_confidence_interval_two_step
from mutation_model_simulator import MutationModel
import numpy as np
import sys

if __name__ == '__main__':
    mutation_rate = 0.01
    scale_factor = 0.1
    confidence = 0.95
    k_mer_length = 21
    num_k_mers = 1000000000
    
    jaccard_indices_filename = 'all_jaccard_indices.txt'
    f = open(jaccard_indices_filename, 'r')
    lines = f.readlines()
    
    scaled_jaccard_indices = [ float( line.rstrip() ) for line in lines ]
    print (scaled_jaccard_indices)
    
    ci_one_step_list = compute_confidence_interval_one_step(scaled_jaccard_indices, num_k_mers,
                                                                            k_mer_length, confidence, scale_factor)
                    
    hit_count = 0
    miss_count = 0
    above_mid_point_count = 0
    below_mid_point_count = 0
    for ci in ci_one_step_list:
        C_scale = ci[3]
        plow = ci[4]
        phigh = ci[5]
        p_point = ci[6]
        #print(C_scale, mutation_rate, plow, phigh)
        miss = False
        if mutation_rate >= plow and mutation_rate <= phigh:
            hit_count += 1
            miss = False
        else:
            miss_count += 1
            miss = True
        mid_point = (plow + phigh) / 2.0
        if mutation_rate >= mid_point and miss:
            above_mid_point_count += 1
        elif mutation_rate < mid_point and miss:
            below_mid_point_count += 1
    print(k_mer_length, mutation_rate, num_k_mers, scale_factor, 1.0*hit_count/(hit_count+miss_count), hit_count, miss_count, above_mid_point_count, below_mid_point_count)
