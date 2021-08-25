from p_from_scaled_jaccard import compute_confidence_interval_one_step, compute_confidence_interval_two_step
from mutation_model_simulator import MutationModel
import numpy as np
import sys

if __name__ == '__main__':
    mutation_rates = [0.01]
    scale_factors = [0.1]
    confidence = 0.95
    k_mer_lengths = [21]
    num_k_mers_list = [10000000]
    num_simulations = 1000

    for scale_factor in scale_factors:
        max_trials = int(1.0/scale_factor)+1
        for num_k_mers in num_k_mers_list:
            for mutation_rate in mutation_rates:
                for k_mer_length in k_mer_lengths:    
                    mutation_model = MutationModel(num_k_mers+k_mer_length-1, k_mer_length, mutation_rate, scale_factor)
                    scaled_jaccard_indices = []
                    for iter_count in range(num_simulations):
                        counter = 0
                        while counter < max_trials:
                            mutation_model.generate()
                            candidate_scaled_jaccard_index = mutation_model.count_scaled_jaccard()
                            if candidate_scaled_jaccard_index <= 0.0001 or candidate_scaled_jaccard_index >= 0.9999:
                                counter += 1
                            else:
                                break
                        scaled_jaccard_indices.append(candidate_scaled_jaccard_index)
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
                        if mutation_rate >= plow and mutation_rate <= phigh:
                            hit_count += 1
                        else:
                            miss_count += 1
                        mid_point = (plow + phigh) / 2.0
                        if mutation_rate >= mid_point:
                            above_mid_point_count += 1
                        else:
                            below_mid_point_count += 1
                    print(k_mer_length, mutation_rate, num_k_mers, scale_factor, 1.0*hit_count/(hit_count+miss_count), hit_count, miss_count, above_mid_point_count, below_mid_point_count)