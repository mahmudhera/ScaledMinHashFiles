from mutation_model_simulator import MutationModel
from math import sqrt, log, exp, erf
import kmer_mutation_formulas_thm5 as thm5
from scipy.stats import norm as scipy_norm
from third_moment_calculator import *
import third_moment_calculator as third

def probit(p):
    return scipy_norm.ppf(p)

def r1_to_q(k,r1):
	q = 1-(1-r1)**k
	return float(q)

def exp_n_mutated(L,k,r1):
	q = r1_to_q(k,r1)
	return L*q

def var_n_mutated(L,k,r1):
	if (r1 == 0): return 0.0
	r1 = float(r1)
	q = r1_to_q(k,r1)
	varN = L*(1-q)*(q*(2*k+(2/r1)-1)-2*k) \
	     + k*(k-1)*(1-q)**2 \
	     + (2*(1-q)/(r1**2))*((1+(k-1)*(1-q))*r1-q)
	assert (varN>=0.0)
	return float(varN)

def exp_n_mutated_squared(L, k, p):
    return var_n_mutated(L, k, p) + exp_n_mutated(L, k, p) ** 2

def f1(L, N, z_alpha, s):
    ex_j = 1.0 * (L-N) / (L+N)
    bf = 1 - (1 - s) ** (L + N)
    var_j = 2.0 * N * (L - N) * (1 - s) / ( s * (L + N)**3 )
    return ex_j + z_alpha * sqrt (var_j / (bf ** 2))

def derivative_f1(L, N, z_alpha, s):
    return -(2*L)/(L+N)**2 + z_alpha * sqrt((1-s)/s) * sqrt(1.0*(L+N)/(2*N*(L-N))) * 1.0 * (L*L + N*N - 6*N*L)/(L+N)**3

def f2(L, N, z_alpha, s):
    ex_j = 1.0 * (L-N) / (L+N)
    bf = 1 - (1 - s) ** (L + N)
    var_j = 2.0 * N * (L - N) * (1 - s) / ( s * (L + N)**3 )
    return ex_j - z_alpha * sqrt (var_j / (bf ** 2))

def func(L, N):
    return L*L + N*N - 6*N*L


def variance_scaled_jaccard(L, p, k, s):
    exp_n_mut = thm5.exp_n_mutated(L, k, p)
    exp_n_mut_squared = third.exp_n_mutated_squared(L, k, p)
    exp_n_mut_cubed = third.exp_n_mutated_cubed(L, k, p)
    bias_factor = 1 - (1 - s) ** ( int(L + exp_n_mut) )

    factor1 = (1-s)/(s * bias_factor**2)
    factor2 = (2 * L * exp_n_mut - 2 * exp_n_mut_squared) / (L ** 3 + 3*L*exp_n_mut_squared + 3*L*L*exp_n_mut + exp_n_mut_cubed)
    term1 = factor1 * factor2
    term2 = (L**2 - 2 * L * exp_n_mut + exp_n_mut_squared) / (L**2 + 2 * L * exp_n_mut + exp_n_mut_squared)
    term3 = ((L - exp_n_mut) / (L + exp_n_mut))**2

    return term1 + term2 - term3

def exp_probability_path_case_taylor(L, k, p, s):
    c = exp_n_mutated(L, k, p)
    term1 = (1-s)**(L-c)
    term2 = -1.0 * (1-s)**(L-c) * log(1-s) * (exp_n_mutated(L, k, p) - c)
    term3 = 0.5 * (1-s)**(L-c) * ( log(1-s) )**2 * (exp_n_mutated_squared(L, k, p) - c**2)
    term4 = -1.0/6.0 * (1-s)**(L-c) * ( log(1-s) )**3 * (exp_n_mutated_cubed(L, k, p) - 3*c*exp_n_mutated_squared(L, k, p) + 2*c**3)
    #print (term1, term2, term3, term4)
    return sum([term1, term2, term3, term4])

def exp_probability_path_case_are_correction(L, k, p, s):
    mu = L * (1-p)**k
    sigma = sqrt(var_n_mutated(L, k, p))
    f = lambda a, b: 0.5 * ( erf( (mu-a)/(sqrt(2)*sigma) ) - erf( (mu-b)/(sqrt(2)*sigma) ) )
    delta = f(-0.5, L+0.5)
    probability = lambda n: f(n-0.5, n+0.5)/delta
    expected_probability = 0.0
    for n in range(L+1):
        expected_probability += ((1-s)**n)*probability(n)
    return expected_probability

def exp_probability_path_case_david(L, k, p, s):
    return (1-s)**(L*(1-p)**k)

def exp_probability_path_case(L, k, p, s):
    try:
        t = log (1-s)
    except:
        print (s)
    t1 = t * L * (1-p)**k
    t2 = 0.5 * var_n_mutated(L, k, p) * t * t

    #return exp(t1 + t2)

    mu = L * (1-p)**k
    sigma = sqrt(var_n_mutated(L, k, p))

    if sigma <= 0.0:
        return exp(t1 + t2)/2.0

    value1 = erf( (mu + sigma**2 * t)/(sqrt(2) * sigma) )
    value2 = erf( (mu + sigma**2 * t - L)/(sqrt(2) * sigma) )

    if value1 == value2:
        return 0
    try:
        return (value1 - value2) * exp(t1 + t2) / 2.0
    except:
        print (value1)
        print (value2)

def get_percentage_deviation(estimated, target):
    return 100.0 * (abs(estimated-target))/(abs(target))


def nothing_common_between_sketches(mutation_model):
    mutated_list = []
    for pos in range(len(mutation_model.errorSeq) - mutation_model.kmerSize + 1):
        num_errors_in_k_mer = sum(mutation_model.errorSeq[pos:pos+mutation_model.kmerSize])
        if num_errors_in_k_mer > 0:
            mutated_list.append(1)
        else:
            mutated_list.append(0)
    from random import shuffle
    shuffle(mutated_list)
    scaled_length = int(len(mutated_list) * mutation_model.scale_factor)
    n_mutated_scaled = sum(mutated_list[:scaled_length])
    return n_mutated_scaled==scaled_length


'''
For 10k simulations, calculate if nothing common in the sketches. Count # of times
this happens. Count the % of times this happens. Then compare and make a contrast
with the value we have from formula.
'''
if __name__ == '__main__':
    mutation_rates = [0.2, 0.3]
    scale_factors = [0.1]
    k_mer_lengths = [21, 31]
    num_k_mers_list = [10000]
    num_simulations = 2000

    for scale_factor in scale_factors:
        for num_k_mers in num_k_mers_list:
            for mutation_rate in mutation_rates:
                for k_mer_length in k_mer_lengths:
                    mutation_model = MutationModel(num_k_mers+k_mer_length-1, k_mer_length, mutation_rate, scale_factor)
                    nothing_common_counter = 0
                    for iter in range(num_simulations):
                        mutation_model.generate()
                        if nothing_common_between_sketches(mutation_model):
                            nothing_common_counter += 1
                    estimated_expected_probability = 1.0 * nothing_common_counter / num_simulations
                    print (num_k_mers, k_mer_length, mutation_rate, scale_factor, estimated_expected_probability)
