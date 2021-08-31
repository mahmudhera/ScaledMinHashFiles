from mpl_toolkits import mplot3d
from math import sqrt, log, exp, erf
import kmer_mutation_formulas_thm5 as thm5
from scipy.stats import norm as scipy_norm
from third_moment_calculator import *
from matplotlib import pyplot as plt
import third_moment_calculator as third
from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np

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

def exp_probability_path_case(L, k, p, s):
    try:
        t = log (1-s)
    except:
        print (s)
    t1 = t * L * (1-p)**k
    t2 = 0.5 * var_n_mutated(L, k, p) * t * t
    
    #return exp(t1 + t2)
    
    mu = L * (1.0-p)**k
    sigma = sqrt(var_n_mutated(L, k, p))
    
    #return mu + sigma**2 * t - L
    
    if sigma <= 0.0:
        return exp(t1 + t2)
    
    value1 = erf( (mu + sigma**2 * t + L)/(sqrt(2) * sigma) )
    value2 = erf( (mu + sigma**2 * t - L)/(sqrt(2) * sigma) )
    
    #print (value1, value2)
    
    if value1 == value2:
        return 0
    try:
        return (value1 - value2) * exp(t1 + t2) / 2.0
    except:
        print (value1)
        print (value2)

confidence = 0.95
alpha = 1 - confidence
z_alpha = probit(1-alpha/2)
L = 1000
s = 0.99
k = 21
p = 0.35


probabilities, scale_factors = np.mgrid[0.01:0.99:100j, 0.01:0.99:100j]
k_mer_sizes = [ 21, 31, 51 ]
k_mer_size = 21

'''
V = np.vectorize(exp_probability_path_case)

elevation_angle, z_axis_angle = 10, 150

fig = plt.figure(figsize=(20,8))
ax = fig.add_subplot(1, 3, 1, projection='3d')
ax.view_init(elevation_angle, z_axis_angle)
ax.plot_surface(probabilities, scale_factors, V(L, k_mer_size, probabilities, scale_factors), rstride=1, cstride=1,
                cmap='viridis', edgecolor='none')
ax.set_title('K = 21')
ax.set_xlabel('Mutation rate')
ax.set_ylabel('Scale factor')
ax.set_zlabel('Expected probability')

k_mer_size = 31
ax = fig.add_subplot(1, 3, 2, projection='3d')
ax.view_init(elevation_angle, z_axis_angle)
ax.plot_surface(probabilities, scale_factors, V(L, k_mer_size, probabilities, scale_factors), rstride=1, cstride=1,
                cmap='viridis', edgecolor='none')
ax.set_title('K = 31')
ax.set_xlabel('Mutation rate')
ax.set_ylabel('Scale factor')
ax.set_zlabel('Expected probability')

k_mer_size = 51
ax = fig.add_subplot(1, 3, 3, projection='3d')
ax.view_init(elevation_angle, z_axis_angle)
ax.plot_surface(probabilities, scale_factors, V(L, k_mer_size, probabilities, scale_factors), rstride=1, cstride=1,
                cmap='viridis', edgecolor='none')
ax.set_title('K = 51')
ax.set_xlabel('Mutation rate')
ax.set_ylabel('Scale factor')
ax.set_zlabel('Expected probability')

plt.savefig('combined4.pdf')
'''
# expected probability that there is nothing common after taking the sketch
ks = range(10,100)
ps = [1.0*i/100 for i in range(1,100)]
ss = [1.0*i/100 for i in range(1,100)]
expectations = [ exp_probability_path_case(L, k_mer_size, mut_rate, s) for mut_rate in ps ]

plt.plot(ps, expectations)
plt.show()

'''
# monotone property of scaled jaccard
var_direct = lambda pest: variance_scaled_jaccard(L, pest, k, s)
f = lambda pest: 2.0/(2- (1-pest)**k ) - 1 - z_alpha * sqrt(var_direct(pest))

p_values = [x/1000.0 for x in range(225, 500)]
f_values = [f(p) for p in p_values]
plt.plot(p_values, f_values)
plt.xlabel('Mutation rate, p')
plt.ylabel('$J_{scale}$')
plt.savefig('foo.pdf')
'''

'''
bias_factor = 1 - (1 - s) ** L
term_1 = (1.0-s) / (s * L**3 * bias_factor**2)
term_2 = lambda pest: L * exp_n_mutated(L, k, pest) - exp_n_mutated_squared(L, k, pest)
term_3 = lambda pest: var_n_mutated(L, k, pest) / (L**2)
    
var_direct = lambda pest: term_1 * term_2(pest) + term_3(pest)
f = lambda pest: (1-pest)**k - z_alpha * sqrt(var_direct(pest))

p_values = [x/1000.0 for x in range(0,1001)]
f_values = [f(p) for p in p_values]



plt.plot(p_values, f_values)
plt.show()
'''


