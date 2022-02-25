from scipy.special import hyp2f1
from math import log, sqrt, exp, factorial
from matplotlib import pyplot as plt
import numpy as np
from scipy.stats import norm

def r1_to_q(k,r1):
	#return 1-(1-r1)**k
	r1 = float(r1)
	q = 1-(1-r1)**k
	return float(q)

def exp_n_mutated(L,k,r1):
	q = r1_to_q(k,r1)
	return L*q

def var_n_mutated(L,k,r1,q=None):
	if (r1 == 0): return 0.0
	r1 = float(r1)
	if (q == None):
		q = r1_to_q(k,r1)
	varN = L*(1-q)*(q*(2*k+(2/r1)-1)-2*k) \
	     + k*(k-1)*(1-q)**2 \
	     + (2*(1-q)/(r1**2))*((1+(k-1)*(1-q))*r1-q)
	assert (varN>=0.0)
	return float(varN)

def third_moment_nmut_exact(L,k,p):
    t1 = (L * (-2 + 3*L) * p**2 + 3 * (1 - p)**(2*k) * (2 + (-1 + k - L) * p * (2 + k * p - L * p)) - (1 - p)**k * (6 + p * (-6 + L * (-6 + p + 6 * L * p))))/(p**2)
    t2 = (-2 + 2 * k - L) * (-1 + 2 * k - L) * (2 * k - L) * (-1 + (1 - p)**k)**3
    t3 = (1/(p**3))*(-6 * (-1 + k)**2 * (k - L) * p**3 + 6 * (1 - p)**(3 * k) * (2 + (-2 + 2 * k - L) * p) + (1 - p)**(2 * k) * (-12 + 6 * (2 * k + L) * p + 6 * (4 * k**2 + 2 * (1 + L) - 3 * k * (2 + L)) * p**2 - (-1 + k) * k * (-2 + 4 * k - 3 * L) * p**3) + 6 * (-1 + k) * (1 - p)**k * p * (-2 + p * (2 - k + 2 * L + (k * (-2 + 3 * k - 3 * L) + L) * p)))
    t4 = 6 * (-1 + (1 - p)**k) * ((k + k**2 - 2 * k * L + (-1 + L) * L) * (-1 + 2 * (1 - p)**k) * hyp2f1(1, 2 + k - L, k - L, 1) + (k + k**2 - 2 * k * L + (-1 + L) * L) * (1 - p)**k * (-1 + p) * hyp2f1(1, 2 + k - L, k - L, 1 - p) - (-2 * k + 4 * k**2 + L - 4 * k * L + L**2) * ((-1 + 2 * (1 - p)**k) * hyp2f1(1, 1 + 2 * k - L, -1 + 2 * k - L, 1)- (-1 + p)**(2 * k) * hyp2f1(1, 1 + 2 * k - L, -1 + 2 * k - L, 1 - p)))
    return t1+t2+t3+t4

def exp_n_mutated_squared(L, k, p):
    return var_n_mutated(L, k, p) + exp_n_mutated(L, k, p) ** 2

def exp_n_mutated_cubed(L, k, p):
    return third_moment_nmut_exact(L, k, p)

def exp_n_mutated_to_the_fourth_power(L, k, p):
    return fourth_moment_using_normal(L, k, p)

def var_n_mutated_squared(L, k, p):
    return exp_n_mutated_to_the_fourth_power(L, k, p) - exp_n_mutated_squared(L, k, p) ** 2

def third_moment_nmut_using_normal(L,k,p):
    mu = exp_n_mutated(L, k, p)
    var = var_n_mutated(L, k, p)
    third_moment = mu**3 + 3 * mu * var
    return third_moment

def exp_probability_nothing_common(L, k, p, s):
    c = exp_n_mutated(L, k, p)
    print("exp_n_mut = " + str(c))
    term1 = (1-s)**(L-c)
    term2 = -1.0 * (1-s)**(L-c) * log(1-s) * (exp_n_mutated(L, k, p) - c)
    term3 = 0.5 * (1-s)**(L-c) * ( log(1-s) )**2 * (exp_n_mutated_squared(L, k, p) - c**2)
    term4 = -1.0/6.0 * (1-s)**(L-c) * ( log(1-s) )**3 * (exp_n_mutated_cubed(L, k, p) - 3*c*exp_n_mutated_squared(L, k, p) + 2*c**3)
    ret_val = max(round(sum([term1, term2, term3, term4]),9),0.0)
    print([term1, term2, term3, term4])
    return min(max(ret_val, 0.0), 1.0)

def exp_probability_all_common(L, k, p, s):
    c = exp_n_mutated(L, k, p)
    print("exp_n_mut = " + str(c))
    term1 = (1-s)**(2*c)
    term2 = 1.0 * (1-s)**(2*c) * log(1-s) * (exp_n_mutated(L, k, p) - c)
    term3 = 2 * (1-s)**(2*c) * ( log(1-s) )**2 * (exp_n_mutated_squared(L, k, p) - c**2)
    term4 = 4.0/3.0 * (1-s)**(2*c) * ( log(1-s) )**3 * (exp_n_mutated_cubed(L, k, p) - 3*c*exp_n_mutated_squared(L, k, p) + 2*c**3)
    ret_val = max(round(sum([term1, term2, term3, term4]),9),0.0)
    print([term1, term2, term3, term4])
    return min(max(ret_val, 0.0), 1.0)

def f_test_plot(L, k, p, s):
	c = exp_n_mutated(L, k, p)
	print("exp_n_mut = " + str(c))
	term1 = (1-s)**(L-c)
	term2 = -1.0 * (1-s)**(L-c) * log(1-s) * (exp_n_mutated(L, k, p) - c)
	term3 = 0.5 * (1-s)**(L-c) * ( log(1-s) )**2 * (exp_n_mutated_squared(L, k, p) - c**2)
	term4 = -1.0/6.0 * (1-s)**(L-c) * ( log(1-s) )**3 * (exp_n_mutated_cubed(L, k, p) - 3*c*exp_n_mutated_squared(L, k, p) + 2*c**3)
	#return ( log(1-s) )**3 * (1-s)**(L-c) * (exp_n_mutated_cubed(L, k, p) - 3*c*exp_n_mutated_squared(L, k, p) + 2*c**3)
	return sum([term1, term2, term3, term4])
	#return exp_n_mutated_cubed(L, k, p) + 2 * c**3
	#return exp_n_mutated_cubed(L, k, p) + 2 * c**3 - 3*c*exp_n_mutated_squared(L, k, p)
	#return -1.0/6.0 * (1-s)**(L-c) * ( log(1-s) )**3

def f_test_all_common(L, k, p, s):
	c = exp_n_mutated(L, k, p)
	term1 = (1-s)**(2*c)
	term2 = 1.0 * (1-s)**(2*c) * log(1-s) * (exp_n_mutated(L, k, p) - c)
	term3 = 2 * (1-s)**(2*c) * ( log(1-s) )**2 * (exp_n_mutated_squared(L, k, p) - c**2)
	return sum([term1, term2, term3])

def get_probabilities(L, k, p):
	mu = exp_n_mutated(L, k, p)
	if int(mu) == L:
		return [0]*L + [1]
	sigma = sqrt(var_n_mutated(L, k, p))
	pdf_ = norm(mu, sigma).pdf
	probabilities = [ pdf_(nMut) for nMut in range(0, L+1)]
	sum_ = sum(probabilities)
	return [p/sum_ for p in probabilities]

def f_test_pdf(L, k, p, s):
	probabilities = get_probabilities(L, k, p)
	total = 0.0
	for x in range(0, L+1):
		if probabilities[x] == 0:
			continue
		term = (L-x) * log(1-s) + log(probabilities[x])
		total += exp(term)
	return total

if __name__ == '__main__':
	L = 10000
	k = 21
	p = 0.2
	s = 0.1

	sigma = sqrt(var_n_mutated(L, k, p))
	mu = exp_n_mutated(L, k, p)
	const = sigma * log(1.0-s)/sqrt(2.0)
	lst = []
	print(const, mu, L)
	for i in [2*x for x in range(1,1000)]:
		logarithm = (L-mu) * log(1.0-s) - sum( log(y) for y in range(1,i/2+1) )
		lst.append(exp(logarithm)*const**i)
		if i==2:
			print((L-mu) * log(1.0-s), logarithm)
	print(sum(lst))


	'''
	print('testing nothing common')
	print(exp_probability_nothing_common(L, k, p, s))
	print('testing all common')
	print(exp_probability_all_common(L, k, p, s))
	'''

	# Below, we will do some experimentations on exp_probability_nothing_common
	# 1. Plot the individual terms
	# 2. Plot the raw sum
	# 3. Plot clipped sum
	# 4. Find out why not so sensitive w.r.t. p

	'''
	p_values = np.linspace(0.00001, 1.0, 100)
	# f_test_plot f_test_pdf
	f_values = np.vectorize(f_test_plot)(L, k, p_values, s)
	plt.plot(p_values, f_values)
	plt.xlabel('mutation rate, p')
	plt.ylabel('f_value')
	plt.title('L = %d, k = %d, s = %f' %(L, k, s))
	plt.savefig('f_test_3.pdf')
	'''
