import numpy as np
from scipy.special import hyp2f1

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

def exp_n_mutated_squared(L, k, p):
    return var_n_mutated(L, k, p) + exp_n_mutated(L, k, p) ** 2

def third_moment_nmut_exact(L,k,p):
    t1 = (L * (-2 + 3*L) * p**2 + 3 * (1 - p)**(2*k) * (2 + (-1 + k - L) * p * (2 + k * p - L * p)) - (1 - p)**k * (6 + p * (-6 + L * (-6 + p + 6 * L * p))))/(p**2)
    t2 = (-2 + 2 * k - L) * (-1 + 2 * k - L) * (2 * k - L) * (-1 + (1 - p)**k)**3
    t3 = (1/(p**3))*(-6 * (-1 + k)**2 * (k - L) * p**3 + 6 * (1 - p)**(3 * k) * (2 + (-2 + 2 * k - L) * p) + (1 - p)**(2 * k) * (-12 + 6 * (2 * k + L) * p + 6 * (4 * k**2 + 2 * (1 + L) - 3 * k * (2 + L)) * p**2 - (-1 + k) * k * (-2 + 4 * k - 3 * L) * p**3) + 6 * (-1 + k) * (1 - p)**k * p * (-2 + p * (2 - k + 2 * L + (k * (-2 + 3 * k - 3 * L) + L) * p)))
    t4 = 6 * (-1 + (1 - p)**k) * ((k + k**2 - 2 * k * L + (-1 + L) * L) * (-1 + 2 * (1 - p)**k) * hyp2f1(1, 2 + k - L, k - L, 1) + (k + k**2 - 2 * k * L + (-1 + L) * L) * (1 - p)**k * (-1 + p) * hyp2f1(1, 2 + k - L, k - L, 1 - p) - (-2 * k + 4 * k**2 + L - 4 * k * L + L**2) * ((-1 + 2 * (1 - p)**k) * hyp2f1(1, 1 + 2 * k - L, -1 + 2 * k - L, 1)- (-1 + p)**(2 * k) * hyp2f1(1, 1 + 2 * k - L, -1 + 2 * k - L, 1 - p)))
    return t1+t2+t3+t4

def perform_exhaustive_counting(k, p, num_bases):
    arr = np.zeros((2, num_bases+1))
    # arr[i,j] indicates c[k,i,j]
    #arr[0, k] = 1 # (1-p)**k
    arr[0, k] = (1-p)**k
    for i in range(1,k+1):
        #arr[1, k-i] = 2**(i-1) # p * (1-p)**(k-i)
        arr[1, k-i] = p * (1-p)**(k-i)
    return arr

def get_nmut_pdf(num_bases, k, p):
    # here, L means # of bases, length of the whole string
    arr = perform_exhaustive_counting(k, p, num_bases)
    curr = np.zeros((num_bases+1, num_bases+1))
    next = np.zeros((num_bases+1, num_bases+1))
    curr[0] = arr[0]
    curr[1] = arr[1]

    for l in range(k+1, num_bases+1):
        next[next>0] = 0
        for x in range(num_bases):
            for z in range(num_bases):
                if z < k-1:
                    next[x+1, z+1] += curr[x, z] * (1-p)
                else:
                    next[x, z+1] += curr[x, z] * (1-p)
                next[x+1, 0] += curr[x, z] * p
        next, curr = curr, next
        print('done for l = ' + str(l))
    return np.sum(curr, axis=1)
