import numpy as np

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

def perform_exhaustive_counting(k, p, L):
    arr = np.zeros((2, L+1))
    # arr[i,j] indicates c[k,i,j]
    #arr[0, k] = 1 # (1-p)**k
    arr[0, k] = (1-p)**k
    for i in range(1,k+1):
        #arr[1, k-i] = 2**(i-1) # p * (1-p)**(k-i)
        arr[1, k-i] = p * (1-p)**(k-i)
    return arr

def get_nmut_pdf(L, k, p):
    # here, L means # of bases, length of the whole string
    arr = perform_exhaustive_counting(k, p, L)
    curr = np.zeros((L+1, L+1))
    next = np.zeros((L+1, L+1))
    curr[0] = arr[0]
    curr[1] = arr[1]

    for l in range(k+1, L+1):
        next[next>0] = 0
        for x in range(L):
            for z in range(L):
                if z < k-1:
                    next[x+1, z+1] += curr[x, z] * (1-p)
                else:
                    next[x, z+1] += curr[x, z] * (1-p)
                next[x+1, 0] += curr[x, z] * p
        next, curr = curr, next
    #print(curr)
    #print(np.sum(curr, axis=0))
    #print(int(sum(sum(curr))))
    #print(2**L)
    return np.sum(curr, axis=1)

L = 50
k = 21
p = 0.1
q = 1 - (1-p)**k
print('Expectation from formula:')
print(L*q)
print('Expectation from pdf:')
pdf = get_nmut_pdf(L+k-1, k, p)
expectation = sum( [i*pdf[i] for i in range(L+1)] )
print(expectation)

print('2nd moment from formula:')
print( exp_n_mutated_squared(L,k,p) )
print('2nd moment from pdf:')
pdf = get_nmut_pdf(L+k-1, k, p)
expectation = sum( [i**2*pdf[i] for i in range(L+1)] )
print(expectation)

'''
print(sum(sum(arr)))
print(2**k)

print(arr)
L = k+1
arr_L = np.zeros((k+1, k+2))
for x in range(L-k+1):
    for z in range(k+1):
        if z < k-1:
            arr_L[x+1, z+1] += arr[x,z]
        else:
            arr_L[x, z+1] += arr[x,z]
        arr_L[x+1, 0] += arr[x,z]
        print(x, z)
print(arr_L)
print(sum(sum(arr_L)))
'''
'''
[[  0.   0.   0.   0.   0.   0.   1.]
 [  1.   0.   0.   0.   0.   1.   0.]
 [ 31.  16.   8.   4.   2.   0.   0.]
 [  0.   0.   0.   0.   0.   0.   0.]
 [  0.   0.   0.   0.   0.   0.   0.]
 [  0.   0.   0.   0.   0.   0.   0.]]
 [ 0.12157665,  0.02701703,  0.03151987,  0.0366898,   0.04261946,  0.18449904, 0.0707005, 0.07409356,  0.07605214,  0.07600764,  0.08826228,  0.05290512, 0.04346563,  0.03291675,  0.02219589,  0.01276791,  0.00671071]
'''
