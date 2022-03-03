from math import log, exp

def r1_to_q(k,r1):
	#return 1-(1-r1)**k
	r1 = float(r1)
	q = 1-(1-r1)**k
	return float(q)

def exp_n_mutated(L,k,r1):
	q = r1_to_q(k,r1)
	return L*q

def get_expected_log_probability(L, k, p, s):
    exp_nmut = exp_n_mutated(L, k, p)
    try:
        return (L - exp_nmut) * log(1.0 - s)
    except:
        return float('-inf')

def get_exp_probability_measure(L, k, p, s):
    return exp( get_expected_log_probability(L, k, p, s) )

if __name__ == '__main__':
    for (L, k, p, s) in [ (10000, 21, 0.2, 0.1), (10000, 31, 0.2, 0.1), (10000, 21, 0.3, 0.1), (10000, 31, 0.3, 0.1) ]:
        print( L, k, p, s, get_exp_probability_measure(L, k, p, s) )
