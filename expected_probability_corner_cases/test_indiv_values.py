from expected_probability_corner_cases import f_test_pdf

ks = [21, 31]
ps = [0.2, 0.3, 0.4]
L = 10000
s = 0.1
for p in ps:
	for k in ks:
		print(L, k, p, s, f_test_pdf(L, k, p, s))