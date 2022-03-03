from pdf_nmut import get_nmut_pdf

def exp_probability_using_pdf(num_kmers, kmer_len, mut_rate, scale_factor):
    L = num_kmers
    k = kmer_len
    p = mut_rate
    s = scale_factor

    pdf = get_nmut_pdf(L+k-1, k, p)
    expectation = 0.0
    for nmut in range(L+1):
        expectation += (1-s)**(L-nmut) * pdf[nmut]
    return expectation

L = 100
k = 13
p = 0.1
s = 0.1

print ( exp_probability_using_pdf(L, k, p, s) )
