from pdf_nmut import get_nmut_pdf

def exp_probability_using_pdf(num_kmers, kmer_len, mut_rate, scale_factor, pdf=None):
    L = num_kmers
    k = kmer_len
    p = mut_rate
    s = scale_factor

    if pdf==None:
        pdf = get_nmut_pdf(L+k-1, k, p)
    expectation = 0.0
    for nmut in range(L+1):
        expectation += (1-s)**(L-nmut) * pdf[nmut]
    return expectation

if __name__ == '__main__':
    mutation_rates = [0.2, 0.3]
    scale_factors = [0.1]
    k_mer_lengths = [21, 31]
    num_k_mers_list = [1000]
    num_simulations = 2000

    for num_k_mers in num_k_mers_list:
        for mutation_rate in mutation_rates:
            for k_mer_length in k_mer_lengths:
                pdf = get_nmut_pdf(num_k_mers+k_mer_length-1, k_mer_length, mutation_rate)
                for scale_factor in scale_factors:
                    estimated_expected_probability = exp_probability_using_pdf(num_k_mers, k_mer_length,
                                                                            mutation_rate, scale_factor, pdf)
                    print (num_k_mers, k_mer_length, mutation_rate, scale_factor, estimated_expected_probability)
