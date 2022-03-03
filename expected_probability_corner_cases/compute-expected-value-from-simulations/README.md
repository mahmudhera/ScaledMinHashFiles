In this experiment we will simulate point mutations using Binomial, then take random samples, and then return true/false if there is nothing common between the sketches. We will perform this multiple times to generate expected number of times we see nothing common.

Using 2000 runs, the following few values were observed:

(L, k, p, s, exp_probability)

L = 10k
(10000, 21, 0.2, 0.1, 0.001)
(10000, 31, 0.2, 0.1, 0.4915)
(10000, 21, 0.3, 0.1, 0.6275)
(10000, 31, 0.3, 0.1, 0.988)

L = 1k
(1000, 21, 0.2, 0.1, 0.502)
(1000, 31, 0.2, 0.1, 0.9335)
(1000, 21, 0.3, 0.1, 0.96)
(1000, 31, 0.3, 0.1, 0.9965)
