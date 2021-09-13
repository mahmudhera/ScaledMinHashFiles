## Overview
The codes in this directory are used to calculate the expected probability of the cases where there is nothing common between the sketches.

To get nothing common in two sketches, there are $`L-N_{mut}`$ k-mers that are not mutated by the mutation process, and all of these need to be thrown out while taking a sketch. The probability of that happening is $`(1-s)^{L - N_{mut}}`$. Next, using mathematics, we can calculate the expectation of this probability.

The python file implements this expected probability. The pdf files are sample plots.

### plot_expected_probability_of_pathological_cases.python

This is the code used to plot the expected probability. The other files are used to make the whole thing run.

#### Results

The results from running this script are the PDF files in this directory (from different viewpoints). For example, if we look at `combined2.pdf`, we will see three figures for different k-mer sizes. In each of these figures, for a fixed `L = 10^9`, we see the expected probability for all scale-factors and mutation rates. Because this is a 3D plot, we change the viewing direction in the other files.

1. combined.pdf: L = 10^9
1. combined2.pdf: L = 10^9
1. combined3.pdf: L = 10^9
1. combined4.pdf: L = 10^12

These plots reveal the range of scale-factors and mutation-rates for which our confidence intervals are not going to have any trouble. It shows that with increased scale-factor, the probability of finding nothing common in the two sketches decreases. The plots also show that with larger mutation rate, we would need to use a larger scale-factor. We can also see that with increased k-mer size and decreased L, the range where our intervals work well shrinks, as expected.

### check_expected_probability_against_simulations.py
In this script, the code to run the following experiment resides: since we do not know if the formula (for expected probability of having nothing common in the sketches) is correct, we calculate the probability of finding nothing in common in the two sketches (before and after the mutation process) from multiple simulations, and compare the value with what our formula spits out.

#### Results
For 10k simulations, we calculated if nothing common in the sketches. We count the # of times
this happens. We also counted the percentage of times this happens from all the experiments. Then we compared and made a contrast
with the value we have from formula. The results are as follows:

|L|k|scale_factor|mutation_rate|estimated_from_experiments|estimated_from_formula
|---|---|---|---|---|---|
|10000|21|0.1|0.2|0.001|0.00299|
|10000|31|0.1|0.2|0.491|0.302|
|10000|21|0.1|0.3|0.6292|0.434|
|10000|31|0.1|0.3|0.9877|0.520|
|10000|21|0.1|0.4|0.9811|0.544|
|10000|31|0.1|0.4|1.0|0.504|

With these results, we can conclude that this formula has something seriously wrong. We need to figure that out.

##### Problems:
1. Discrepancy in the figures in this table
1. The problem we have is: we have an expression for the expected probability of nothing common in sketches. The formula has some issues. At scale_factor = 0, the probability rises, then falls. Conceptually, it should just be 1.0, as we are making sure that there is nothing common in the sketches.

#### Ideas

1. Try the PDF
1. Consider two random processes at the same time
1. Try area correction
1. Taylor expansion
1. Check literature for PMF of sum of m-dependent Bernouli


### Area correction
Sep 8, 2021

When we approximate a discrete distribution (for example, N_mut) with a continuous distribution (such as the normal distribution), we often have to do area correction. By using that, probability(N_mut = n) = intergration of continuous PDF from n-0.5 to n+0.5. The range of N_mut is [0, L]. Therefore, we scale the PMF by scaling with a factor of the integration from -0.5 to L+0.5.

This is implemented the file `check_expected_probability_against_simulations.py` as follows:

```
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
```

Implementation notes: `f = lambda a, b: 0.5 * ( erf( (mu-a)/(sqrt(2)*sigma) ) - erf( (mu-b)/(sqrt(2)*sigma) ) )` is the integral from a to b.

After implementing this, we calculated the expected probability from multiple simulations, and compared with what we get from the formula. The results are as follows.

|L|k|scale_factor|mutation_rate|estimated_from_experiments|estimated_from_formula
|---|---|---|---|---|---|
|10000|21|0.1|0.2|0.001|0.003|
|10000|31|0.1|0.2|0.491|0.362|
|10000|21|0.1|0.3|0.6292|0.529|
|10000|31|0.1|0.3|0.9877|0.942|
|10000|21|0.1|0.4|0.9811|0.940|
|10000|31|0.1|0.4|1.0|0.999|

This shows that there is still room for improvement. Still, the formula looks way better than the previous, which never got values close to 1.0 (stuck at around 0.5).

### Getting a lower bound

This has been proved by David. Check the notebook PDF on dropbox for the proof. In short, it has been mathematically proved that the expected probability is at least (>=) (1-s)^[L*(1-p)**k]. Implemented this is: `check_expected_probability_against_simulations.py`.

```
def exp_probability_path_case_david(L, k, p, s):
    return (1-s)**(L*(1-p)**k)
```

Indeed, by matching with simulations, we see the following:

|L|k|scale_factor|mutation_rate|estimated_from_experiments|estimated_lower_bound
|---|---|---|---|---|---|
|10000|21|0.1|0.2|0.001|6e-5|
|10000|31|0.1|0.2|0.491|0.352|
|10000|21|0.1|0.3|0.6292|0.555|
|10000|31|0.1|0.3|0.9877|0.983|
|10000|21|0.1|0.4|0.9811|0.977|
|10000|31|0.1|0.4|1.0|0.9998|

At least for this very brief set of experiments, the lower bound seems to be working.