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

L, k, scale_factor, mutation_rate, estimated_from_experiments, estimated_from_formula
10000, 21, 0.1, 0.2, 0.001, 0.0029871921240774667
10000, 31, 0.1, 0.2, 0.491, 0.30206161823753114
10000, 21, 0.1, 0.3, 0.6292, 0.43434678390200443
10000, 31, 0.1, 0.3, 
10000, 21, 0.1, 0.4, 
10000, 31, 0.1, 0.4, 
